// 	$Id: Image.cpp,v 1.30 2011/02/09 19:22:23 garyb Exp $	
// Implementation code for Image class and related classes.
#include "Image.h"
#include <cstring>
#include <cstdio>
#include <cctype>

using namespace img;
using namespace std;

// Case-raising & blank-stripping function for keywords.
string img::KeyFormat(const string input) {
  string output;
  //strip leading blanks
  int i=0;
  while (i<input.size() && isspace(input[i])) i++;
  //Convert to upper case, stop at whitespace
  while (i<input.size() && !isspace(input[i]))
    output += toupper(input[i++]);
  return output;
}

string
HdrRecordBase::writeCard() const {
  string vv=getValueString();
  string card=keyword;
  for (int i=card.size(); i<8; i++) card += " ";
  if (!vv.empty()) {
    card += "= ";
    for (int i=vv.size(); i<20; i++) card += " ";
    card += vv;
    if (!comment.empty() || !units.empty())
      card += " /";
  }
  card += " ";
  if (!units.empty()) card+= "[" + units + "] ";
  card += comment;
  return card;
}

namespace img {
  //specializations for bool
  template<>
  string
  HdrRecord<bool>::getValueString() const {
    if (val) return "T";
    else return "F";
  }

  template<>
  bool 
  HdrRecord<bool>::setValueString(const string _v) {
      istringstream iss(_v.c_str());
      string s;
      string leftover;
      // Should be only one word, T or F:
      if (!(iss >> s) || (iss >> leftover)) return true;
      if (s=="T" || s=="t" || s=="true") val=true;
      else if (s=="F" || s=="f" || s=="false") val=false;
      else return true;
      return false;
  };

  //and for string - enclose in quotes
  template<>
  string 
  HdrRecord<string>::getValueString() const {
    return "'" + val + "'";
  }

  //and for double - use capital E, many digits...
  template<>
  string 
  HdrRecord<double>::getValueString() const {
    ostringstream oss;
    oss << right << setw(20) << uppercase << showpos
	<< setprecision(12) << val;
      return oss.str();
  }
}
/////////////////////////////////////////////////////////////////////
// Constructor for out-of-bounds that has coordinate info
/////////////////////////////////////////////////////////////////////
void
ImageBounds::FormatAndThrow(const string &m, 
			    const int min, 
			    const int max, 
			    const int tried) {
  char sbuff[120];
  sprintf(sbuff,"Attempt to access %8s number %5d, range is %5d to %5d", 
	  m.c_str(), tried, min, max);
  string mm(sbuff);
  throw ImageError(mm);
}

void
ImageBounds::FormatAndThrow(const int x, 
			    const int y, 
			    const Bounds<int> b) {
  char sbuff[120];
  if (x<b.getXMin() || x>b.getXMax()) {
    sprintf(sbuff,"Attempt to access column number %5d, range is %5d to %5d", 
	  x, b.getXMin(), b.getXMax());
  } else if (y<b.getYMin() || y>b.getYMax()) {
    sprintf(sbuff,"Attempt to access row number %5d, range is %5d to %5d", 
	  y, b.getYMin(), b.getYMax());
  } else 
    strcpy(sbuff,"Cannot find bounds violation ???");
  string mm(sbuff);
  throw ImageError(mm);
}

/////////////////////////////////////////////////////////////////////
// Routines for the underlying data structure ImageData
/////////////////////////////////////////////////////////////////////

// build or destroy the data array
template <class T>
void ImageData<T>::acquireArrays(Bounds<int> inBounds) {
  bounds = inBounds;

  int	xsize, ysize;
  T *dataArray, *dptr, **rPtrs;
  xsize = bounds.getXMax() - bounds.getXMin() + 1;
  ysize = bounds.getYMax() - bounds.getYMin() + 1;
  
  dataArray  = new T[xsize*ysize]; 

  rPtrs = new T*[ysize]; //??? also be sure to delete dataArray if this fails
  rowPointers = rPtrs - bounds.getYMin();
  dptr = dataArray - bounds.getXMin();
  for (int i=bounds.getYMin(); i<=bounds.getYMax(); i++) {
    rowPointers[i] = dptr;
    dptr += xsize;
  }

  ownDataArray = true;
  ownRowPointers = true;
  isContiguous = true;
}

template <class T>
void
ImageData<T>::discardArrays() {
  if (ownDataArray) {
    //Free the data array(s)
    if (isContiguous) {
      delete [] location(bounds.getXMin(),bounds.getYMin());
    } else {
      T *dptr;
      for (int i=bounds.getYMin(); i<=bounds.getYMax(); i++) {
	dptr = rowPointers[i]+bounds.getXMin();
	delete [] dptr;
      }
    }
  }
  if (ownRowPointers) {
    // Free the row pointer array:
    T **rptr;
    rptr = rowPointers + bounds.getYMin();
    delete [] rptr;
  }
}

// image with unspecified data values:
template <class T>
ImageData<T>::ImageData(const Bounds<int> inBounds): parent(0) {
  acquireArrays(inBounds);
}

// image filled with a scalar:
template <class T>
ImageData<T>::ImageData(const Bounds<int> inBounds, 
			const T initValue): parent(0) {
  acquireArrays(inBounds);
  // Initialize the data
  Assert(isContiguous);
  long int i=0;
  long npix = bounds.getXMax() - bounds.getXMin() + 1;
  npix *= bounds.getYMax() - bounds.getYMin() + 1;

  for (T *initptr=location(bounds.getXMin(), bounds.getYMin()); 
       i<npix;
       i++, initptr++)
    *initptr = initValue;
}

// image for which the data array has been set up by someone else:
template <class T>
ImageData<T>::ImageData(const Bounds<int> inBounds, 
			T** rptrs,
			bool _contig): bounds(inBounds),
				       parent(0),
				       ownDataArray(false),
				       ownRowPointers(true),
				       rowPointers(rptrs),
				       isContiguous(_contig) {}

// image which will be a subimage of a parent:
// Note that there is no assumption that parent is contiguous, its
// storage state could change.
template <class T>
ImageData<T>::ImageData(const Bounds<int> inBounds, 
			const ImageData<T>* _parent): 
  bounds(inBounds),
  parent(_parent),
  ownDataArray(false),
  ownRowPointers(false),
  rowPointers(_parent->rowPointers),
  isContiguous(false) {}

// Destructor:  Be sure to free all memory
template <class T>
ImageData<T>::~ImageData() {
  if (!children.empty() ) 
    throw ImageError("Destroying ImageData that still has children");
  if (parent) parent->unlinkChild(this);
  discardArrays();
}

template <class T>
void
ImageData<T>::unlinkChild(const ImageData<T>* child) const {
  typename list<const ImageData<T>*>::iterator cptr =
    find(children.begin(), children.end(), child);
  if (cptr==children.end() && !uncaught_exception()) 
    throw ImageError("ImageData::unlinkChild cannot find the child");
  children.erase(cptr);
}

template <class T>
void
ImageData<T>::linkChild(const ImageData<T>* child) const {
  children.push_back(child);
}

template <class T>
void
ImageData<T>::deleteSubimages() const {
  while (!children.empty()) {
    delete *(children.begin());
    children.erase(children.begin());
  }
}

// Get a subimage of this one
template <class T>
ImageData<T>* ImageData<T>::subimage(const Bounds<int> bsub) {
  if (!bounds.includes(bsub)) 
    throw ImageError("Attempt to create subimage outside of ImageData"
		     " bounds");
  ImageData<T>* child = new ImageData<T>(bsub, this);
  linkChild(child);
  return child;
}

// Get a read-only subimage of this one
template <class T>
const ImageData<T>* ImageData<T>::const_subimage(const Bounds<int> bsub) const {
  if (!bounds.includes(bsub)) 
    throw ImageError("Attempt to create subimage outside of ImageData"
		     " bounds");
  const ImageData* child = new ImageData(bsub, this);
  linkChild(child);
  return child;
}
  
// Get a read-only subimage of this one (if this is const)
template <class T>
const ImageData<T>* ImageData<T>::subimage(const Bounds<int> bsub) const {
  return const_subimage(bsub);
}

// Create a new subimage that is duplicate of this ones data
template <class T>
ImageData<T>* 
ImageData<T>::duplicate() const {
  ImageData<T>* dup = new ImageData<T>(bounds);
  Assert(dup->isContiguous);

  // Copy the data from old array to new array
  T  *inptr, *dptr;
  long int xsize,ysize;
  xsize = bounds.getXMax() - bounds.getXMin() + 1;
  ysize = bounds.getYMax() - bounds.getYMin() + 1;
  
  if (isContiguous) {
    // Can be done is a single large copy if both contiguous
    inptr = location(bounds.getXMin(),bounds.getYMin());
    dptr = dup->location(bounds.getXMin(),bounds.getYMin());
    memmove( (void *)dptr, (void *)inptr, sizeof(T)*xsize*ysize);
  } else {
    // Do copy by rows
    for (int i=bounds.getYMin(); i<=bounds.getYMax(); i++) {
      inptr = location(bounds.getXMin(), i);
      dptr = dup->location(bounds.getXMin(), i);
      memmove( (void *)dptr, (void *)inptr, sizeof(T)*xsize);
    }
  }
  return dup;
}

// Change image size (flushes data), if data are under this object's
// control and there are no subimages to invalidate.
template <class T>
void
ImageData<T>::resize(const Bounds<int> newBounds) {
  if (newBounds==bounds) return;
  if (!ownDataArray || !ownRowPointers)
    throw ImageError("Cannot ImageData::resize() when data are"
		     " not owned by object");
  if (!children.empty())
    throw ImageError("Attempt to ImageData::resize() with subimages"
		     " in use.");
  discardArrays();
  acquireArrays(newBounds);
}

// Make this ImageData be a copy of another
template <class T>
void
ImageData<T>::copyFrom(const ImageData<T>& rhs) {
  resize(rhs.getBounds());	//???Exception thrown will have misleading msg
  // Copy the data from old array to new array
  T  *dptr;
  const T* inptr;
  long int xsize = bounds.getXMax() - bounds.getXMin() + 1;
  // Do copy by rows
  for (int i=bounds.getYMin(); i<=bounds.getYMax(); i++) {
    dptr = location(bounds.getXMin(), i);
    inptr = rhs.location(bounds.getXMin(), i);
    memmove( (void *)dptr, (void *)inptr, sizeof(T)*xsize);
  }
}

// Change origin of array, saving data
template <class T>
void
ImageData<T>::shift(int x0, int y0) {
  if (!ownDataArray || !ownRowPointers)
    throw ImageError("Cannot ImageData::shift() when data are"
		     " not owned by object");
  if (!children.empty())
    throw ImageError("Attempt to ImageData::shift() with subimages"
		     " in use.");
  int dx = x0 - bounds.getXMin();
  int dy = y0 - bounds.getYMin();
  for (int i=bounds.getYMin(); i<=bounds.getYMax(); i++)
    rowPointers[i] -= dx;
  rowPointers -= dy;

  bounds.setXMin(bounds.getXMin() + dx);
  bounds.setXMax(bounds.getXMax() + dx);
  bounds.setYMin(bounds.getYMin() + dy);
  bounds.setYMax(bounds.getYMax() + dy);
}
  
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// Image<> class operations
/////////////////////////////////////////////////////////////////////////////

// Create subimage of given image.  
template <class T>
Image<T>
Image<T>::subimage(const Bounds<int> bsub) {
  // Make a subimage of current pixel array:
  ImageData<T>* subD=D->subimage(bsub);
  // New image gets this subimage along with the full header
  // of the parent.  ??? do something to keep subimage from writing
  // to header???
  return Image(subD,H, new int(0), hcount);
}

// Create subimage of given image, this time const.  
template <class T>
const Image<T>
Image<T>::subimage(const Bounds<int> bsub) const {
  // Make a subimage of current pixel array, return :
  ImageData<T>* subD=D->subimage(bsub);
  // New image gets this subimage along with the full header
  // of the parent.  ??? do something to keep subimage from writing
  // to header???
  return Image(subD,H, new int(0), hcount);
}
  
// Make a fresh copy of the image (with fresh header)
template <class T>
Image<T>
Image<T>::duplicate() const {
  ImageData<T>* dD=D->duplicate();
  ImageHeader* dH=H->duplicate();
  return Image(dD, dH);
}

// Image arithmetic operations
template <class T>
void
Image<T>::operator+=(const Image<T> rhs) {
  if (!rhs.getBounds().includes(getBounds())) 
    throw ImageError("+= with smaller image");
  transform_pixel(*this, rhs, plus<T>());
}
template <class T>
void
Image<T>::operator-=(const Image<T> rhs) {
  if (!rhs.getBounds().includes(getBounds())) 
    throw ImageError("-= with smaller image");
  transform_pixel(*this, rhs, minus<T>());
}
template <class T>
void
Image<T>::operator*=(const Image<T> rhs) {
  if (!rhs.getBounds().includes(getBounds())) 
    throw ImageError("*= with smaller image");
  transform_pixel(*this, rhs, multiplies<T>());
}
template <class T>
void
Image<T>::operator/=(const Image<T> rhs) {
  if (!rhs.getBounds().includes(getBounds())) 
    throw ImageError("/= with smaller image");
  transform_pixel(*this, rhs, divides<T>());
}

template <class T>
Image<T>
Image<T>::operator+(const Image<T> rhs) const {
  Bounds<int> b=getBounds()&rhs.getBounds();
  if (!b) throw ImageError("no intersection for binary operation");
  Image<T> result(b);
  transform_pixel(result, *this, rhs, plus<T>());
  return result;
}
template <class T>
Image<T>
Image<T>::operator-(const Image<T> rhs) const {
  Bounds<int> b=getBounds()&rhs.getBounds();
  if (!b) throw ImageError("no intersection for binary operation");
  Image<T> result(b);
  transform_pixel(result, *this, rhs, minus<T>());
  return result;
}
template <class T>
Image<T>
Image<T>::operator*(const Image<T> rhs) const {
  Bounds<int> b=getBounds()&rhs.getBounds();
  if (!b) throw ImageError("no intersection for binary operation");
  Image<T> result(b);
  transform_pixel(result, *this, rhs, multiplies<T>());
  return result;
}
template <class T>
Image<T>
Image<T>::operator/(const Image<T> rhs) const {
  Bounds<int> b=getBounds()&rhs.getBounds();
  if (!b) throw ImageError("no intersection for binary operation");
  Image<T> result(b);
  transform_pixel(result, *this, rhs, divides<T>());
  return result;
}

// instantiate for expected types
template class Image<double>;
template class Image<float>;
template class Image<int>;
template class Image<short>;
template class ImageData<double>;
template class ImageData<float>;
template class ImageData<int>;
template class ImageData<short>;
