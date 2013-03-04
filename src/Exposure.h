#ifndef EXPOSURE_H
#define EXPOSURE_H
#include "Std.h"
#include "Astrometry.h"
#include "Match.h"
#include <map>
#include "NameIndex.h"
#include "PixelMapCollection.h"

class Exposure;

class Field {
  // *** careful, the orient is not destroyed here...
public:
  string name;
  Orientation* orient;		// defines the TangentPlane coords used in this field
  std::map<int, MCat> affinities;	// All MCats for this field
  vector<int> exposures;	// Which exposures are in this field
};

class Instrument {
public:
  Instrument(string name_):
    name(name_), nExtensions(0) {}
  string name;
  int nExtensions;
  NameIndex extensionNames;	// Names of all extensions that exist for this instrument
  vector<int> exposures;	// Which exposures use this Instrument
  vector<PixelMapChain> maps;	// Instrument parts of PixelMaps for each extension - key chains
  vector<PixelMap*> pixelMaps;	// Instrument parts of PixelMaps for each extension - callable map
  // Keep track of range of pixel coords per extension
  vector<Bounds<double> > domains;	// Rectangles bounding pixel coords of objects
  void addExtension(string extName) {
    extensionNames.append(extName);
    maps.push_back(PixelMapChain());
    pixelMaps.push_back(0);
    domains.push_back(Bounds<double>());
    nExtensions = extensionNames.size();
    Assert(maps.size()==nExtensions);
    Assert(domains.size()==nExtensions);
  }
  void addDetection(const Detection& d) {
    if (d.extension<0 || d.extension>=nExtensions) {
	cerr << "Instrument " << name
	     << " does not have extension number " << d.extension
	     << endl;
	exit(1);
    }
    domains[d.extension] += Position<double>(d.xpix, d.ypix);
  }
  // Everything cleaned up in destructor:
  ~Instrument() {}
private:
  // Hide:
  Instrument(const Instrument& rhs) {}
  void operator=(const Instrument& rhs) {}
};

class Exposure {
public:
  Exposure(): orient(0), reproject(-1) {}
  string name;
  int field;
  int instrument;
  Orientation *orient;	// Telescope pointing for this one
  PixelMapKey reproject; // Map from exposure's TangentPlane to field's TangentPlane
  PixelMapChain warp;		// Exposure portion of PixelMap
  vector<SubMap*> extensionMaps; // Compounded maps for each extension (owned by PixelMapCollection)
  vector<PixelMap*> startpm;  // Input PixelMap for this extension (owned by this class)
  // Everything cleaned up in destructor:
  ~Exposure() {
    if (orient) delete orient;
    for (int i=0; i<startpm.size(); i++)
      if (startpm[i]) delete startpm[i];
  }
private:
  // Hide:
  Exposure(const Exposure& rhs) {}
  void operator=(const Exposure& rhs) {}
};

#endif
