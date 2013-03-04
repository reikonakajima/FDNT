// Read header information from ASCII files formatted
// to like like FITS cards (as written by Bertin's software).
// Note this is not too sophisticated about parsing the cards.
// Reads from the stream until end or until an END keyword shows up.
#ifndef HDRFROMSTREAM_H
#define HDRFROMSTREAM_H

#include "Image.h"
namespace img {
  ImageHeader HeaderFromStream(std::istream& is);
}
  // This one writes all header cards to stream:
ostream& operator<<(ostream& os, const img::ImageHeader& h);
#endif  // HDRFROMSTREAM_H
