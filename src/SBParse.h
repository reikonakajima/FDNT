// Routines to parse a string into an SBProfile

#ifndef SBPARSE_H
#define SBPARSE_H
#include "SBProfile.h"

namespace sbp {
  extern SBProfile* SBParse(string in_text, Shear& out_shear);

  // overload function, for when we don't need to specify out_shear
  inline SBProfile* SBParse(string in_text)
  {
    Shear s;         // temporary variable that holds the output shear
    return SBParse(in_text, s);
  }
}

#endif // SBPARSE_H
