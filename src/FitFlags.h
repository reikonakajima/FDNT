// $Id: FitFlags.h,v 1.4 2011/08/25 18:09:42 dgru Exp $

// Separate include file that defines flags to use for different fits
#ifndef FITFLAGS_H
#define FITFLAGS_H

namespace laguerre {
    const int Edge =            1;	// 1: mask goes over edge of image
    const int BDegeneracy=      1 << 1;	// 2: b Solution has low singular values
    const int MaskMismatch=     1 << 2; // 4: Mask poorly fits solution
    const int CentroidMismatch= 1 << 3; // 8: Max likelihood center too far from phase center or mask
    const int SizeOutOfRange=   1 << 4; // 16: Substantial likelihood of galaxy size out of bounds.
    const int SizeFixed=        1 << 5; // 32: Galaxy size was frozen by algorithms
    const int UnderSampled=     1 << 6; // 64: e prob distribution sampled with few points
    const int DidNotConverge=   1 << 7; // 128: past maximum iterations
    const int Singularity=      1 << 8;	// 256: Singular matrix or other failure
    const int OutOfBounds=      1 << 9;	// 512: Fitting region off image
    const int TooLarge=         1 << 10;// 1024: Object size too large
    const int TooElliptical=    1 << 11;// 2048: Solver went too elliptical
    const int GLFailure=        1 << 12;// 4096: Fatal error during GL fitting step of larger process
}
#endif
