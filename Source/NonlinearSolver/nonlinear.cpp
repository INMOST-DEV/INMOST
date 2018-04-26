#include "inmost_nonlinear.h"

#if defined(_WIN32) || defined(_WIN64) || defined(WIN32) || defined(WIN64)
__declspec( dllexport ) void nonlinear_stub(){} //to avoid LNK4221 warning
#endif
