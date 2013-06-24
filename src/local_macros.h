#ifndef __SRC_LOCAL_MACROS_H
#define __SRC_LOCAL_MACROS_H

#include "mpreal.h"

constexpr int PREC=2048;
constexpr int NGRID=14;
const static mpfr::mpreal PI= atan(mpfr::mpreal("1", PREC)) * mpfr::mpreal("4", PREC);

#endif
