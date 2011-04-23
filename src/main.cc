//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include "mpreal.h"
#include "local_macros.h"

extern void check(const mpfr::mpreal, const mpfr::mpreal);
extern void generate();

int main() {
  mpfr::mpreal::set_default_prec(PREC);

  generate();
//check("200","2e-5");
}

