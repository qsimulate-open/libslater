//
// BAGEL - Parallel electron correlation program.
// Filename: srootlist.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_RYSINT_SROOTLIST_H
#define __SRC_RYSINT_SROOTLIST_H

#include <functional>
#include <slatermem.h>

namespace bagel {

struct SRootList  {
  private:
    std::function<void (const double*, const double*, double*, double*, const int)> srfunc[14];

    static void root1(const double*, const double*, double*, double*, const int);
    static void root2(const double*, const double*, double*, double*, const int);
    static void root3(const double*, const double*, double*, double*, const int);
    static void root4(const double*, const double*, double*, double*, const int);
    static void root5(const double*, const double*, double*, double*, const int);
    static void root6(const double*, const double*, double*, double*, const int);
    static void root7(const double*, const double*, double*, double*, const int);
    static void root8(const double*, const double*, double*, double*, const int);
    static void root9(const double*, const double*, double*, double*, const int);
    static void root10(const double*, const double*, double*, double*, const int);
    static void root11(const double*, const double*, double*, double*, const int);
    static void root12(const double*, const double*, double*, double*, const int);
    static void root13(const double*, const double*, double*, double*, const int);

  public:
    SRootList() {
      srfunc[1] = &root1;
      srfunc[2] = &root2;
      srfunc[3] = &root3;
      srfunc[4] = &root4;
      srfunc[5] = &root5;
      srfunc[6] = &root6;
      srfunc[7] = &root7;
      srfunc[8] = &root8;
      srfunc[9] = &root9;
      srfunc[10] = &root10;
      srfunc[11] = &root11;
      srfunc[12] = &root12;
      srfunc[13] = &root13;
    }

    void srootfunc_call(const unsigned int i, const double* a0, const double* a1, double* a2, double* a3, const int a4) const {
      return (srfunc[i])(a0, a1, a2, a3, a4);
    }

};

}

#endif

