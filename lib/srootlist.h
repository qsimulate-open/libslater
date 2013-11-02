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
#include <memory>
#include <vector>

namespace bagel {

struct SRootList  {
  private:
    std::function<void (const double*, const double*, double*, double*, const int)> srfunc[14];

    void root1(const double*, const double*, double*, double*, const int);
    void root2(const double*, const double*, double*, double*, const int);
    void root3(const double*, const double*, double*, double*, const int);
    void root4(const double*, const double*, double*, double*, const int);
    void root5(const double*, const double*, double*, double*, const int);
    void root6(const double*, const double*, double*, double*, const int);
    void root7(const double*, const double*, double*, double*, const int);
    void root8(const double*, const double*, double*, double*, const int);
    void root9(const double*, const double*, double*, double*, const int);
    void root10(const double*, const double*, double*, double*, const int);
    void root11(const double*, const double*, double*, double*, const int);
    void root12(const double*, const double*, double*, double*, const int);
    void root13(const double*, const double*, double*, double*, const int);

    struct SlaterMem {
      private:
        std::vector<std::unique_ptr<double[]>> datax_; 
        std::vector<std::unique_ptr<double[]>> dataw_; 
        void fill1(); void fill2(); void fill3(); void fill4(); void fill5();
        void fill6(); void fill7(); void fill8(); void fill9(); void fill10();
        void fill11(); void fill12(); void fill13();
    
      public:
        SlaterMem() {
          for (int i = 1; i != 14; ++i) {
            datax_.push_back(std::unique_ptr<double[]>(new double[19600*i]));
            dataw_.push_back(std::unique_ptr<double[]>(new double[19600*i]));
          }
          fill1(); fill2(); fill3(); fill4(); fill5(); fill6(); fill7();
          fill8(); fill9(); fill10(); fill11(); fill12(); fill13();
        }
        const double* x(const int i) const { return datax_[i-1].get(); }
        const double* w(const int i) const { return dataw_[i-1].get(); }
    } slatermem__;

  public:
    SRootList() : slatermem__() {
      srfunc[1] = std::bind(&SRootList::root1, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
      srfunc[2] = std::bind(&SRootList::root2, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
      srfunc[3] = std::bind(&SRootList::root3, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
      srfunc[4] = std::bind(&SRootList::root4, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
      srfunc[5] = std::bind(&SRootList::root5, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
      srfunc[6] = std::bind(&SRootList::root6, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
      srfunc[7] = std::bind(&SRootList::root7, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
      srfunc[8] = std::bind(&SRootList::root8, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
      srfunc[9] = std::bind(&SRootList::root9, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
      srfunc[10] = std::bind(&SRootList::root10, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
      srfunc[11] = std::bind(&SRootList::root11, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
      srfunc[12] = std::bind(&SRootList::root12, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
      srfunc[13] = std::bind(&SRootList::root13, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
    }

    void srootfunc_call(const unsigned int i, const double* a0, const double* a1, double* a2, double* a3, const int a4) const {
      return (srfunc[i])(a0, a1, a2, a3, a4);
    }

};

}

#endif

