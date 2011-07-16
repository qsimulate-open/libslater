//
// author : Toru Shiozaki
// date   : June 2009
//
#include <vector>
#include "mpreal.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <sstream>
#include <string>
#include "local_macros.h"
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace mpfr;
using namespace boost;

extern void find_root_weight(const mpreal T, const mpreal U, vector<mpreal>& dx, vector<mpreal>& dw, const int rank, const bool);
extern vector<vector<double> > get_C(const int, const int, int);

void generate() {
  const int npoints = NGRID;
  const int tstart = 0;
  const int tfence = 10;
  const int ustart = -7;
  const int ufence = 3;
  const mpreal tiny = "1.0e-100";
  const mpreal zero = "0.0e0";

  const int nboxes = (tfence - tstart) * (ufence - ustart);
  const int npoints_sq = npoints * npoints;

for (int rank = 1; rank != 14; ++rank) {
  ofstream ofs;
  stringstream ofs2;
  stringstream list;
  const string filename = "_root" + lexical_cast<string>(rank) + ".f90";
  ofs.open(filename.c_str());
  assert(ofs.is_open());

  const string indent = "      ";
  const string contline = "     &";
  const string srank = lexical_cast<string>(2 * rank);
  const string snpoints = lexical_cast<string>(npoints);
  const string snboxes = lexical_cast<string>(nboxes);

  ofs <<
"!\n\
! Author : Toru Shiozaki\n\
! Date   : June 2009\n\
!\n\
      subroutine fill" << rank << "\n\
      use rank" << rank << endl;

  ofs2 << 
"!\n\
!\n\
      subroutine root" << rank <<  "(ta, ua, rr, ww, n)\n\
      use rank" << rank << "\n\
      implicit none\n\
      integer i, j, k, n, offset, it, iu, of(" << rank << ")\n\
      double precision t, u, t2, u2, tt, uu, d(" + srank + "), g(" + srank + ")\n\
      double precision ta(*), ua(*), rr(*), ww(*)\n\
      double precision im(" + srank + ", " + snpoints + "), imc(" + srank + ", " + snpoints + ")\n\
      double precision pi, o7, o14\n\
      parameter (pi = 3.141592653589793238462643d0)\n\
      parameter (o7  = (2.0d0 / " + snpoints + ".0d0))\n\
      parameter (o14 = (1.0d0 / " + snpoints + ".0d0))\n";

  int nblock = 0;
  int index = 0;
  for (int tbase = tstart; tbase != tfence; ++tbase) { // sqrt(t)
    for (int ubase = ustart; ubase != ufence; ++ubase, ++nblock) { // log(u) 
      list << "\
!\n\
      subroutine fill" << rank << "_" << nblock <<  "\n\
      use rank" << rank << "\n\
      implicit none\n\
";
      ofs << "\
      call fill" << rank << "_" << nblock << endl;

      const vector<vector<double> > c_all = get_C(tbase, ubase, rank);
      assert(c_all.size() == rank * 2);

      for (int i = 0; i != rank; ++i, ++index) {
        const int ii = 2 * i;
        const vector<double> x = c_all[ii]; 
        const vector<double> w = c_all[ii + 1]; 

        int cnt = 0;
        for (vector<double>::const_iterator iter = x.begin(); iter != x.end(); ++iter, ++cnt)
          list << indent << "x(" << index * npoints_sq + 1 + cnt << ") = " << scientific << setprecision(17) << setw(25) << 
                 (fabs(*iter) < tiny ? 0.0 : *iter) << "\n"; 
        cnt = 0;
        for (vector<double>::const_iterator iter = w.begin(); iter != w.end(); ++iter, ++cnt) {
          list << indent << "w(" << index * npoints_sq + 1 + cnt << ") = " << scientific << setprecision(17) << setw(25) << 
                 (fabs(*iter) < tiny ? 0.0 : *iter) << "\n";
        }
      }
      list << indent << "end\n";
    }
  }
  ofs << indent << "end" << endl;

  ofs << list.str() << endl;


  ofs2<< "\
      do i = 1, n ! loop over parameter set\n\
        t = ta(i) \n\
        if (t .lt. 0.0d0) cycle \n\
        if (t .gt. 19682.99d0) t = 19682.99d0 \n\
        u = ua(i) \n\
        tt = dsqrt(t) \n\
        if (tt .gt. 1.0d0) tt = dlog(t) * 0.9102392266268373d0 + 1.0d0 ! log(3)+1 \n\
        uu = dlog10(u)\n\
\n\
        it = tt \n\
        tt = tt - it\n\
        tt = 2.0d0 * tt - 1.0d0 \n\
\n\
        iu = uu + 7 ! 0 <= iu <= 9\n\
        if (iu < 0 .or. iu > 10) then\n\
          print *, \"current implementation assumes 1.0e-7 < U < 1.0e3\"\n\
          stop  \n\
        else\n\
          uu = uu - (iu - 7)\n\
          uu = 2.0d0 * uu - 1.0d0 \n\
        endif \n\
\n\
        offset = " << (rank * npoints_sq) << " * (iu + it * " + lexical_cast<string>(ufence - ustart) + ")\n\
        u2 = uu * 2.0d0\n\
        do j = 1, " + snpoints + "\n";
  for (int i = 0; i != rank * 2; ++i) { 
    ofs2<< "\
          d(" << (i + 1) << ") = 0.0d0" << endl;
    ofs2<< "\
          g(" << (i + 1) << ") = 0.0d0" << endl;
  }
  ofs2<< "\
          of(1) = offset + (j - 1) * " << npoints << endl;
  for (int i = 2; i <= rank; ++i)
    ofs2<< "\
          of(" << i << ") = of(" << (i - 1) << ") + " << npoints_sq << endl;
  for (int k = npoints; k >= 2; --k) {
    if (k % 2 == 0 ) {
      for (int i = 0; i != rank; ++i) { 
        ofs2<< "\
          g(" << (i + 1) << ") = u2 * d(" << (i + 1) <<  ") - g(" << (i + 1) << ") + x(of(" << (i + 1) << ") + " << k << ")" << endl;
      } 
      for (int i = rank; i != rank * 2; ++i) {
        ofs2<< "\
          g(" << (i + 1) << ") = u2 * d(" << (i + 1) <<  ") - g(" << (i + 1) << ") + w(of(" << (i - rank + 1) << ") + " << k << ")" << endl;
      }
    } else {
      for (int i = 0; i != rank; ++i) { 
        ofs2<< "\
          d(" << (i + 1) << ") = u2 * g(" << (i + 1) <<  ") - d(" << (i + 1) << ") + x(of(" << (i + 1) << ") + " << k << ")" << endl;
      } 
      for (int i = rank; i != rank * 2; ++i) {
        ofs2<< "\
          d(" << (i + 1) << ") = u2 * g(" << (i + 1) <<  ") - d(" << (i + 1) << ") + w(of(" << (i - rank + 1) << ") + " << k << ")" << endl;
      }
    }
  }
  for (int i = 0; i != rank; ++i) { 
    ofs2<< "\
          im(" << (i + 1) << ", j) = uu * g(" << (i + 1) <<  ") - d(" << (i + 1) << ") + x(of(" << (i + 1) << ") + 1) * 0.5d0" << endl;
  } 
  for (int i = rank; i != rank * 2; ++i) {
    ofs2<< "\
          im(" << (i + 1) << ", j) = uu * g(" << (i + 1) <<  ") - d(" << (i + 1) << ") + w(of(" << (i - rank + 1) << ") + 1) * 0.5d0" << endl;
  }
  ofs2<< "\
        enddo" << endl;
  for (int j = 1; j <= npoints; ++j) {
    for (int i = 1; i <= rank * 2; ++i) 
      ofs2<< "\
        d(" << i << ") = 0.0d0" << endl;

      mpreal pij = PI;
      pij *= (j - 1);
      const mpreal mp28 = npoints * 2;
      for (int k = 1; k <= npoints; ++k) {
        for (int i = 1; i <= rank * 2; ++i) {
          const double mpcos = static_cast<double>(mpfr::cos(pij * (k + k - 1) / mp28));
          if (fabs(mpcos) < tiny) {
          ofs2<< "\
!       d(" << i << ") = d(" << i <<  ")" << endl;
          } else if (fabs(mpcos-1.0) < tiny) {
          ofs2<< "\
        d(" << i << ") = d(" << i <<  ") + im(" << i << ", " << k << ")" << endl;
          } else if (fabs(mpcos+1.0) < tiny) {
          ofs2<< "\
        d(" << i << ") = d(" << i <<  ") - im(" << i << ", " << k << ")" << endl;
          } else {
          ofs2<< "\
        d(" << i << ") = d(" << i <<  ") + im(" << i << ", " << k << ") * ("
                << setw(25) << setprecision(17) << scientific <<  mpcos << ")" << endl;
          }
       }
    }
    for (int i = 1; i <= rank * 2; ++i) {
      ofs2<< "\
        imc(" << i << ", " << j << ") = o7 * d(" << i << ")" << endl;
    }
  }
  ofs2<< endl << "\
        t2 = tt * 2.0d0" << endl;
  for (int i = 1; i <= rank * 2; ++i)
    ofs2<< "\
        d(" << i << ") = 0.0d0" << endl;
  for (int i = 1; i <= rank * 2; ++i) 
    ofs2<< "\
        g(" << i << ") = 0.0d0" << endl;

  for (int j = npoints; j >= 2; --j) {
    if (j % 2 == 0) {
      for (int i = 1; i <= rank * 2; ++i)
        ofs2<< "\
        g(" << i << ") = t2 * d(" << i << ") - g(" << i << ") + imc(" << i << ", " << j << ")" << endl;
    } else {
      for (int i = 1; i <= rank * 2; ++i)
        ofs2<< "\
        d(" << i << ") = t2 * g(" << i << ") - d(" << i << ") + imc(" << i << ", " << j << ")" << endl;
    }
  }
  ofs2<< endl;
  ofs2<< "\
        offset = (i - 1) * " << rank << endl;
  ofs2<< "\
        uu = 1.0d0 / dsqrt(u)" << endl;
  for (int i = 1; i <= rank; ++i)
    ofs2<< "\
        rr(offset + " << i << ") = tt * g(" << i << ") - d(" << i << ") + imc(" << i << ", 1) * 0.5d0" << endl;
  for (int i = 1; i <= rank; ++i)
    ofs2<< "\
        ww(offset + " << i << ") = (tt * g(" << (i + rank) << ") - d(" << (i + rank) << ") + imc(" << (i + rank) << ", 1) * 0.5d0) * uu" << endl;

  ofs2<< "\
      enddo" << endl << endl;
  ofs2<< "\
      return\n\
      end" << endl;

  ofs << ofs2.str() << endl;
  ofs.close();
}

}
