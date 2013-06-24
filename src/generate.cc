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
  const string filename = "_root" + lexical_cast<string>(rank) + ".cc";
  ofs.open(filename.c_str());
  assert(ofs.is_open());

  const string indent = "  ";
  const string srank = lexical_cast<string>(2 * rank);
  const string snpoints = lexical_cast<string>(npoints);
  const string snboxes = lexical_cast<string>(nboxes);

  ofs << "\
//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: " + filename + "\n\
// Copyright (C) 2013 Toru Shiozaki\n\
//\n\
// Author: Toru Shiozaki <shiozaki@northwestern.edu>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the BAGEL package.\n\
//\n\
// The BAGEL package is free software; you can redistribute it and/or modify\n\
// it under the terms of the GNU Library General Public License as published by\n\
// the Free Software Foundation; either version 2, or (at your option)\n\
// any later version.\n\
//\n\
// The BAGEL package is distributed in the hope that it will be useful,\n\
// but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
// GNU Library General Public License for more details.\n\
//\n\
// You should have received a copy of the GNU Library General Public License\n\
// along with the BAGEL package; see COPYING.  If not, write to\n\
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.\n\
//\n\
\n\
#include <cmath>\n\
#include <stdexcept>\n\
#include <slatermem.h>\n\
#include <srootlist.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
void SlaterMem::fill" << rank << "() {\n\
  double* x = datax_[" << rank-1 << "].get();\n\
  double* w = dataw_[" << rank-1 << "].get();" << endl;

  ofs2 << 
"\n\
\n\
void SRootList::root" << rank <<  "(const double* ta, const double* ua, double* rr, double* ww, const int n) {\n\
  const double* x = slatermem__.x(" << rank << ");\n\
  const double* w = slatermem__.w(" << rank << ");\n\
  int of[" << rank << "];\n\
  double d[" + srank + "], g[" + srank + "];\n\
  double im[" + snpoints + "][" + srank + "], imc[" + snpoints + "][" + srank + "];\n\
  constexpr double pi = 3.141592653589793238462643;\n\
  constexpr double o7 = 2.0 / " + snpoints + ";\n\
  constexpr double o14 = 0.5 * o7;" << endl;

  int nblock = 0;
  int index = 0;
  for (int tbase = tstart; tbase != tfence; ++tbase) { // sqrt(t)
    for (int ubase = ustart; ubase != ufence; ++ubase, ++nblock) { // log(u) 
      const vector<vector<double> > c_all = get_C(tbase, ubase, rank);
      assert(c_all.size() == rank * 2);

      for (int i = 0; i != rank; ++i, ++index) {
        const int ii = 2 * i;
        const vector<double> x = c_all[ii]; 
        const vector<double> w = c_all[ii + 1]; 

        int cnt = 0;
        for (vector<double>::const_iterator iter = x.begin(); iter != x.end(); ++iter, ++cnt)
          ofs << indent << "x[" << index * npoints_sq + cnt << "] = " << scientific << setprecision(17) << setw(25) << 
                 (fabs(*iter) < tiny ? 0.0 : *iter) << ";\n"; 
        cnt = 0;
        for (vector<double>::const_iterator iter = w.begin(); iter != w.end(); ++iter, ++cnt) {
          list << indent << "w[" << index * npoints_sq + cnt << "] = " << scientific << setprecision(17) << setw(25) << 
                 (fabs(*iter) < tiny ? 0.0 : *iter) << ";\n";
        }
      }
    }
  }

  ofs << list.str() << endl;
  ofs << "}" << endl;


  ofs2<< "\
  for (int i = 0; i != n; ++i) {// loop over parameter set\n\
    double t = ta[i];\n\
    if (t < 0.0) continue;\n\
    if (t > 19682.99) t = 19682.99;\n\
    const double u = ua[i];\n\
    double tt = sqrt(t);\n\
    if (tt > 1.0)\n\
      tt = log(t) * 0.9102392266268373 + 1.0; // log(3)+1 \n\
    double uu = log10(u);\n\
\n\
    const int it = tt;\n\
    tt = tt - it;\n\
    tt = 2.0 * tt - 1.0;\n\
\n\
    const int iu = uu + 7; // 0 <= iu <= 9\n\
    if (iu < 0 || iu > 10) {\n\
      throw runtime_error(\"current implementation assumes 1.0e-7 < U < 1.0e3\");\n\
    } else {\n\
      uu = uu - (iu - 7);\n\
      uu = 2.0 * uu - 1.0;\n\
    }\n\
\n\
    int offset = " << (rank * npoints_sq) << " * (iu + it * " + lexical_cast<string>(ufence - ustart) + ");\n\
    const double u2 = uu * 2.0;\n\
    for (int j = 0; j != " + snpoints + "; ++j) {\n";
  for (int i = 0; i != rank * 2; ++i) { 
    ofs2<< "\
      d[" << i << "] = 0.0;" << endl;
    ofs2<< "\
      g[" << i << "] = 0.0;" << endl;
  }
  ofs2<< "\
      of[0] = offset + j * " << npoints << ";" << endl;
  for (int i = 2; i <= rank; ++i)
    ofs2<< "\
      of[" << i-1 << "] = of[" << (i - 2) << "] + " << npoints_sq << ";" << endl;
  for (int k = npoints; k >= 2; --k) {
    if (k % 2 == 0 ) {
      for (int i = 0; i != rank; ++i) { 
        ofs2<< "\
      g[" << i << "] = u2 * d[" << i << "] - g[" << i << "] + x[of[" << i << "] + " << k-1 << "];" << endl;
      } 
      for (int i = rank; i != rank * 2; ++i) {
        ofs2<< "\
      g[" << i << "] = u2 * d[" << i << "] - g[" << i << "] + w[of[" << i - rank << "] + " << k-1 << "];" << endl;
      }
    } else {
      for (int i = 0; i != rank; ++i) { 
        ofs2<< "\
      d[" << i << "] = u2 * g[" << i << "] - d[" << i << "] + x[of[" << i << "] + " << k-1 << "];" << endl;
      } 
      for (int i = rank; i != rank * 2; ++i) {
        ofs2<< "\
      d[" << i << "] = u2 * g[" << i << "] - d[" << i << "] + w[of[" << i - rank << "] + " << k-1 << "];" << endl;
      }
    }
  }
  for (int i = 0; i != rank; ++i) { 
    ofs2<< "\
      im[j][" << i << "] = uu * g[" << i <<  "] - d[" << i << "] + x[of[" << i << "]] * 0.5;" << endl;
  } 
  for (int i = rank; i != rank * 2; ++i) {
    ofs2<< "\
      im[j][" << i << "] = uu * g[" << i <<  "] - d[" << i << "] + w[of[" << (i - rank) << "]] * 0.5;" << endl;
  }
  ofs2<< "\
    }" << endl;
  for (int j = 1; j <= npoints; ++j) {
    for (int i = 1; i <= rank * 2; ++i) 
      ofs2<< "\
    d[" << i-1 << "] = 0.0;" << endl;

      mpreal pij = PI;
      pij *= (j - 1);
      const mpreal mp28 = npoints * 2;
      for (int k = 1; k <= npoints; ++k) {
        for (int i = 1; i <= rank * 2; ++i) {
          const double mpcos = (mpfr::cos(pij * (k + k - 1) / mp28)).toDouble();
          if (fabs(mpcos) < tiny) {
          ofs2<< "\
//  d[" << i-1 << "] = d[" << i-1 << "];" << endl;
          } else if (fabs(mpcos-1.0) < tiny) {
          ofs2<< "\
    d[" << i-1 << "] = d[" << i-1 << "] + im[" << k-1 << "][" << i-1 << "];" << endl;
          } else if (fabs(mpcos+1.0) < tiny) {
          ofs2<< "\
    d[" << i-1 << "] = d[" << i-1 << "] - im[" << k-1 << "][" << i-1 << "];" << endl;
          } else {
          ofs2<< "\
    d[" << i-1 << "] = d[" << i-1 << "] + im[" << k-1 << "][" << i-1 << "] * ("
                << setw(25) << setprecision(17) << scientific <<  mpcos << ");" << endl;
          }
       }
    }
    for (int i = 1; i <= rank * 2; ++i) {
      ofs2<< "\
    imc[" << j-1 << "][" << i-1 << "] = o7 * d[" << i-1 << "];" << endl;
    }
  }
  ofs2<< endl << "\
    const double t2 = tt * 2.0;" << endl;
  for (int i = 1; i <= rank * 2; ++i)
    ofs2<< "\
    d[" << i-1 << "] = 0.0;" << endl;
  for (int i = 1; i <= rank * 2; ++i) 
    ofs2<< "\
    g[" << i-1 << "] = 0.0;" << endl;

  for (int j = npoints; j >= 2; --j) {
    if (j % 2 == 0) {
      for (int i = 1; i <= rank * 2; ++i)
        ofs2<< "\
    g[" << i-1 << "] = t2 * d[" << i-1 << "] - g[" << i-1 << "] + imc[" << j-1 << "][" << i-1 << "];" << endl;
    } else {
      for (int i = 1; i <= rank * 2; ++i)
        ofs2<< "\
    d[" << i-1 << "] = t2 * g[" << i-1 << "] - d[" << i-1 << "] + imc[" << j-1 << "][" << i-1  << "];" << endl;
    }
  }
  ofs2<< endl;
  ofs2<< "\
    offset = (i - 1) * " << rank << ";" << endl;
  ofs2<< "\
    uu = 1.0 / sqrt(u);" << endl;
  for (int i = 1; i <= rank; ++i)
    ofs2<< "\
    rr[offset + " << i-1 << "] = tt * g[" << i-1 << "] - d[" << i-1 << "] + imc[0][" << i-1 << "] * 0.5;" << endl;
  for (int i = 1; i <= rank; ++i)
    ofs2<< "\
    ww[offset + " << i-1 << "] = (tt * g[" << (i + rank-1) << "] - d[" << (i + rank-1) << "] + imc[0][" << (i + rank-1) << "] * 0.5) * uu;" << endl;

  ofs2<< "\
  }" << endl << endl;
  ofs2<< "\
}" << endl;

  ofs << ofs2.str() << endl;
  ofs.close();
}

}
