//
// Toru Shiozaki
// June 2009
//

#define USE_DOUBLE

#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "mpreal.h"
#include "local_macros.h"
#include <algorithm>
#include <vector>
#include <map>
#include "f77.h"

using namespace std;
using namespace mpfr;

extern void find_root_weight(const mpreal T, const mpreal U, vector<mpreal>& dx, vector<mpreal>& dw, const int rank, const bool);


vector<mpreal> chebft(int n) {
  mpfr::mpreal::set_default_prec(PREC);
  vector<mpreal> out(n);
  for (int k = 0; k != n; ++k) {
     const mpreal y = mpfr::cos(PI * static_cast<mpreal>(k + 0.5) / n);
     out[k] = y;
  }
  return out;
}


vector<vector<double> > get_C(const int tbase, const int ubase, int rank) {
  mpfr::mpreal::set_default_prec(PREC);
  const int n = NGRID;

  const mpreal zero = "0.0";
  const mpreal half = "0.5";
  vector<mpreal> cheb = chebft(n); 

  const mpreal Umin = ubase;
  const mpreal Umax = Umin + 1;
  const mpreal Up = half * (Umin + Umax);
  const mpreal Tmin = tbase;
  const mpreal Tmax = Tmin + 1;
  const mpreal Tp = half * (Tmin + Tmax);

  vector<mpreal> Upoints(n);
  vector<mpreal> Tpoints(n);
  for (int i = 0; i != n; ++i) { 
    Upoints[i] = half * cheb[i] + Up; 
    Tpoints[i] = half * cheb[i] + Tp;
  }

  vector<pair<vector<mpreal>, vector<mpreal> > > table_reserve;
  table_reserve.resize(n * n);
  int count = 0;
  for (vector<mpreal>::const_iterator titer = Tpoints.begin(); titer != Tpoints.end(); ++titer) { 
    const int iall = Upoints.size();
    #pragma omp parallel for
    for (int i = 0; i < iall; ++i ) { 
      vector<mpreal>::const_iterator uiter = Upoints.begin() + i;
      vector<mpreal> dx(rank);
      vector<mpreal> dw(rank);
      find_root_weight(tbase==0 ? *titer * *titer : mpfr::pow("3", *titer-1.0), mpfr::pow("10", *uiter), dx, dw, rank, true);
      table_reserve[count + i] = make_pair(dx, dw);
    }
    count += n;
  }

  vector<vector<double> > c;
  for (int ii = 0; ii != rank; ++ii) {
    vector<double> tc(n * n);
    vector<double> tc2(n * n);

    vector<mpreal> cdx, cdw;
    for (int j = 0; j != n * n; ++j) {
      cdx.push_back((table_reserve[j].first)[ii]);
      cdw.push_back((table_reserve[j].second)[ii]);
    }  

    const mpreal two = "2.0";
    for (int i = 0; i != n; ++i) {
      const mpreal fac = two / n; 
      const mpreal pi = PI;
      for (int j = 0; j != n; ++j) {
        mpreal sum = 0.0;
        mpreal sum2 = 0.0;
        for (int k = 0; k != n; ++k) {
          sum += cdx[i * n + k] * cos(pi * j * (k + 0.5) / n);
          sum2 += cdw[i * n + k] * cos(pi * j * (k + 0.5) / n);
        }
        tc[i * n + j] = sum * fac;
        tc2[i * n + j] = sum2 * fac;
      }
    }
    c.push_back(tc);
    c.push_back(tc2);
  }
  
  return c;
}

