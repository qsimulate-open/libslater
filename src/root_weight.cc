//
// Author: Toru Shiozaki
// Date  : June 2009
//

// absolutely for reference purpose
// therefore I will not pursue the performance.

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

#define SQRTTYPE mpfr::sqrt

using namespace std;
using namespace mpfr;


void find_root_weight(const mpreal T, const mpreal U, vector<mpreal>& dx, vector<mpreal>& dw, const int quad_rank, const bool mul) { 
  const mpreal zero = "0.0";
  const mpreal one  = "1.0";
  const mpreal two  = "2.0";
  const mpreal half = "0.5";

  mpfr::mpreal::set_default_prec(PREC);
  vector<mpreal> w(quad_rank); 
  vector<mpreal> x(quad_rank); 

  const int gtuend = 2 * quad_rank + 1;
  vector<mpreal> g_tu(gtuend + 1);
  const mpreal pi = PI;

  assert(T > zero);
  assert(U != zero);

  // some parameter
  const mpreal sqrtt = sqrt(T);
  const mpreal sqrtu = sqrt(U);
  const mpreal kappa  = - sqrtt + sqrtu;
  const mpreal lambda = sqrtt + sqrtu;

  // target Gm(T, U)
  const mpreal prefactor = exp(-T) * half * half * sqrt(pi); 
  g_tu[0] = prefactor / sqrt(U) * (exp(kappa * kappa) * erfc(kappa) + exp(lambda * lambda) * erfc(lambda)); 
  g_tu[1] = prefactor / sqrt(T) * (exp(kappa * kappa) * erfc(kappa) - exp(lambda * lambda) * erfc(lambda)); 

  const mpreal mone = g_tu[0] * sqrt(U);
  const mpreal expmt = exp(-T);
  const mpreal twoU = U + U;
  const mpreal halfpT = half / T;
  for (int i = 2; i != gtuend + 1; ++i)
    g_tu[i] = halfpT * (static_cast<mpreal>(2 * i - 3) * g_tu[i - 1] + twoU * g_tu[i - 2] - expmt);


  // Chebyshev algorithm   
  const int n = quad_rank; 
  vector<mpreal> mlp(2 * n + 2);
  vector<mpreal> sigma1 = g_tu;
  vector<mpreal> sigma0(2 * n + 2, zero);

  mlp[0] = one / sigma1[0];
  x[0] = sigma1[1] * mlp[0];
  w[0] = zero;
 
  for (int k = 0; k <= n - 2; k += 2) {
    for (int l = k; l <= 2 * n - k - 3; ++l) {
      sigma0[l + 1] = sigma1[l + 2] - x[k] * sigma1[l + 1] - w[k] * sigma0[l + 1] ;
    }
    mlp[k + 1] = one / sigma0[k + 1];
    x[k + 1] = - sigma1[k + 1] * mlp[k] + sigma0[k + 2] * mlp[k + 1];
    w[k + 1] = sigma0[k + 1] * mlp[k];

    if(k != n - 2) { 
      for (int l = k + 1; l <= 2 * n - k - 2; ++l)  
        sigma1[l + 1] = sigma0[l + 2] - x[k + 1] * sigma0[l + 1] - w[k + 1] * sigma1[l + 1];
      mlp[k + 2] = one / sigma1[k + 2];
      x[k + 2] = - sigma0[k + 2] * mlp[k + 1] + sigma1[k + 3] * mlp[k + 2];
      w[k + 2] = sigma1[k + 2] * mlp[k + 1];
    }
  }
  // projecting to mpreal if needed;
  dx[0] = x[0];
  for (int i = 1; i != n; ++i) {
    const mpreal tmp = w[i];
    dw[i - 1] = SQRTTYPE(tmp);
    dx[i] = x[i]; 
  }

//////////
  // solve tri-diagonal linear equation 
    
  vector<mpreal> lp(2 * n);
  lp[0] = one;
  for (int i = 1; i <= n + n - 2; ++i) lp[i] = zero;

  for (int l = 0; l <= n - 1; ++l) {
    int iter = 0;
line1:
    int mm;
    for (mm = l; mm <= n - 2; ++mm) { // ?????
      mpreal dd = fabs(dx[mm]) + fabs(dx[mm + 1]);
      if (fabs(dw[mm]) + dd == dd) goto line2;
    }
    mm = n - 1; // ??
line2:
    if(mm != l) {
      ++iter;
      mpreal g = (dx[l + 1] - dx[l]) / (dw[l] * 2);
      const mpreal r = SQRTTYPE(g * g + one);
      g = dx[mm] - dx[l] + dw[l] / (g + g >= 0 ? fabs(r) : -fabs(r));
      mpreal s = one;
      mpreal c = one;
      mpreal p = zero;
      for (int i = mm - 1; i >= l; --i) {
        mpreal f = s * dw[i];
        mpreal bb = c * dw[i];
        mpreal r = SQRTTYPE(f * f + g * g);
        dw[i + 1] = r;
        if(r == zero) {
          dx[i + 1] = dx[i + 1] - p;
          dw[mm] = zero;
          goto line1;
        }
        s = f / r;
        c = g / r;
        g = dx[i + 1] - p;
        mpreal td = c * bb;
        r = (dx[i] - g) * s + td * 2;
        p = s * r;
        dx[i + 1] = g + p;     
        g = c * r - bb;
        f= lp[i + 1];
        lp[i+1] = s * lp[i] + c * f;
        lp[i] = c * lp[i] - s * f;
      }
      dx[l] = dx[l] - p;
      dw[l] = g;
      dw[mm] = zero;
      goto line1;
    }
  }
  for (int i = 0; i <= n - 1; ++i) dw[i] = lp[i] * lp[i] * (mul ? mone : one);

  map<mpreal, mpreal> rw;
  for (int i = 0; i != n; ++i) rw.insert(make_pair(dx[i], dw[i]));
  int id = 0;
  for (map<mpreal, mpreal>::const_iterator miter = rw.begin(); miter != rw.end(); ++miter, ++id) {
    dx[id] = miter->first;
    dw[id] = miter->second;
  }
}
