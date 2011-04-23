//
// Author : Toru Shiozaki
// Date   : June 2009
//

#include "f77.h"
#include "mpreal.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace mpfr;

extern void find_root_weight(const mpreal T, const mpreal U, vector<mpreal>& dx, vector<mpreal>& dw, const int rank, const bool);

// testing...
void check(const mpreal Tvalue, const mpreal Uvalue) {

  const int rank = 2;
  fill2();

  vector<mpreal> resulta;
  vector<mpreal> resultb;
  {
    const double tarray = static_cast<double>(Tvalue);
    const double uarray = static_cast<double>(Uvalue);

    vector<double> rarray(rank);
    vector<double> warray(rank);
    const int unit = 1;

    root2(&tarray, &uarray, &rarray[0], &warray[0], unit);
    cout << setprecision(16) << scientific << setw(20) << rarray[0] << "  " << setw(20) << warray[0] << "  "
                                           << setw(20) << rarray[1] << "  " << setw(20) << warray[1] << endl;

    cout << "moments" << endl;
    for (int i = 0; i != rank * 2; ++i) {
      mpreal result;
      for (int j = 0; j != rank; ++j) result += mpfr::pow(rarray[j], i) * warray[j];
      resulta.push_back(result);
    } 
    for (int i = 0; i != rank; ++i) resultb.push_back(rarray[i]);
  }

  // verification code
  cout << "T" << " " << Tvalue << endl;
  cout << "U" << " " << Uvalue << endl;
  cout << "moments" << endl;

  {
    vector<mpreal> dx(rank);
    vector<mpreal> dw(rank);
    find_root_weight(Tvalue, Uvalue, dx, dw, rank, true);

    vector<mpreal> diff_moments(rank * 2);
    vector<mpreal> diff_position(rank);

    for (int i = 0; i != rank * 2; ++i) {
      mpreal result;
      for (int j = 0; j != rank; ++j)
        result += mpfr::pow(dx[j], i) * dw[j] / mpfr::sqrt(Uvalue);

      diff_moments[i] = abs(result - resulta[i]);
      cout << i << ":  " << setprecision(15) << fixed
           << setw(20) << (double)(result) << " "
           << setw(20) << (double)(result - resulta[i]) << " "
           << setw(20) << (double)((result - resulta[i]) / result)
           << endl; 
    } 
    for (int i = 0; i != rank; ++i) diff_position[i] = abs(dx[i] - resultb[i]); 

    cout << "max diff moment  : " << fixed << setprecision(15)
         << setw(20) << (double)(*max_element(diff_moments.begin(), diff_moments.end()))
         << endl 
         << endl;
  }

}

