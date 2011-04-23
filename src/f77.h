#ifndef __src_rysint_f77_h
#define __src_rysint_f77_h

extern "C" {
 void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
 void dcopy_(const int*, const double*, const int*, double*, const int*);

 void root7_(const double*, const double*, double*, double*, const int*); 
};
void fill2();
void root2(const double*, const double*, double*, double*, const int); 

#endif
