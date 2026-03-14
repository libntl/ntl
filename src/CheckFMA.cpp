
#include <NTL/ctools.h>

#include <cstdlib>
#include <iostream>


#include <NTL/simde_fma.h>


#if (NTL_BITS_PER_LONG != 64 || NTL_BITS_PER_INT != 32 || NTL_DOUBLE_PRECISION != 53)
#error "FMA not supported"
// sanity check -- code that uses this feature also relies on this
#endif

#ifndef NTL_HAVE_ALIGNED_ARRAY
#error "FMA not supported"
#endif

using namespace std;


void fun(double * x, const double *a, const double *b)
{
   __m256d xvec, avec, bvec, cvec;

   avec = _mm256_load_pd(a);
   bvec = _mm256_load_pd(b);
   xvec = _mm256_load_pd(x);

   xvec = _mm256_fmadd_pd(avec, bvec, xvec);

   _mm256_store_pd(x, xvec);
}

double power2(long k)
{
   long i;
   double res;

   res = 1;

   for (i = 1; i <= k; i++)
      res = res * 2;

   return res;
}


int main()
{
   NTL_AVX_LOCAL_ARRAY(vp, double, 12);

   double *a = vp + 0*4;
   double *b = vp + 1*4;
   double *x = vp + 2*4;

   a[0] = _ntl_nofold(1 + power2(NTL_DOUBLE_PRECISION-1));
   a[1] = _ntl_nofold(2);
   a[2] = _ntl_nofold(3);
   a[3] = _ntl_nofold(4);

   b[0] = _ntl_nofold(1 + power2(NTL_DOUBLE_PRECISION-1));
   b[1] = _ntl_nofold(3);
   b[2] = _ntl_nofold(4);
   b[3] = _ntl_nofold(5);

   x[0] = _ntl_nofold(-(1 + power2(NTL_DOUBLE_PRECISION-2))*power2(NTL_DOUBLE_PRECISION));
   x[1] = _ntl_nofold(4);
   x[2] = _ntl_nofold(5);
   x[3] = _ntl_nofold(6);

   // a[0] == 1 + 2^{52} 
   // b[0] == 1 + 2^{52}
   // x[0] == -(2^{53} + 2^{104})
   // a[0]*b[0] + x[0] == 1 if FMA
   //                  == 0 if not FMA

   fun(x, a, b);

   if (x[0] == 1 && x[1] == 10 && x[2] == 17 && x[3] == 26)
      return 0;
   else
      return -1;
}



