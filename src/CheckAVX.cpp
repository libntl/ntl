#include <NTL/ctools.h>

#include <cstdlib>
#include <immintrin.h>
#include <iostream>


#if (!defined(__GNUC__) || !defined(__x86_64__) || !defined(__AVX__))
#error "AVX not supported"
#endif

#if (NTL_BITS_PER_LONG != 64 || NTL_BITS_PER_INT != 32 || NTL_DOUBLE_PRECISION != 53)
#error "AVX not supported"
// sanity check -- code that uses this feature also relies on this
#endif

#ifndef NTL_HAVE_ALIGNED_ARRAY
#error "AVX not supported"
#endif

using namespace std;

void fun(double * x, const double *a, const double *b)
{
   __m256d xvec, avec, bvec, cvec;

   avec = _mm256_load_pd(a);
   bvec = _mm256_load_pd(b);
   xvec = _mm256_load_pd(x);

   xvec = _mm256_add_pd(_mm256_mul_pd(avec, bvec), xvec);

   _mm256_store_pd(x, xvec);
}

int main()
{
   NTL_AVX_LOCAL_ARRAY(vp, double, 12);

   double *a = vp + 0*4;
   double *b = vp + 1*4;
   double *x = vp + 2*4;

   a[0] = atoi("1");
   a[1] = atoi("2");
   a[2] = atoi("3");
   a[3] = atoi("4");

   b[0] = atoi("2");
   b[1] = atoi("3");
   b[2] = atoi("4");
   b[3] = atoi("5");

   x[0] = atoi("3");
   x[1] = atoi("4");
   x[2] = atoi("5");
   x[3] = atoi("6");

   fun(x, a, b);

   if (x[0] == 5 && x[1] == 10 && x[2] == 17 && x[3] == 26)
      return 0;
   else
      return -1;
}



