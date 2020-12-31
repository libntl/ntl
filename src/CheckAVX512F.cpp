
#include <NTL/ctools.h>

#include <cstdlib>
#include <immintrin.h>
#include <iostream>

// Ths actually checks for AVX512F+DQ+VL


#if (!defined(__GNUC__) || !defined(__x86_64__) || !defined(__AVX512F__))
#error "AVX512F not supported"
#endif

#if (!defined(__AVX512VL__) || !defined(__AVX512DQ__))
#error "AVX512F not supported"
#endif

#if (NTL_BITS_PER_LONG != 64 || NTL_BITS_PER_INT != 32 || NTL_DOUBLE_PRECISION != 53)
#error "AVX512F not supported"
// sanity check -- code that uses this feature also relies on this
#endif

#ifndef NTL_HAVE_ALIGNED_ARRAY
#error "AVX512F not supported"
#endif

using namespace std;


void fun(double * x, const double *a, const double *b)
{
   __m512d xvec, avec, bvec, cvec;

   avec = _mm512_load_pd(a);
   bvec = _mm512_load_pd(b);
   xvec = _mm512_load_pd(x);

   xvec = _mm512_fmadd_pd(avec, bvec, xvec);

   _mm512_store_pd(x, xvec);
}

void fun1(double *x, const long *p)
{
   __m256i a = _mm256_load_si256((const __m256i*)p);
   _mm256_store_pd(x,  _mm256_cvtepi64_pd(a));
}


int main()
{
   NTL_AVX512_LOCAL_ARRAY(vp, double, 24);

   double *a = vp + 0*8;
   double *b = vp + 1*8;
   double *x = vp + 2*8;

   a[0] = atoi("1");
   a[1] = atoi("2");
   a[2] = atoi("3");
   a[3] = atoi("4");
   a[4] = atoi("5");
   a[5] = atoi("6");
   a[6] = atoi("7");
   a[7] = atoi("8");

   b[0] = atoi("2");
   b[1] = atoi("3");
   b[2] = atoi("4");
   b[3] = atoi("5");
   b[4] = atoi("6");
   b[5] = atoi("7");
   b[6] = atoi("8");
   b[7] = atoi("9");

   x[0] = atoi("3");
   x[1] = atoi("4");
   x[2] = atoi("5");
   x[3] = atoi("6");
   x[4] = atoi("7");
   x[5] = atoi("8");
   x[6] = atoi("9");
   x[7] = atoi("10");

   fun(x, a, b);

   NTL_AVX_LOCAL_ARRAY(lp, long, 4);
   NTL_AVX_LOCAL_ARRAY(dp, double, 4);

   lp[0] = atoi("1");
   lp[1] = atoi("2");
   lp[2] = atoi("3");
   lp[3] = atoi("4");
  
   fun1(dp, lp);

   if (x[0] ==  5 && x[1] == 10 && x[2] == 17 && x[3] == 26 &&
       x[4] == 37 && x[5] == 50 && x[6] == 65 && x[7] == 82 &&
       dp[0] == 1 && dp[1] == 2 && dp[2] == 3 && dp[3] == 4)
      return 0;
   else
      return -1;
}



