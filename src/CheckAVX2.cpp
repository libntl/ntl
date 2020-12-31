#include <NTL/ctools.h>

#include <cstdlib>
#include <immintrin.h>
#include <iostream>


#if (!defined(__GNUC__) || !defined(__x86_64__) || !defined(__AVX2__))
#error "AVX2 not supported"
#endif

#if (NTL_BITS_PER_LONG != 64 || NTL_BITS_PER_INT != 32 || NTL_DOUBLE_PRECISION != 53)
#error "AVX2 not supported"
// sanity check -- code that uses this feature also relies on this
#endif


#ifndef NTL_HAVE_ALIGNED_ARRAY
#error "AVX2 not supported"
#endif

using namespace std;

void fun(unsigned int* x, const unsigned int *a)
{
   __m256i xvec, avec;

   avec = _mm256_loadu_si256((const __m256i *) a);
   xvec = _mm256_slli_epi32(avec, 5);

   _mm256_storeu_si256((__m256i *)x, xvec);
}

int main()
{
   unsigned int a[8];
   unsigned int x[8];
 
   for (long i = 0; i < 7; i++) {
      a[i] = atoi("0") + i;
   }


   fun(x, a);

   for (long i = 0; i < 7; i++) {
      if (x[i] != 32*i) return -1;
   }

   return 0;
}



