#include <NTL/ctools.h>

#include <cstdlib>
#include <wmmintrin.h>
#include <iostream>

#if (!defined(__GNUC__) || !defined(__x86_64__) || !defined(__AVX__))
#error "PCLMUL not supported"
#endif

#if (NTL_BITS_PER_LONG != 64)
#error "PCLMUL not supported"
#endif

// NOTE: gcc and clang define __PCLMUL__, but icc does not

using namespace std;



static inline void
pclmul_mul1 (unsigned long *c, unsigned long a, unsigned long b)
{
   __m128i aa = _mm_setr_epi64( _mm_cvtsi64_m64(a), _mm_cvtsi64_m64(0));
   __m128i bb = _mm_setr_epi64( _mm_cvtsi64_m64(b), _mm_cvtsi64_m64(0));
   _mm_storeu_si128((__m128i*)c, _mm_clmulepi64_si128(aa, bb, 0));
}

int main()
{
   unsigned long a = ((unsigned long) atoi("15")) << (NTL_BITS_PER_LONG-4);
   unsigned long b = atoi("4");
   unsigned long c[2];

   pclmul_mul1(c, a, b);

   unsigned long c0 = ((unsigned long) atoi("3")) << (NTL_BITS_PER_LONG-2);
   unsigned long c1 = atoi("3");

   if (c[0] == c0 && c[1] == c1) 
      return 0;
   else
      return -1;
}
