#include <NTL/ctools.h>

#include <cstdlib>
#include <emmintrin.h>
#include <tmmintrin.h>
#include <iostream>


#if (!defined(__GNUC__) || !defined(__x86_64__) || !defined(__SSSE3__))
#error "SSSE3 not supported"
#endif

#if (NTL_BITS_PER_LONG != 64 || NTL_BITS_PER_INT != 32 || NTL_DOUBLE_PRECISION != 53)
#error "SSSE3 not supported"
// sanity check -- code that uses this feature also relies on this
#endif


#ifndef NTL_HAVE_ALIGNED_ARRAY
#error "SSSE3 not supported"
#endif

#define ROL_VEC_8(x)	_mm_shuffle_epi8(x,_mm_set_epi8(14,13,12,15,10,9,8,11,6,5,4,7,2,1,0,3))
using namespace std;

void fun(unsigned int* x, const unsigned int *a)
{
   __m128i xvec, avec;

   avec = _mm_loadu_si128((const __m128i *) a);
   xvec = ROL_VEC_8(avec);

   _mm_storeu_si128((__m128i *)x, xvec);
}

int main()
{
   unsigned int a[4];
   unsigned int x[4];
 
   for (long i = 0; i < 4; i++) {
      a[i] = atoi("0") + i;
   }


   fun(x, a);

   for (long i = 0; i < 4; i++) {
      if (x[i] != 256*i) return -1;
   }

   return 0;
}



