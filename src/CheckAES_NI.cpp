#include <NTL/ctools.h>

#include <cstdlib>
#include <iostream>

#include <wmmintrin.h>
#include <emmintrin.h>
#include <tmmintrin.h>

using namespace std;

int main()
{
   __m128i out, rkeys[16] = {0}, nv;
   __m128i temp = _mm_xor_si128(nv, rkeys[0]);
   int i;

   for (i = 1 ; i < 14 ; i++) {
      temp = _mm_aesenc_si128(temp, rkeys[i]);
   }
   temp = _mm_aesenclast_si128(temp, rkeys[14]);
   _mm_store_si128(&out, temp);
   return 0;
}
