#include <NTL/ctools.h>

#include <cstdlib>
#include <iostream>

#include <wmmintrin.h>
#include <emmintrin.h>
#include <tmmintrin.h>

using namespace std;

#if (NTL_BITS_PER_LONG != 64)
#error "NTL_BITS_PER_LONG != 64"
#endif

int main()
{
  __m128i a=_mm_cvtsi64x_si128(atol("17"));
  __m128i key=_mm_cvtsi64x_si128(atol("42"));
  a = _mm_aesenclast_si128(a, key);
  long x = _mm_cvtsi128_si64x(a);
  if (x != atol("7161677110969590696")) return -1;

  return 0;
}
