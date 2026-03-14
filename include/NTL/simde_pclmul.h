#ifndef NTL_SIMDE_PCLMUL__H
#define NTL_SIMDE_PCLMUL__H

#if (defined(__GNUC__) && defined(__x86_64__) && defined(__AVX__))

#include <immintrin.h>

namespace {

inline void
pclmul_mul1 (unsigned long *c, unsigned long a, unsigned long b)
{
   __m128i aa = _mm_setr_epi64( _mm_cvtsi64_m64(a), _mm_cvtsi64_m64(0));
   __m128i bb = _mm_setr_epi64( _mm_cvtsi64_m64(b), _mm_cvtsi64_m64(0));
   _mm_storeu_si128((__m128i*)c, _mm_clmulepi64_si128(aa, bb, 0));
}

}



#elif (defined(__GNUC__) && defined(__aarch64__)) 

#include <arm_neon.h>
#include <cstring>

namespace {

inline void
pclmul_mul1 (unsigned long *c, unsigned long a, unsigned long b)
{
    poly64_t poly_a, poly_b;
    std::memcpy(&poly_a, &a, sizeof(poly64_t));
    std::memcpy(&poly_b, &b, sizeof(poly64_t));
    poly128_t product = vmull_p64(poly_a, poly_b);
    std::memcpy(c, &product, sizeof(product));
}

}


#else
#error "PCLMUL not available on this platform"
#endif

#endif
