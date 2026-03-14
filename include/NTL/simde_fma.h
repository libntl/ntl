#ifndef NTL_SIMDE_FMA__H
#define NTL_SIMDE_FMA__H

#if (defined(__GNUC__) && defined(__x86_64__) && defined(__AVX2__))

// On x86, just include the real AVX headers
#include <immintrin.h>


#elif (defined(__GNUC__) && defined(__aarch64__)) 

#include <NTL/simde_avx.h>

namespace {

// Fused multiply-add: a*b + c
inline __m256d _mm256_fmadd_pd(__m256d a, __m256d b, __m256d c) {
    __m256d result;
    // vfmaq_f64(c, a, b) computes c + a*b
    result.lo = vfmaq_f64(c.lo, a.lo, b.lo);
    result.hi = vfmaq_f64(c.hi, a.hi, b.hi);
    return result;
}

}

#else
#error "FMA not available on this platform"
#endif


#endif
