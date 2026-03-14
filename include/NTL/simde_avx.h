#ifndef NTL_SIMDE_AVX__H
#define NTL_SIMDE_AVX__H

#if (defined(__GNUC__) && defined(__x86_64__) && defined(__AVX__))

// On x86, just include the real AVX headers
#include <immintrin.h>


#elif (defined(__GNUC__) && defined(__aarch64__)) 

#include <arm_neon.h>

namespace {

// Mimic __m256d with a struct of two NEON registers
typedef struct {
    float64x2_t lo;
    float64x2_t hi;
} __m256d;

// Load 4 aligned doubles
inline __m256d _mm256_load_pd(const double* ptr) {
    __m256d result;
    result.lo = vld1q_f64(ptr);
    result.hi = vld1q_f64(ptr + 2);
    return result;
}

// Add 4 doubles
inline __m256d _mm256_add_pd(__m256d a, __m256d b) {
    __m256d result;
    result.lo = vaddq_f64(a.lo, b.lo);
    result.hi = vaddq_f64(a.hi, b.hi);
    return result;
}

// Multiply 4 doubles
inline __m256d _mm256_mul_pd(__m256d a, __m256d b) {
    __m256d result;
    result.lo = vmulq_f64(a.lo, b.lo);
    result.hi = vmulq_f64(a.hi, b.hi);
    return result;
}

// Store 4 aligned doubles
inline void _mm256_store_pd(double* ptr, __m256d v) {
    vst1q_f64(ptr, v.lo);
    vst1q_f64(ptr + 2, v.hi);
}

// Broadcast a single double to all 4 elements
inline __m256d _mm256_broadcast_sd(const double* ptr) {
    __m256d result;
    result.lo = vdupq_n_f64(*ptr);
    result.hi = vdupq_n_f64(*ptr);
    return result;
}

}


#else
#error "AVX not available on this platform"
#endif

#endif 
