#ifndef NTL_ZZ_limbs__H
#define NTL_ZZ_limbs__H

#include <NTL/ZZ.h>


// NOTE: unlike other NTL header files, this one needs access
// to GMP's header file, which means that C++ files that include
// this file will need to ensure that the compiler has the 
// right "include path" to get at GMP's header file.




#ifdef NTL_GMP_LIP
#include <gmp.h>

typedef mp_limb_t _ntl_limb_t;

#else

typedef unsigned long _ntl_limb_t;

// #define NTL_BITS_PER_LIMB_T NTL_BITS_PER_LONG
// This is already defined in ZZ.h 


#endif

void _ntl_glimbs_set(const _ntl_limb_t *p, long n, _ntl_gbigint *x);

// DIRT: This exposes some internals that should be in lip.cpp,
// but are here to make it inline.
inline 
const _ntl_limb_t * _ntl_glimbs_get(_ntl_gbigint p)
   { return p ? ((_ntl_limb_t *) (((long *) (p)) + 2)) : 0; }


NTL_OPEN_NNS

typedef _ntl_limb_t ZZ_limb_t;


inline 
void ZZ_limbs_set(ZZ& x, const ZZ_limb_t *p, long n)
{
   _ntl_glimbs_set(p, n, &x.rep);
}

inline
const ZZ_limb_t * ZZ_limbs_get(const ZZ& a)
{
   return _ntl_glimbs_get(a.rep);
}


NTL_CLOSE_NNS


#endif
