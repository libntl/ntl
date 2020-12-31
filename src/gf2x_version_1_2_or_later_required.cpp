#include <NTL/config.h>

#if (defined(NTL_GF2X_LIB) && defined(NTL_THREADS))
// we require v1.2 or later

#include <gf2x.h>

#ifndef GF2X_VERSION_MAJOR
// versions after v1.2 should define GF2X_VERSION_MAJOR

extern "C" {

struct gf2x_ternary_fft_info_s;

typedef const struct gf2x_ternary_fft_info_s * gf2x_ternary_fft_info_srcptr;

int gf2x_ternary_fft_compatible(gf2x_ternary_fft_info_srcptr o1, gf2x_ternary_fft_info_srcptr o2);

}

// This will fail to link for versions prior to v1.2.
void fun()
{
   gf2x_ternary_fft_compatible(0, 0);
}

#endif

#endif

int main()
{
   return 0;
}
