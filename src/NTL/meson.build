# Generated NTL/ headers live here so that NTL's source code's
# `#include <NTL/config.h>` (and friends) resolves to the build-tree file.
# include_directories('.') from src/meson.build makes this path visible.

mach_desc_h = custom_target('mach_desc.h',
    output: 'mach_desc.h',
    command: [python, files('../meson/run-makedesc.py'), makedesc],
    capture: true,
    install: true,
    install_dir: get_option('includedir') / 'NTL',
)

# gmp_aux.h is produced by a small Python generator rather than NTL's
# native gen_gmp_aux.cpp executable. The C++ generator runs on the build
# host, links against the build host's GMP, and aborts under cross-compile
# whenever the target's bits_per_long differs from the host's. Computing
# the same macros from the target-side `cc.sizeof('mp_limb_t')` value works
# correctly for both native and cross builds.
gmp_aux_writer = custom_target('gmp_aux.h',
    output: 'gmp_aux.h',
    command: [
        python, files('../meson/gen-gmp-aux.py'),
        limb_size_bits.to_string(), abi['bits_per_long'],
    ],
    capture: true,
    install: true,
    install_dir: get_option('includedir') / 'NTL',
)

config_h = configure_file(
    input: '../config.h.in',
    output: 'config.h',
    configuration: cfg,
    install: true,
    install_dir: get_option('includedir') / 'NTL',
)

# ALL_FEATURES.h #includes NTL/HAVE_<feature>.h for 16 features. The
# Makefile build generates these via MakeCheckFeatures (which compiles
# and runs Check<feature>.cpp probes — not cross-compile-safe). We do
# compile-time probes here via cc.compiles() and hand the discovered
# features to gen-have-headers.py. COPY_TRAITS1 and CHRONO_TIME are
# unconditionally marked present inside the script (required by C++11
# semantics that NTL_SAFE_VECTORS depends on).
have_extra_args = []

# LL_TYPE: NTL uses __int128 as NTL_ULL_TYPE on GCC 4+. The probe
# mirrors NTL's own ctools.h gate.
if cpp.compiles(
    '__int128 x = 1; int main(void) { (void)x; return 0; }',
    name: 'compiler has __int128 (LL_TYPE)',
)
    have_extra_args += ['--present', 'LL_TYPE']
endif

# BUILTIN_CLZL: GCC / clang builtin.
if cpp.compiles(
    'int main(void) { return __builtin_clzl(1UL); }',
    name: 'compiler has __builtin_clzl',
)
    have_extra_args += ['--present', 'BUILTIN_CLZL']
endif

# ALIGNED_ARRAY is intentionally NOT enabled. NTL's _ntl_make_aligned
# implementation casts char* through NTL_UPTRINT_T, which is
# unsigned long unless NTL_BIG_POINTERS is set in mach_desc.h. On
# LLP64 targets like x86_64-w64-mingw32, long is 32-bit and pointers
# are 64-bit — the cast loses precision. NTL_BIG_POINTERS detection is
# driven by MakeDesc running on the build host, which sees the build
# host's pointer/long relationship, not the target's. Until we plumb
# target-specific NTL_BIG_POINTERS via the ABI table, ALIGNED_ARRAY
# stays off to keep the build correct on cross-LLP64 targets.

# POSIX_TIME, MACOS_TIME: build-time time-source probes. Intentionally
# not enabled — would need bootstrapping through ctools.h which itself
# depends on mach_desc.h. A polish-phase follow-up adds these via
# header-only probes that bypass ctools.h.

# COPY_TRAITS2 is intentionally NOT enabled. NTL's NTL_SAFE_VECTORS
# code path uses COPY_TRAITS1 when both are available; enabling
# COPY_TRAITS2 in addition would force the older SFINAE-based form
# which trips on certain modern compilers' strict template handling.

have_headers = [
    'HAVE_ALIGNED_ARRAY.h', 'HAVE_BUILTIN_CLZL.h', 'HAVE_LL_TYPE.h',
    'HAVE_SSSE3.h', 'HAVE_AVX.h', 'HAVE_PCLMUL.h',
    'HAVE_AVX2.h', 'HAVE_FMA.h', 'HAVE_AVX512F.h',
    'HAVE_COPY_TRAITS1.h', 'HAVE_COPY_TRAITS2.h',
    'HAVE_CHRONO_TIME.h', 'HAVE_MACOS_TIME.h', 'HAVE_POSIX_TIME.h',
    'HAVE_AES_NI.h', 'HAVE_KMA.h',
]
have_target = custom_target('have-headers',
    output: have_headers,
    command: [
        python, files('../meson/gen-have-headers.py'),
        meson.current_build_dir(),
        have_extra_args,
    ],
    install: true,
    install_dir: get_option('includedir') / 'NTL',
)
