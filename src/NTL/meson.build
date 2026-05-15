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

# ALIGNED_ARRAY: C++11 alignas/alignof. Required by C++11; assume
# present when cpp_std=c++11+.
have_extra_args += ['--present', 'ALIGNED_ARRAY']

# POSIX_TIME: clock_gettime, monotonic clock. Present on Linux,
# Darwin, FreeBSD; absent on bare Windows (MinGW has it via winpthreads).
if cpp.has_header_symbol('time.h', 'CLOCK_MONOTONIC')
    have_extra_args += ['--present', 'POSIX_TIME']
endif

# MACOS_TIME: mach_absolute_time. Darwin only.
if cpp.has_header('mach/mach_time.h')
    have_extra_args += ['--present', 'MACOS_TIME']
endif

# COPY_TRAITS2: __has_trivial_copy SFINAE form, older GCC. Modern
# compilers also support it. cc.compiles() probe matches NTL's
# CheckCOPY_TRAITS2.cpp's bare-essential check.
if cpp.compiles(
    '#include <type_traits>\nint main(void) { return __has_trivial_copy(int) ? 0 : 0; }',
    name: 'compiler has __has_trivial_copy',
)
    have_extra_args += ['--present', 'COPY_TRAITS2']
endif

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
