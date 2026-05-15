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
# and runs Check<feature>.cpp probes — not cross-compile-safe). For MVP
# we emit empty stubs (every feature absent). A compile-time probing
# pass via cc.compiles() is a polish-phase follow-up.
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
    command: [python, files('../meson/gen-have-headers.py'),
              meson.current_build_dir()],
    install: true,
    install_dir: get_option('includedir') / 'NTL',
)
