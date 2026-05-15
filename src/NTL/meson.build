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

if use_gmp
    # gen_gmp_aux includes <NTL/mach_desc.h>, so its sources list must
    # carry mach_desc_h as a build-order dependency, and its
    # include_directories must reach the build-tree's NTL/ dir.
    # `meson.current_build_dir()/..` is `<builddir>/src`, which contains
    # the NTL/ subdir holding the generated headers.
    gen_gmp_aux_exe = executable('gen_gmp_aux',
        sources: ['../gen_gmp_aux.cpp', mach_desc_h],
        include_directories: ['..', '../../include'],
        cpp_args: ['-DNTL_GMP_LIP'],
        native: true,
        dependencies: [gmp_dep],
        install: false,
    )

    gmp_aux_writer = custom_target('gmp_aux.h',
        output: 'gmp_aux.h',
        input: mach_desc_h,
        command: [gen_gmp_aux_exe],
        capture: true,
        install: true,
        install_dir: get_option('includedir') / 'NTL',
    )
else
    # Fallback: emit a minimal stub so consumers' #include <NTL/gmp_aux.h>
    # resolves but defines no GMP-specific macros.
    gmp_aux_writer = custom_target('gmp_aux.h',
        output: 'gmp_aux.h',
        command: [python, files('../meson/gen-gmp-aux.py'), limb_size_bits.to_string()],
        capture: true,
        install: true,
        install_dir: get_option('includedir') / 'NTL',
    )
endif

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
