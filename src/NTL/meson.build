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
    # gen_gmp_aux needs mach_desc.h available (it includes <NTL/mach_desc.h>).
    # Run it via a python wrapper that:
    #   1. cd's into a tempdir
    #   2. exposes the build-tree NTL/ dir on the include path via env or
    #      symlink — but gen_gmp_aux only needs include resolution at compile
    #      time, which is already provided to the executable target above.
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
