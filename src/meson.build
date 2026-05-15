# T022: Library, generated headers, and tests for NTL under Meson.
#
# Consumes:
#   - The `abi` dict from meson.build (resolved by src/meson/pick-abi.py).
#   - meson.options toggles.
# Produces:
#   - mach_desc.h via a custom_target that runs MakeDesc on the build host.
#   - gmp_aux.h and config.h via configure_file.
#   - libntl shared library (and static if requested) from the sources listed
#     in src/meson/sources.txt (mechanically synced from src/mfile).
#   - The QuickTest / BerlekampTest / ZZTest programs under meson test.
#   - An ntl.pc pkg-config file (FR-006).

# -------------------------------------------------------------------
# Native MakeDesc executable. Built with native: true so it always runs on the
# build host regardless of cross status. FORCE_BPL / FORCE_NO_FMA are passed
# according to the resolved ABI table (R-002).
# -------------------------------------------------------------------

makedesc_args = [
    '-DNTL_FORCE_BPL=' + abi['bits_per_long'],
]
if abi['fma_policy'] == 'off'
    makedesc_args += ['-DNTL_FORCE_NO_FMA']
endif

# Note: native: true tells Meson to build with the build-host compiler, not the
# host (target) compiler. This is the cross-compile-safe path documented in the
# roadmap (Phase 1 §c: build-host code generation).
makedesc = executable('MakeDesc',
    sources: ['MakeDesc.cpp', 'MakeDescAux.cpp'],
    include_directories: ['../include'],
    cpp_args: makedesc_args,
    native: true,
    install: false,
)

# Generated headers (mach_desc.h, gmp_aux.h, config.h) are declared in
# src/NTL/meson.build so the build path matches NTL's `#include <NTL/foo.h>`
# convention. They share state with src/meson.build via the variables set
# above (cfg, limb_size_bits, makedesc, python).

# -------------------------------------------------------------------
# gmp_aux.h — produced via cc.sizeof('mp_limb_t') (compile-time, cross-safe).
# -------------------------------------------------------------------

if use_gmp
    limb_size_bits = cpp.sizeof('mp_limb_t',
        prefix: '#include <gmp.h>',
        dependencies: gmp_dep,
    ) * 8
    if limb_size_bits < 1
        error('Could not determine sizeof(mp_limb_t). Is GMP installed for the target?')
    endif

    # gen_gmp_aux is declared inside src/NTL/meson.build, AFTER mach_desc.h
    # has been declared as a custom_target, so its sources list can include
    # mach_desc.h and Meson schedules the generation before the compile.
else
    limb_size_bits = abi['bits_per_long'].to_int()
endif

# (Generator targets moved to src/NTL/meson.build — see above.)

# -------------------------------------------------------------------
# config.h — produced from src/config.h.in via configure_file.
# -------------------------------------------------------------------

cfg = configuration_data()

# configure_file substitutes @KEY@ with the value; the `#if @KEY@` idiom
# in config.h.in evaluates 1 / 0 directly.

cfg.set('NTL_LEGACY_NO_NAMESPACE', get_option('legacy_no_namespace') ? 1 : 0)
cfg.set('NTL_LEGACY_INPUT_ERROR',  get_option('legacy_input_error')  ? 1 : 0)
cfg.set('NTL_TLS_HACK',            (get_option('tls_hack').enabled() or
                                    (get_option('tls_hack').auto() and abi['tls_hack'] == 'true')) ? 1 : 0)
cfg.set('NTL_THREADS',             use_threads ? 1 : 0)
cfg.set('NTL_EXCEPTIONS',          get_option('exceptions') ? 1 : 0)
cfg.set('NTL_THREAD_BOOST',        (get_option('thread_boost') and use_threads) ? 1 : 0)
cfg.set('NTL_GMP_LIP',             use_gmp ? 1 : 0)
cfg.set('NTL_GF2X_LIB',            use_gf2x ? 1 : 0)
cfg.set('NTL_STD_CXX11',           1)  # NTL is C++11+ era; sanity checks in cfile require it.
cfg.set('NTL_STD_CXX14',           0)
cfg.set('NTL_DISABLE_MOVE_ASSIGN', 0)
cfg.set('NTL_DISABLE_MOVE',        0)
cfg.set('FLAG_UNSIGNED_LONG_LONG_TYPE', 0)
cfg.set('NTL_UNSIGNED_LONG_LONG_TYPE',  '')
cfg.set('NTL_CLEAN_INT',           0)
cfg.set('NTL_CLEAN_PTR',           0)
cfg.set('NTL_SAFE_VECTORS',        get_option('safe_vectors') ? 1 : 0)
# NTL_ENABLE_AVX_FFT requires either AVX-512F or (AVX2 && FMA) at runtime.
# Probing those reliably is Phase 3+ scope; default off for MVP. ABI table's
# x86_specializations governs the AVX FFT source files included, not this
# define.
cfg.set('NTL_ENABLE_AVX_FFT',      0)
cfg.set('NTL_AVOID_AVX512',        0)
cfg.set('NTL_RANGE_CHECK',         get_option('range_check') ? 1 : 0)
cfg.set('NTL_NO_INIT_TRANS',       0)
cfg.set('NTL_X86_FIX',             0)
cfg.set('NTL_NO_X86_FIX',          0)
cfg.set('NTL_LEGACY_SP_MULMOD',    0)
disable_ld = (get_option('disable_longdouble').enabled() or
              (get_option('disable_longdouble').auto() and abi['long_double'] == 'disable'))
cfg.set('NTL_DISABLE_LONGDOUBLE',  disable_ld ? 1 : 0)
cfg.set('NTL_DISABLE_LONGLONG',    0)
cfg.set('NTL_DISABLE_LL_ASM',      0)
cfg.set('NTL_MAXIMIZE_SP_NBITS',   0)
cfg.set('NTL_SPMM_ULL',            0)
cfg.set('NTL_FFT_BIGTAB',          0)
cfg.set('NTL_FFT_LAZYMUL',         0)
cfg.set('NTL_AVOID_BRANCHING',     0)
cfg.set('NTL_TBL_REM',             0)
# NTL_CRT_ALTCODE: activates _ntl_crt_struct_tbl in lip.cpp — a table-
# driven CRT path that's faster on x86 and only viable when LL_TYPE is
# available with NTL_BITS_PER_LIMB_T == NTL_BITS_PER_LONG. The Makefile
# build's `./configure` enables this by default on x86; mirror that
# heuristic via the ABI table's x86_specializations field.
cfg.set('NTL_CRT_ALTCODE',         (abi['x86_specializations'] == 'true') ? 1 : 0)
cfg.set('NTL_CRT_ALTCODE_SMALL',   0)
cfg.set('NTL_GF2X_ALTCODE',        0)
cfg.set('NTL_GF2X_ALTCODE1',       0)
cfg.set('NTL_GF2X_NOINLINE',       0)
cfg.set('NTL_RANDOM_AES256CTR',    0)
# Wizard never runs under Meson — emit a comment in its place.
cfg.set('WIZARD_HACK',
        '/* Meson build: WIZARD_HACK intentionally empty (Wizard not run). */')

# cfg dict is now built; the actual configure_file lives in src/NTL/meson.build.
subdir('NTL')

# -------------------------------------------------------------------
# Library sources.
# -------------------------------------------------------------------

src_list_text = run_command(
    python, files('../tools/sync-sources.py'),
    check: true,
).stdout().strip().split('\n')

# Resolve each name to an absolute file under src/.
library_sources = []
foreach src_name : src_list_text
    library_sources += files(src_name)
endforeach

# GetTime and GetPID are selected at configure time by the Makefile build's
# MakeGetTime / MakeGetPID probe scripts. For Meson we hard-pick:
#   - GetTime5.cpp: std::chrono based, requires C++11 (which NTL guarantees).
#   - GetPID1.cpp:  POSIX getpid(); use GetPID2.cpp on Windows where no posix.
# Future targets (Windows in Phase 7) will adjust this selection via the ABI
# table. For Phase 3 MVP (Linux native) the chrono+POSIX combo works.
if os_name == 'windows'
    library_sources += files('GetTime5.cpp', 'GetPID2.cpp')
else
    library_sources += files('GetTime5.cpp', 'GetPID1.cpp')
endif

# Public include dir: ../include from this file's perspective.
ntl_inc = include_directories('../include')

# -------------------------------------------------------------------
# libntl target.
# -------------------------------------------------------------------

library_link_args = []
if abi['shlib_style'] == 'dll'
    # Coarse-grained export (R-004) — refine later via NTL_API macro work.
    library_link_args += ['-Wl,--export-all-symbols']
endif

library_deps = []
if use_gmp
    library_deps += [gmp_dep]
endif
if use_gf2x
    library_deps += [gf2x_dep]
endif
if use_threads
    library_deps += [threads_dep]
endif

# Build the library. The generated headers (mach_desc.h, gmp_aux.h, config.h,
# HAVE_*.h) are listed as sources so Meson schedules them before any .cpp
# compilation.
libntl = library('ntl',
    library_sources + [mach_desc_h, gmp_aux_writer, config_h, have_target],
    include_directories: [ntl_inc, include_directories('.')],
    dependencies: library_deps,
    cpp_args: ['-DNTL_BUILD'],
    link_args: library_link_args,
    install: true,
    soversion: '0',
)

# -------------------------------------------------------------------
# pkg-config (FR-006).
# -------------------------------------------------------------------

pkg = import('pkgconfig')
pkg.generate(libntl,
    name: 'ntl',
    description: 'NTL: a Library for doing Number Theory',
    version: meson.project_version(),
    # GMP often ships without a .pc file (e.g. Debian's libgmp-dev). Avoid
    # `requires_private: ['gmp']` and just emit `-lgmp` instead.
    libraries_private: use_gmp ? ['-lgmp'] : [],
)

# -------------------------------------------------------------------
# Tests (FR-007).
# -------------------------------------------------------------------

ntl_test_dep = declare_dependency(
    link_with: libntl,
    include_directories: [ntl_inc, include_directories('.')],
    dependencies: library_deps,
    # Generated headers must exist before any consumer .cpp compiles.
    # Listing them as `sources` of the dependency makes Meson schedule
    # the generator targets first.
    sources: [mach_desc_h, gmp_aux_writer, config_h, have_target],
)

# Translation of exec_mode → Meson `should_run` and exe_wrapper handling.
tests_should_run = abi['exec_mode'] != 'cross-only'

# Tests come in two shapes:
#   - simple: program is run with no stdin, just exit code matters.
#     (e.g. QuickTest, ZZTest)
#   - golden-diff: program reads <name>In from stdin, output is compared
#     with `diff -b` against <name>Out. (BerlekampTest, and many others
#     in src/TestScript.)
#
# We BUILD QuickTest and ZZTest binaries so users can run them locally
# (matching what NTL's `make check` does), but we DON'T register them
# under `meson test`. QuickTest in particular is a benchmark that loops
# at sizes up to 2^18 doubling iteration counts until each measurement
# runs >=0.5s — typical wall time is 30-60min natively and several hours
# under qemu user emulation, which isn't a CI fit. BerlekampTest is a
# real algorithmic correctness check (factors a degree-128 polynomial
# over GF(2)) that runs in seconds on every target.
build_only_tests = ['QuickTest', 'ZZTest']
golden_tests     = ['BerlekampTest']

# Wrapper for golden-diff tests: redirect stdin from <Name>In, capture
# stdout, diff against <Name>Out.
golden_runner = files('meson/run-golden-test.sh')

foreach tname : build_only_tests
    executable(tname,
        sources: tname + '.cpp',
        dependencies: ntl_test_dep,
        install: false,
    )
endforeach

foreach tname : golden_tests
    tprog = executable(tname,
        sources: tname + '.cpp',
        dependencies: ntl_test_dep,
        install: false,
    )
    if tests_should_run
        test(tname, golden_runner,
            args: [tprog, files(tname + 'In')[0], files(tname + 'Out')[0]],
            timeout: 1800,
        )
    else
        test(tname + ' (skipped, cross-only target)', python,
             args: ['-c', 'import sys; sys.exit(77)'],
             should_fail: false)
    endif
endforeach
