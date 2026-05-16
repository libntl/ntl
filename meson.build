# T020: Top-level Meson build for NTL.
#
# Cohabits with the legacy Perl ./configure + Makefile build. No file outside
# the new tree (meson.build, meson.options, src/meson/, src/config.h.in,
# src/MakeDesc.cpp [single small edit], .github/, tools/, ci/, doc/build-meson.txt,
# CHANGELOG.md) is modified by adopting Meson. The Makefile path continues
# to work exactly as before. See doc/build-meson.txt and
# specs/001-meson-cross-compile/ for design context.

project('ntl', 'cpp',
    version: files('version.txt'),
    license: 'LGPL-2.1-or-later',
    meson_version: '>= 1.1.0',
    default_options: [
        'cpp_std=c++11',
        'buildtype=release',
        'warning_level=1',
    ],
)

# -------------------------------------------------------------------
# Resolve the tune-table source.
#
# `-Dtune=default` (the build default):
#     x86 family → 'x86', s390x → 'linux-s390x', else 'generic'.
# `-Dtune=generic|x86|linux-s390x`: use the named static table from
#     src/meson/tune-tables/<name>.ini (FR-005's static fallback).
# `-Dtune=host`: read src/meson/tune-tables/host-tuned.ini (produced
#     by the Python Wizard `ntl-wizard`). Fails if absent.
# `-Dtune=auto`: rejected explicitly with a pointer to ntl-wizard.
# -------------------------------------------------------------------

if get_option('tune') == 'auto'
    error(
        'The legacy auto-tuning Wizard (TUNE=auto) is gone. ' +
        'Run `pip install ./tools && ntl-wizard` to produce ' +
        'src/meson/tune-tables/host-tuned.ini, then build with ' +
        '`-Dtune=host`. For a static table without running the Wizard, ' +
        'pass `-Dtune=generic`, `-Dtune=x86`, or `-Dtune=linux-s390x`.',
    )
endif

# -------------------------------------------------------------------
# Resolve target triplet and load its ABI table entry.
# -------------------------------------------------------------------

cpu_family = host_machine.cpu_family()
os_name    = host_machine.system()
# Best-effort canonical triplet. Cross-files SHOULD set host_machine.cpu /
# system / endian correctly; libc is inferred from suffix below.
triplet = get_option('abi_triplet')
if triplet == ''
    # Derive a default triplet from host_machine.
    if os_name == 'linux'
        libc = 'gnu'  # cross-files override via -Dabi_triplet
        if cpu_family == 'x86'
            triplet = 'i686-linux-' + libc
        elif cpu_family == 'arm'
            triplet = 'armv7l-linux-gnueabihf-musl'
        else
            triplet = cpu_family + '-linux-' + libc
        endif
    elif os_name == 'darwin'
        triplet = cpu_family + '-apple-darwin'
    elif os_name == 'windows'
        triplet = cpu_family + '-w64-mingw32'
    elif os_name == 'freebsd'
        triplet = cpu_family + '-unknown-freebsd'
    else
        error('Cannot derive default triplet for system=' + os_name + '; set -Dabi_triplet=<triplet> explicitly.')
    endif
endif

python = find_program('python3', required: true)
fs = import('fs')
pick_abi = files('src/meson/pick-abi.py')

abi_result = run_command(
    python, pick_abi, '--triplet', triplet,
    check: false,
)
if abi_result.returncode() != 0
    error(
        'pick-abi.py failed for triplet ' + triplet + ':\n' +
        abi_result.stderr(),
    )
endif

# Parse the key=value lines emitted by pick-abi.py into a dict. Our schema
# guarantees values never contain '=', so a plain split suffices.
abi = {}
foreach line : abi_result.stdout().strip().split('\n')
    parts = line.split('=')
    if parts.length() >= 1 and parts[0] != ''
        key = parts[0]
        # Reassemble value (defensive in case the value ever contains '=').
        value = ''
        i = 1
        foreach piece : parts
            if i == parts.length() and i >= 2
                value = value + piece
            elif i >= 2
                value = value + piece + '='
            endif
            i += 1
        endforeach
        abi += { key: value }
    endif
endforeach

# Strict assertion that every required key landed; pick-abi.py validates this
# but defending the Meson side too makes drift detection robust.
foreach key : ['bits_per_long', 'arith_right_shift', 'fma_policy', 'long_double',
               'x86_specializations', 'tune_table', 'exec_mode', 'exe_wrapper',
               'shlib_style', 'threading', 'tls_hack']
    if not abi.has_key(key)
        error('ABI table for ' + triplet + ' is missing key: ' + key)
    endif
endforeach

message('Resolved ABI table for ' + triplet + ': bits_per_long=' +
        abi['bits_per_long'] + ', shlib_style=' + abi['shlib_style'])

# -------------------------------------------------------------------
# Resolve and load the tune-parameter table.
#
# `tune == 'default'` is mapped to one of the static tables based on
# the resolved cpu_family. `tune == 'host'` consumes the artifact
# produced by `ntl-wizard`.
# -------------------------------------------------------------------

tune_opt = get_option('tune')
if tune_opt == 'default'
    if cpu_family == 'x86_64' or cpu_family == 'x86'
        tune_resolved = 'x86'
    elif cpu_family == 's390x'
        tune_resolved = 'linux-s390x'
    else
        tune_resolved = 'generic'
    endif
else
    tune_resolved = tune_opt
endif

if tune_resolved == 'host'
    tune_artifact_path = get_option('tune_artifact')
    if tune_artifact_path == ''
        tune_artifact_path = meson.project_source_root() + '/src/meson/tune-tables/host-tuned.ini'
    endif
    if not fs.exists(tune_artifact_path)
        error(
            'No host-tuned.ini found at ' + tune_artifact_path + '. ' +
            'Run `pip install ' + meson.project_source_root() + '/tools && ntl-wizard` ' +
            'first, or pick a static tune: ' +
            '`-Dtune=default|generic|x86|linux-s390x`.'
        )
    endif
else
    tune_artifact_path = meson.project_source_root() + '/src/meson/tune-tables/' + tune_resolved + '.ini'
endif

read_tune = files('src/meson/read-tune-table.py')
tune_result = run_command(python.full_path(), read_tune, tune_artifact_path, check: false)
if tune_result.returncode() != 0
    error('Failed to read tune table ' + tune_artifact_path + ':\n' + tune_result.stderr())
endif
tune_flags = []
foreach line : tune_result.stdout().strip().split('\n')
    if line != ''
        tune_flags += line
    endif
endforeach
message('Tune table ' + tune_resolved + ' loaded: ' + tune_flags.length().to_string() + ' parameter flags')
add_project_arguments(tune_flags, language: 'cpp')

# -------------------------------------------------------------------
# Dependencies.
# -------------------------------------------------------------------

cpp = meson.get_compiler('cpp')

gmp_opt = get_option('gmp')
gmp_dep = dependency('gmp', required: false)
if not gmp_dep.found() and gmp_opt.allowed()
    # Debian's libgmp-dev does not ship gmp.pc; fall back to manual probing.
    if cpp.has_header('gmp.h')
        gmp_dep = cpp.find_library('gmp', required: gmp_opt.enabled())
    else
        if gmp_opt.enabled()
            error('GMP requested (-Dgmp=enabled) but gmp.h not found')
        endif
        gmp_dep = dependency('', required: false)  # empty placeholder
    endif
endif
use_gmp = gmp_dep.found() and gmp_opt.allowed()

gf2x_opt = get_option('gf2x')
gf2x_dep = dependency('gf2x', required: false)
if not gf2x_dep.found() and gf2x_opt.allowed()
    if cpp.has_header('gf2x.h')
        gf2x_dep = cpp.find_library('gf2x', required: gf2x_opt.enabled())
    else
        if gf2x_opt.enabled()
            error('GF2X requested (-Dgf2x=enabled) but gf2x.h not found')
        endif
        gf2x_dep = dependency('', required: false)
    endif
endif
use_gf2x = gf2x_dep.found() and gf2x_opt.allowed()

threads_dep = dependency('threads', required: get_option('threads'))
use_threads = get_option('threads') and threads_dep.found() and abi['threading'] != 'none'

# -------------------------------------------------------------------
# Build subdirectory.
# -------------------------------------------------------------------

subdir('src')

# -------------------------------------------------------------------
# Summary.
# -------------------------------------------------------------------

summary({
    'target triplet'      : triplet,
    'bits_per_long'       : abi['bits_per_long'],
    'shlib style'         : abi['shlib_style'],
    'long double'         : abi['long_double'],
    'tune table'          : get_option('tune'),
    'threads'             : use_threads,
    'exceptions'          : get_option('exceptions'),
    'GMP'                 : use_gmp,
    'GF2X'                : use_gf2x,
}, section: 'NTL build configuration')
