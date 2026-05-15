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

# The Wizard (TUNE=auto) is intentionally unsupported by the Meson path; see
# FR-014. The tune option's combo type already rejects 'auto', but we also
# guard explicitly for clarity in error messages.
if get_option('tune') == 'auto'
    error(
        'The auto-tuning Wizard is not supported under the Meson build. ' +
        'Use the Makefile build (./configure TUNE=auto) if you need ' +
        'Wizard-tuned parameters. The Meson build accepts only the static ' +
        'tune tables: generic, x86, linux-s390x.',
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
