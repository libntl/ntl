#!/usr/bin/env python3
"""Emit NTL gmp_aux.h to stdout.

Replaces the role of src/gen_gmp_aux.cpp for the Meson build. The C++
program executes at build time, which means under cross-compile it
either won't run at all or (worse) runs on the build host with the
host's GMP and aborts when its consistency checks see a mismatch
against the target's expected bits-per-long.

This script consumes two values that are known to the build system at
configure time:

  bits_per_limb : size of mp_limb_t for the *target*, in bits, from
                  cc.sizeof('mp_limb_t', prefix: '#include <gmp.h>').
                  Compile-time, works in cross mode.
  bits_per_long : the target's bits-per-long, from the ABI table.

and emits the same set of macros gen_gmp_aux.cpp would have written:

  NTL_ZZ_NBITS, NTL_BITS_PER_LIMB_T, NTL_ZZ_FRADIX, and optionally
  NTL_SMALL_MP_SIZE_T (only when sizeof(mp_size_t) < sizeof(long),
  which is a 32-on-64 oddity; we conservatively omit it).

The output matches gen_gmp_aux.cpp's output byte-for-byte on the
mainstream case (mp_bits_per_limb == bits_per_long, nail_bits == 0).

Usage:
    gen-gmp-aux.py <bits_per_limb> <bits_per_long>
"""

from __future__ import annotations

import sys


def print2k(k: int, bpl: int) -> str:
    """Express 2^k as a product of `((double)(1L<<l))` factors with l < bpl.

    Mirrors src/gen_gmp_aux.cpp's print2k() byte-for-byte (in the bpl > 0
    case). When k == 0, gen_gmp_aux.cpp emits the literal `((double) 1.0)`.
    """
    if k <= 0:
        return "((double) 1.0)"

    m = bpl - 2
    pieces: list[str] = []
    while k > 0:
        l = m if k > m else k
        k -= l
        pieces.append(f"((double)(1L<<{l}))")
    return "(" + "*".join(pieces) + ")"


def main() -> int:
    if len(sys.argv) != 3:
        print(
            "usage: gen-gmp-aux.py <bits_per_limb> <bits_per_long>",
            file=sys.stderr,
        )
        return 2
    bits_per_limb = int(sys.argv[1])
    bits_per_long = int(sys.argv[2])

    # Sanity check the same way gen_gmp_aux.cpp does (less strictly —
    # we can't query GMP_NAIL_BITS from Python without compiling).
    if bits_per_limb not in (bits_per_long, 2 * bits_per_long):
        print(
            f"WARNING: bits_per_limb ({bits_per_limb}) is not "
            f"bits_per_long ({bits_per_long}) or 2x that. The target's "
            f"GMP may behave unexpectedly.",
            file=sys.stderr,
        )
    ntl_zz_nbits = bits_per_limb  # assumes nail_bits == 0

    out = []
    # gen_gmp_aux.cpp does not wrap its output in #ifndef guards. The
    # file is included transitively but always at the same depth, and
    # NTL has its own guards elsewhere. We preserve that behavior.
    out.append(f"#define NTL_ZZ_NBITS ({ntl_zz_nbits})\n")
    out.append(f"#define NTL_BITS_PER_LIMB_T ({bits_per_limb})\n")
    out.append(f"#define NTL_ZZ_FRADIX {print2k(ntl_zz_nbits, bits_per_long)}\n")

    sys.stdout.write("".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
