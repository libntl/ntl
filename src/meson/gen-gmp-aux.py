#!/usr/bin/env python3
"""Emit a minimal NTL gmp_aux.h to stdout.

Used as a custom_target generator from src/NTL/meson.build. The sole varying
value is NTL_GMP_LIMB_T_SIZE_BITS, passed as the only argument.
"""

import sys


def main() -> int:
    if len(sys.argv) != 2:
        print(
            "usage: gen-gmp-aux.py <limb-size-bits>", file=sys.stderr
        )
        return 2
    limb_bits = int(sys.argv[1])
    sys.stdout.write(
        "#ifndef NTL_gmp_aux__H\n"
        "#define NTL_gmp_aux__H\n"
        f"#define NTL_GMP_LIMB_T_SIZE_BITS {limb_bits}\n"
        "#endif\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
