"""Enable `python -m ntl_wizard` invocation."""
from .cli import main

if __name__ == "__main__":  # pragma: no cover
    import sys
    sys.exit(main())
