"""Allow ``python -m igvplot`` to run the CLI."""
from igvplot.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
