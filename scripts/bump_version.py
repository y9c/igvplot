#!/usr/bin/env python3
"""Bump the igvplot version across pyproject.toml and igvplot/__init__.py.

Usage::

    python scripts/bump_version.py --patch   # 0.0.1 -> 0.0.2
    python scripts/bump_version.py --minor   # 0.0.1 -> 0.1.0
    python scripts/bump_version.py --major   # 0.0.1 -> 1.0.0
    python scripts/bump_version.py --set 0.2.0
    python scripts/bump_version.py --patch --dry-run      # show, don't write
    python scripts/bump_version.py --patch --tag          # bump, commit, tag v0.1.1

Keeps the version in one source of truth (pyproject.toml) and mirrors it into
``igvplot/__init__.py``. With ``--tag`` it commits both files and creates an
annotated ``v<new>`` tag (which triggers the .github/workflows/release.yml CI).
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = ROOT / "pyproject.toml"
INIT = ROOT / "igvplot" / "__init__.py"

VERSION_RE = re.compile(r'(?m)^(\s*version\s*=\s*["\'])([^"\']+)(["\'])')
INIT_RE = re.compile(r'(?m)^(__version__\s*=\s*["\'])([^"\']+)(["\'])')


def current_version() -> str:
    m = VERSION_RE.search(PYPROJECT.read_text())
    if not m:
        raise SystemExit("no `version = \"...\"` found in pyproject.toml")
    return m.group(2).strip()


def _parts(v: str):
    return [int(x) if x.isdigit() else 0 for x in v.split(".")]


def bump(v: str, which: str) -> str:
    p = _parts(v)
    p = p + [0, 0, 0]
    if which == "major":
        p = [p[0] + 1, 0, 0]
    elif which == "minor":
        p = [p[0], p[1] + 1, 0]
    else:
        p = [p[0], p[1], p[2] + 1]
    return ".".join(str(x) for x in p[:3])


def _sub(text: str, regex, new: str) -> str:
    return regex.sub(lambda m: m.group(1) + new + m.group(3), text, count=1)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--patch", action="store_true", help="bump patch (0.0.1 -> 0.0.2)")
    g.add_argument("--minor", action="store_true", help="bump minor (0.0.1 -> 0.1.0)")
    g.add_argument("--major", action="store_true", help="bump major (0.0.1 -> 1.0.0)")
    g.add_argument("--set", metavar="X.Y.Z", help="set an exact version")
    ap.add_argument("--dry-run", action="store_true", help="print changes without writing")
    ap.add_argument("--tag", action="store_true", help="commit both files and tag v<new>")
    args = ap.parse_args(argv)

    cur = current_version()
    if args.set:
        new = args.set
    elif args.patch or args.minor or args.major:
        which = "major" if args.major else ("minor" if args.minor else "patch")
        new = bump(cur, which)
    else:
        ap.print_help()
        return 1

    if not re.fullmatch(r"\d+\.\d+\.\d+", new):
        raise SystemExit(f"bad version string {new!r}; expected MAJOR.MINOR.PATCH")

    pytext = PYPROJECT.read_text()
    inittext = INIT.read_text()

    if args.dry_run:
        print(f"pyproject.toml          {cur} -> {new}")
        print(f"igvplot/__init__.py     {cur} -> {new}")
        print("(dry run; nothing written)")
        return 0

    PYPROJECT.write_text(_sub(pytext, VERSION_RE, new))
    INIT.write_text(_sub(inittext, INIT_RE, new))
    print(f"bumped {cur} -> {new}")

    if args.tag:
        subprocess.run(["git", "add", "pyproject.toml", "igvplot/__init__.py"], check=True)
        subprocess.run(["git", "commit", "-m", f"Release {new}"], check=True)
        subprocess.run(["git", "tag", "-a", f"v{new}", "-m", f"igvplot {new}"], check=True)
        print(f"committed and tagged v{new}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
