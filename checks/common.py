"""Reporting for the checks of the derivations: each check prints PASS or
FAIL with its name, and a script exits with status 1 if any check failed."""

import sys

_failed = []


def check(name, ok):
    """Records and prints the outcome of one check."""
    ok = bool(ok)
    print(("PASS  " if ok else "FAIL  ") + name)
    if not ok:
        _failed.append(name)


def finish():
    """Ends a script: exit status 1 if any check failed."""
    if _failed:
        print(f"\n{len(_failed)} check(s) failed")
        sys.exit(1)
    print("\nall checks passed")
