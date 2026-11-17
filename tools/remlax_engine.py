"""Locate the Python engine shipped with the installed remlax R package.

Search order: the REMLAX_PY_DIR variable (directory that contains the remlax/
folder), otherwise system.file("python", package = "remlax") asked of Rscript.
The scripts of this repository import it before `import remlax`.
"""
import os
import subprocess
import sys


def engine_dir():
    d = os.environ.get("REMLAX_PY_DIR")
    if not d:
        d = subprocess.run(["Rscript", "-e", 'cat(system.file("python", package = "remlax"))'],
                           capture_output=True, text=True, check=True).stdout.strip()
    if not d or not os.path.isdir(os.path.join(d, "remlax")):
        raise RuntimeError("remlax engine not found: install the remlax R package "
                           "or set REMLAX_PY_DIR")
    return d


ENGINE_DIR = engine_dir()
if ENGINE_DIR not in sys.path:
    sys.path.insert(0, ENGINE_DIR)
