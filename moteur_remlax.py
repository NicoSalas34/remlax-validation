"""Localise le moteur Python livre avec le paquet R remlax installe.

Ordre : variable REMLAX_PY_DIR (repertoire qui contient le dossier remlax/),
sinon system.file("python", package = "remlax") demande a Rscript.
Les scripts de ce depot l'importent avant `import remlax`.
"""
import os
import subprocess
import sys


def repertoire_moteur():
    d = os.environ.get("REMLAX_PY_DIR")
    if not d:
        d = subprocess.run(["Rscript", "-e", 'cat(system.file("python", package = "remlax"))'],
                           capture_output=True, text=True, check=True).stdout.strip()
    if not d or not os.path.isdir(os.path.join(d, "remlax")):
        raise RuntimeError("moteur remlax introuvable : installer le paquet R remlax "
                           "ou definir REMLAX_PY_DIR")
    return d


MOTEUR = repertoire_moteur()
if MOTEUR not in sys.path:
    sys.path.insert(0, MOTEUR)
