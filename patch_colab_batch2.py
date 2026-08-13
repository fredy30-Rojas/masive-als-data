#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Anade ligandos_batch2_plano.tar.gz a la lista URLS del notebook de Colab"""
import json, io, sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

F = "MASIVE-ALS_Colab_GPU_VinaGPU.ipynb"
nb = json.load(open(F, encoding="utf-8"))

# La lista URLS esta en la celda 3 (docking). Buscar la linea fda_pdbqt_plano y anadir batch2 despues.
changed = False
for c in nb["cells"]:
    src = "".join(c.get("source", []))
    if "fda_pdbqt_plano.tar.gz" in src and "ligandos_batch2_plano.tar.gz" not in src:
        lines = src.split("\n")
        out = []
        for ln in lines:
            out.append(ln)
            if "fda_pdbqt_plano.tar.gz" in ln and not changed:
                out.append("    'https://raw.githubusercontent.com/fredy30-Rojas/masive-als-data/main/ligandos_batch2_plano.tar.gz',")
                changed = True
        c["source"] = ["%s\n" % l for l in out]

print("batch2 anadido a URLS:", changed)
if changed:
    json.dump(nb, open(F, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("guardado OK")
else:
    print("AVISO: no se encontro fda_pdbqt o ya tenia batch2")
