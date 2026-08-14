#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Docking Vina-GPU para GPU gratis (Saturn Cloud / Paperspace / Colab / Kaggle).

Uso:
  python dock_nube.py <lote.tar.gz> [n_ligands]
Ej:
  python dock_nube.py ligandos_batch6_plano.tar.gz

- Descarga el binario Vina-GPU, los receptores y el lote de ligandos desde GitHub.
- Acopla cada ligando contra TDP43/SOD1/FUS (3 semillas, se queda con la mejor energia).
- Dedup global contra resultados_vinagpu_total.csv de GitHub (no repite pares ya hechos).
- Escribe resultados_vinagpu_<lote>.csv en ~/resultados.
"""
import csv
import glob
import os
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request

REPO = "https://raw.githubusercontent.com/fredy30-Rojas/masive-als-data/main"
VINA_GPU_URL = REPO + "/vinagpu_linux.tar.gz"
RECEPTORES_URL = REPO + "/colab_receptores_plano.tar.gz"

RECEPTORES = {
    "TDP43": {"centro": [28.3, 43.7, 52.5], "tamano": [25, 25, 25]},
    "SOD1": {"centro": [27.9, 111.8, 64.4], "tamano": [25, 25, 25]},
    "FUS": {"centro": [-14.5, 15.1, -7.8], "tamano": [25, 25, 25]},
}
SEEDS = [42, 2026, 777]
THREAD = 8000

WORK = os.path.expanduser("~/masive_als")
RESDIR = os.path.expanduser("~/resultados")


def log(m):
    print(m, flush=True)


def _download(url, dest):
    if os.path.exists(dest):
        os.remove(dest)
    urllib.request.urlretrieve(url, dest)
    return dest


def _buscar_binario(base):
    if not base or not os.path.exists(base):
        return None
    hits = glob.glob(base + "/**/AutoDock-Vina-GPU-2-1", recursive=True)
    return hits[0] if hits else None


def main():
    if len(sys.argv) < 2:
        log("USO: python dock_nube.py <lote.tar.gz> [n_ligands]")
        sys.exit(1)
    lote = sys.argv[1]
    n_ligands = int(sys.argv[2]) if len(sys.argv) > 2 else 0

    os.makedirs(WORK + "/receptores", exist_ok=True)
    os.makedirs(WORK + "/ligandos", exist_ok=True)
    os.makedirs(RESDIR, exist_ok=True)

    # 0) GPU
    log("=== GPU ===")
    g = subprocess.run(["nvidia-smi"], capture_output=True, text=True, timeout=30)
    log("nvidia-smi rc=%d" % g.returncode)
    if g.returncode == 0:
        log(g.stdout.splitlines()[0] if g.stdout else "nvidia-smi ok")

    # Registrar OpenCL ICD de NVIDIA si no esta (algunos contenedores no lo traen)
    try:
        os.makedirs("/etc/OpenCL/vendors", exist_ok=True)
        icd = "/etc/OpenCL/vendors/nvidia.icd"
        if not os.path.exists(icd):
            with open(icd, "w") as f:
                f.write("libnvidia-opencl.so.1\n")
            log("OpenCL ICD registrado")
    except Exception as ex:
        log("AVISO OpenCL ICD: " + str(ex)[:80])

    # 1) Binario Vina-GPU
    BIN_DIR = os.path.expanduser("~/vinagpu_linux")
    VINA = _buscar_binario(BIN_DIR)
    if not VINA:
        os.makedirs(BIN_DIR, exist_ok=True)
        pkg = os.path.expanduser("~/vinagpu_linux.tar.gz")
        _download(VINA_GPU_URL, pkg)
        with tarfile.open(pkg) as t:
            t.extractall(BIN_DIR)
        VINA = _buscar_binario(BIN_DIR)
    if not VINA:
        log("ERROR: no se encontro binario Vina-GPU")
        sys.exit(2)
    os.chmod(VINA, 0o755)
    BIN_RUN_DIR = os.path.dirname(VINA)
    log("Binario Vina-GPU: " + VINA)

    # 2) Receptores
    rpkg = os.path.expanduser("~/receptores.tar.gz")
    _download(RECEPTORES_URL, rpkg)
    with tarfile.open(rpkg) as t:
        t.extractall(WORK + "/receptores")
    for r in RECEPTORES:
        if not os.path.exists(WORK + "/receptores/" + r + ".pdbqt"):
            log("ERROR: falta receptor " + r)
            sys.exit(3)

    # FUS plano (quitar modelos NMR sobrantes)
    fus = WORK + "/receptores/FUS.pdbqt"
    lines = open(fus).read().splitlines()
    ini = next((i for i, l in enumerate(lines) if l.startswith("MODEL")), None)
    fin = next((i for i, l in enumerate(lines) if l.strip() == "ENDMDL"), None)
    keep = lines[ini + 1:fin] if (ini is not None and fin is not None and fin >= ini) else lines
    if not any(l.startswith(("ATOM", "HETATM")) for l in keep):
        keep = lines
    open(fus, "w").write("\n".join(keep) + "\n")

    # 3) Ligandos del lote
    lpkg = os.path.expanduser("~/ligandos.tar.gz")
    _download(REPO + "/" + lote, lpkg)
    with tarfile.open(lpkg) as t:
        t.extractall(WORK + "/ligandos")

    ligs = sorted(glob.glob(WORK + "/ligandos/*.pdbqt"))
    if n_ligands > 0:
        ligs = ligs[:n_ligands]
    log("Ligandos: %d" % len(ligs))
    if not ligs:
        log("ERROR: lote vacio")
        sys.exit(4)

    # 4) CSV local con reanudacion
    base_lote = lote.replace(".tar.gz", "").replace("ligandos_", "").replace("_plano", "")
    CSV = RESDIR + "/resultados_vinagpu_" + base_lote + ".csv"
    if not os.path.exists(CSV):
        with open(CSV, "w", newline="") as f:
            csv.writer(f).writerow(["ligand", "target", "energy", "seed", "timestamp"])
    hechos = {}
    for r in csv.DictReader(open(CSV)):
        hechos.setdefault(r["ligand"], set()).add(r["target"])

    # Dedup global (GitHub)
    try:
        urllib.request.urlretrieve(REPO + "/resultados_vinagpu_total.csv", "/tmp/total_prev.csv")
        n_glob = 0
        for r in csv.DictReader(open("/tmp/total_prev.csv")):
            hechos.setdefault(r["ligand"], set()).add(r["target"])
            n_glob += 1
        log("Dedup global: %d pares ya hechos" % n_glob)
    except Exception as ex:
        log("AVISO sin dedup global: " + str(ex)[:60])

    def acoplar(lig, target, seed):
        info = RECEPTORES[target]
        cfg = "/tmp/cfg.txt"
        with open(cfg, "w") as f:
            f.write("receptor = %s\n" % (WORK + "/receptores/" + target + ".pdbqt"))
            f.write("ligand = %s\n" % lig)
            f.write("center_x = %s\n" % info["centro"][0])
            f.write("center_y = %s\n" % info["centro"][1])
            f.write("center_z = %s\n" % info["centro"][2])
            f.write("size_x = %s\n" % info["tamano"][0])
            f.write("size_y = %s\n" % info["tamano"][1])
            f.write("size_z = %s\n" % info["tamano"][2])
            f.write("num_modes = 3\n")
            f.write("seed = %d\n" % seed)
            f.write("thread = %d\n" % THREAD)
        try:
            r = subprocess.run([VINA, "--config", cfg], capture_output=True, text=True,
                               timeout=1800, cwd=BIN_RUN_DIR)
            out = (r.stdout or "") + (r.stderr or "")
            if r.returncode != 0:
                return None, "rc=%d %s" % (r.returncode, out[-200:])
            for ln in out.splitlines():
                s = ln.split()
                if len(s) >= 2 and s[0] == "1":
                    try:
                        return round(float(s[1]), 4), None
                    except ValueError:
                        pass
            return None, "sin afinidad: " + out[-150:]
        except Exception as ex:
            return None, str(ex)[:120]

    def tiene_atomos(lig):
        try:
            with open(lig, errors="replace") as f:
                return any(l.startswith(("ATOM", "HETATM")) for l in f)
        except Exception:
            return False

    pendientes = []
    for lig in ligs:
        nombre = os.path.basename(lig).replace(".pdbqt", "")
        hechos_lig = hechos.get(nombre, set())
        for target in RECEPTORES:
            if target not in hechos_lig:
                pendientes.append((lig, nombre, target))
    log("Pares pendientes: %d" % len(pendientes))

    t0 = time.time()
    n_ok = 0
    n_err = 0
    for lig, nombre, target in pendientes:
        if not tiene_atomos(lig):
            continue
        energias = []
        for sd in SEEDS:
            e, er = acoplar(lig, target, sd)
            if e is not None:
                energias.append(e)
            else:
                n_err += 1
                if n_err <= 5:
                    log("ERROR %s %s seed=%s %s" % (nombre, target, sd, er))
        if energias:
            energia = min(energias)
            with open(CSV, "a", newline="") as f:
                csv.writer(f).writerow([nombre, target, energia, "best",
                                        time.strftime("%Y-%m-%d %H:%M:%S")])
            n_ok += 1
        if n_ok % 10 == 0 and n_ok > 0:
            log("[%d acoplados] %.1f min (errores=%d)" % (n_ok, (time.time() - t0) / 60, n_err))

    log("=== LOTE COMPLETADO ===")
    log("Filas en CSV: %d" % (n_ok + (sum(1 for _ in open(CSV)) - 1)))
    log("CSV guardado en: " + CSV)


if __name__ == "__main__":
    main()
