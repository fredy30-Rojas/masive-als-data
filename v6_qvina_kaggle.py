# CELDA 5b: Docking con QVina 2.1 en KAGGLE (dataset + opcional internet)
# Adaptado de la celda 5 del notebook Kaggle (Vina) + v6_qvina.py (Colab).
# - Mismo formato CSV: ligand,target,energy,timestamp
# - Subconjuntos pares/impares, PDBQT reparados del dataset, reanudacion por CSV
# - VALIDAR=True: acopla una muestra de pares YA hechos (seed 42) y compara con resultados_colab.csv
# - QVina2 es CPU multihilo (rapido); Vina-GPU seria para GPU real (requiere compilar).
# - Con internet: descarga binario desde GitHub. Sin internet: lo toma del dataset (qvina_kaggle.tar.gz).
import csv, glob, os, shutil, subprocess, tarfile, time, urllib.request

WORK = '/kaggle/working/masive_als'
IN = '/kaggle/input/masive-als-datos'
QVINA_DIR = '/kaggle/working/qvina_linux'
QVINA_BIN = QVINA_DIR + '/qvina2'
QVINA_URL = 'https://raw.githubusercontent.com/fredy30-Rojas/masive-als-data/main/qvina_kaggle.tar.gz'

# --- 0) Verificar GPU (informativo: QVina2 es CPU; si hay GPU y se quiere Vina-GPU, otro script) ---
print('=== VERIFICACION GPU/CPU ===', flush=True)
gpu = subprocess.run(['nvidia-smi'], capture_output=True, text=True, timeout=30)
print('nvidia-smi rc=%d (GPU disponible si rc=0)' % gpu.returncode, flush=True)
print('CPUs:', os.cpu_count(), flush=True)

# --- 1) Obtener binario QVina + librerias Boost ---
# Con internet: descarga el paquete completo desde GitHub. Sin internet: desde el dataset.
def obtener_qvina():
    os.makedirs(QVINA_DIR, exist_ok=True)
    internet = True
    try:
        urllib.request.urlopen('https://raw.githubusercontent.com', timeout=10).close()
    except Exception:
        internet = False
    print('Internet:', 'SI' if internet else 'NO', flush=True)
    if internet:
        try:
            pkg = QVINA_DIR + '/qvina_kaggle.tar.gz'
            urllib.request.urlretrieve(QVINA_URL, pkg)
            with tarfile.open(pkg) as t:
                t.extractall(QVINA_DIR)
            return True
        except Exception as ex:
            print('Descarga GitHub fallo:', str(ex)[:100], flush=True)
    # Sin internet (o fallo): usar el dataset
    src_tar = IN + '/qvina_kaggle.tar.gz'
    if os.path.exists(src_tar):
        with tarfile.open(src_tar) as t:
            t.extractall(QVINA_DIR)
        return True
    # Ultimo fallback: archivos planos dentro del dataset
    ok = False
    for f in ['qvina2', 'libboost_filesystem.so.1', 'libboost_program_options.so.1', 'libboost_thread.so.1']:
        s = IN + '/' + f
        if os.path.exists(s):
            shutil.copy(s, QVINA_DIR + '/' + f)
            ok = True
    return ok

if not os.path.exists(QVINA_BIN):
    if not obtener_qvina():
        print('ERROR: no se pudo obtener el binario QVina. Avisa a Fredy.', flush=True)
        raise SystemExit(1)
    os.chmod(QVINA_BIN, 0o755)
# LD_LIBRARY_PATH apuntando a las libs boost del paquete
os.environ['LD_LIBRARY_PATH'] = QVINA_DIR + ':' + os.environ.get('LD_LIBRARY_PATH', '')
print('Binario QVina listo:', QVINA_BIN, flush=True)
print(subprocess.run([QVINA_BIN, '--help'], capture_output=True, text=True, timeout=30).stdout[:400], flush=True)

# --- 2) Datos frescos desde /kaggle/input ---
for sub in ['receptores', 'ligandos', 'resultados']:
    os.makedirs(WORK + '/' + sub, exist_ok=True)

def es_receptor(n):
    return n in ('TDP43.pdbqt', 'SOD1.pdbqt', 'FUS.pdbqt')

for d in ('receptores', 'ligandos'):
    for f in glob.glob(WORK + '/' + d + '/*.pdbqt'):
        os.remove(f)
for raiz, _, archivos in os.walk(IN):
    for a in archivos:
        if a.endswith('.pdbqt'):
            destino = WORK + '/receptores/' if es_receptor(a) else WORK + '/ligandos/'
            shutil.copy(os.path.join(raiz, a), destino + a)

# --- 2b) PDBQT reparados (del dataset, sobreescriben los corruptos) ---
REPARADOS = ['CHEMBL1076399', 'CHEMBL1163427', 'CHEMBL1203109', 'CHEMBL1203132',
             'CHEMBL1203140', 'CHEMBL1203155', 'CHEMBL1203199', 'CHEMBL1203224',
             'CHEMBL1203252', 'CHEMBL1204421', 'CHEMBL1207772', 'CHEMBL1208195',
             'CHEMBL152893']
n_cop = 0
for src in sorted(glob.glob(IN + '/ligandos_reparados/**/*.pdbqt', recursive=True)):
    nombre = os.path.basename(src).replace('.pdbqt', '')
    if nombre in REPARADOS:
        shutil.copy(src, WORK + '/ligandos/' + nombre + '.pdbqt')
        n_cop += 1
print('PDBQT reparados copiados:', n_cop, flush=True)
print('Receptores:', len(glob.glob(WORK + '/receptores/*.pdbqt')),
      '| Ligandos:', len(glob.glob(WORK + '/ligandos/*.pdbqt')), flush=True)

# --- 3) Receptores ---
RECEPTORES = {
    'TDP43': {'archivo': WORK + '/receptores/TDP43.pdbqt', 'centro': [28.3, 43.7, 52.5], 'tamano': [25, 25, 25]},
    'SOD1': {'archivo': WORK + '/receptores/SOD1.pdbqt', 'centro': [27.9, 111.8, 64.4], 'tamano': [25, 25, 25]},
    'FUS': {'archivo': WORK + '/receptores/FUS.pdbqt', 'centro': [-14.5, 15.1, -7.8], 'tamano': [25, 25, 25]},
}

# --- 4) FUS plano (sin MODEL/ENDMDL) ---
fus = RECEPTORES['FUS']['archivo']
try:
    lines = open(fus).read().splitlines()
    ini = next((i for i, l in enumerate(lines) if l.startswith('MODEL')), None)
    fin = next((i for i, l in enumerate(lines) if l.strip() == 'ENDMDL'), None)
    keep = lines[ini + 1:fin] if (ini is not None and fin is not None and fin >= ini) else lines
    if not any(l.startswith(('ATOM', 'HETATM')) for l in keep):
        keep = lines
    open(fus, 'w').write('\n'.join(keep) + '\n')
    print('FUS plano: ATOM=%d' % sum(1 for l in keep if l.startswith(('ATOM', 'HETATM'))), flush=True)
except Exception as ex:
    print('ERROR reparando FUS:', str(ex)[:100], flush=True)

# --- 5) CSV (reanuda) ---
CSV = WORK + '/resultados/resultados_qvina_kaggle.csv'
if not os.path.exists(CSV):
    with open(CSV, 'w', newline='') as f:
        csv.writer(f).writerow(['ligand', 'target', 'energy', 'timestamp'])
SALTADOS = WORK + '/resultados/saltados_qvina.txt'
saltados = set()
if os.path.exists(SALTADOS):
    saltados = set(l.strip() for l in open(SALTADOS) if l.strip())
hechos = {}
for r in csv.DictReader(open(CSV)):
    hechos.setdefault(r['ligand'], set()).add(r['target'])

# --- 6) Base Colab (pares ya hechos) ---
BASE = WORK + '/resultados/resultados_colab.csv'
try:
    shutil.copy(IN + '/resultados_colab.csv', BASE)
    for r in csv.DictReader(open(BASE)):
        hechos.setdefault(r['ligand'], set()).add(r['target'])
    print('Base Colab fusionada: %d pares ya hechos' % sum(len(v) for v in hechos.values()), flush=True)
except Exception as ex:
    print('AVISO sin base Colab:', str(ex)[:80], flush=True)

# --- 7) Seleccion de pares (subconjunto + pendientes) ---
SUBCONJUNTO = 'pares'   # 'pares', 'impares' o 'todos'
VALIDAR = True          # True = primero validacion (10 ligandos x 3); False = tanda completa

ligs = sorted(glob.glob(WORK + '/ligandos/*.pdbqt'))
if SUBCONJUNTO in ('pares', 'impares'):
    paridad = 0 if SUBCONJUNTO == 'pares' else 1
    ligs = [l for i, l in enumerate(ligs) if i % 2 == paridad]
print('Ligandos (subconjunto=%s):' % SUBCONJUNTO, len(ligs), flush=True)

pendientes = []
for lig in ligs:
    nombre = os.path.basename(lig).replace('.pdbqt', '')
    if nombre in saltados:
        continue
    para = hechos.get(nombre, set())
    for target in RECEPTORES:
        if target not in para:
            pendientes.append((lig, nombre, target))
print('Pares pendientes:', len(pendientes), flush=True)

# --- Validacion: limitar a muestra de ligandos con las 3 proteinas hechas ---
if VALIDAR and pendientes:
    completos = {n for n in hechos if len(hechos[n]) >= 3}
    muestras = [p for p in pendientes if p[1] in completos][:30]  # 10 ligandos x 3
    if muestras:
        print('VALIDACION: re-acoplo %d pares ya hechos' % len(muestras), flush=True)
        pendientes = muestras

# --- 8) Acoplar con QVina ---
def acoplar(lig, target, exhaustiveness=2):
    info = RECEPTORES[target]
    base = [QVINA_BIN,
            '--receptor', info['archivo'],
            '--ligand', lig,
            '--center_x', str(info['centro'][0]),
            '--center_y', str(info['centro'][1]),
            '--center_z', str(info['centro'][2]),
            '--size_x', str(info['tamano'][0]),
            '--size_y', str(info['tamano'][1]),
            '--size_z', str(info['tamano'][2]),
            '--exhaustiveness', str(exhaustiveness),
            '--num_modes', '3', '--seed', '42']
    variantes = [base + ['--threads', '2'], base]  # con threads primero, fallback sin el
    try:
        r = None
        for cmd in variantes:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
            out = (r.stdout or '') + (r.stderr or '')
            if r.returncode == 0:
                break
        if r is None or r.returncode != 0:
            return None, 'qvina rc=%d %s' % (r.returncode, out[-150:])
        for ln in out.splitlines():
            s = ln.split()
            if len(s) >= 2 and s[0] == '1':
                try:
                    return round(float(s[1]), 4), None
                except ValueError:
                    pass
        return None, 'sin afinidad: ' + out[-120:]
    except Exception as ex:
        return None, str(ex)[:100]

def tiene_atomos(lig):
    try:
        with open(lig, errors='replace') as f:
            return any(l.startswith(('ATOM', 'HETATM')) for l in f)
    except Exception:
        return False

if pendientes:
    t0 = time.time()
    n_ok = 0
    for lig, nombre, target in pendientes:
        if not tiene_atomos(lig):
            print('LIGANDO_VACIO', nombre, '- se salta', flush=True)
            with open(SALTADOS, 'a') as f:
                f.write(nombre + '\n')
            continue
        energia, err = acoplar(lig, target)
        if energia is not None:
            with open(CSV, 'a', newline='') as f:
                csv.writer(f).writerow([nombre, target, energia, time.strftime('%Y-%m-%d %H:%M:%S')])
            n_ok += 1
        else:
            print('ERROR', nombre, target, err, flush=True)
        if n_ok % 5 == 0 and n_ok > 0:
            print('[%d acoplados] %.1f min' % (n_ok, (time.time() - t0) / 60), flush=True)
    print('Acoplados ahora:', n_ok, flush=True)

# --- 9) Comparacion QVina vs Vina 1.2.5 ---
if VALIDAR:
    print()
    print('=== COMPARACION QVina vs Vina 1.2.5 (seed 42) ===', flush=True)
    refs = {}
    for r in csv.DictReader(open(BASE)):
        refs.setdefault(r['ligand'], {})[r['target']] = float(r['energy'])
    diffs = []
    for r in csv.DictReader(open(CSV)):
        q = float(r['energy'])
        ref = refs.get(r['ligand'], {}).get(r['target'])
        if ref is not None:
            diffs.append((r['ligand'], r['target'], ref, q, abs(q - ref)))
    if diffs:
        for lig, tgt, ref, q, d in sorted(diffs, key=lambda x: -x[4])[:15]:
            print('%-22s %-5s Vina=%7.2f QVina=%7.2f diff=%5.2f' % (lig, tgt, ref, q, d), flush=True)
        mean = sum(d[4] for d in diffs) / len(diffs)
        maxd = max(d[4] for d in diffs)
        print('MEDIA diff=%.3f | MAX diff=%.3f | n=%d' % (mean, maxd, len(diffs)), flush=True)
        if maxd > 1.5:
            print('ALERTA: diferencias grandes - NO mezclar motores sin revisar.', flush=True)
        else:
            print('OK: diferencias pequenas, motores compatibles.', flush=True)

print(flush=True)
print('=== TANDAS COMPLETADAS ===', flush=True)
print('Filas en CSV:', len(list(csv.DictReader(open(CSV)))), flush=True)
