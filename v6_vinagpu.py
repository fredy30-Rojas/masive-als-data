# CELDA 5: Docking con Vina-GPU-2.1 (usa la GPU de verdad) - Colab y Kaggle
# - Mismo formato CSV: ligand,target,energy,timestamp
# - Subconjuntos pares/impares, reanudacion por CSV
# - VALIDAR=True: acopla una muestra de pares YA hechos (seed 42) y compara con resultados_colab.csv
# - NOTA: Vina-GPU-2.1 necesita GPU NVIDIA + OpenCL (CUDA toolkit). Colab T4 y Kaggle GPU: OK.
import csv, glob, os, shutil, subprocess, time, urllib.request

# ===== CONFIGURACION =====
EN_KAGGLE = '/kaggle' in os.getcwd() or os.path.exists('/kaggle/working')
if EN_KAGGLE:
    WORK = '/kaggle/working/masive_als'
    IN = '/kaggle/input/masive-als-datos'
    BIN_DIR = '/kaggle/working/vinagpu_linux'
else:
    WORK = '/content/masive_als'
    IN = None
    BIN_DIR = '/content/vinagpu_linux'

VINA_GPU_BIN = BIN_DIR + '/AutoDock-Vina-GPU-2-1'
VINA_GPU_URL = 'https://raw.githubusercontent.com/fredy30-Rojas/masive-als-data/main/vinagpu_linux.tar.gz'
SUBCONJUNTO = 'pares'   # 'pares', 'impares' o 'todos'
VALIDAR = True          # True = primero validacion (10 ligandos x 3); False = tanda completa

# --- 0) Verificar GPU (REQUISITO: Vina-GPU sin GPU no funciona) ---
print('=== VERIFICACION GPU ===', flush=True)
gpu = subprocess.run(['nvidia-smi'], capture_output=True, text=True, timeout=30)
print('nvidia-smi rc=%d' % gpu.returncode, flush=True)
if gpu.returncode != 0:
    print('ERROR: Vina-GPU necesita GPU NVIDIA. Activa GPU T4/P100 en el runtime.', flush=True)
    print('En Colab: Entorno de ejecucion > Cambiar tipo > T4 GPU.', flush=True)
    print('En Kaggle: Settings > Accelerator > GPU T4 x2 o P100.', flush=True)
    raise SystemExit(1)
print([l for l in gpu.stdout.splitlines() if 'Tesla' in l or 'NVIDIA' in l or 'P100' in l or 'T4' in l][:2], flush=True)

# --- 1) Obtener binario Vina-GPU ---
if not os.path.exists(VINA_GPU_BIN):
    os.makedirs(BIN_DIR, exist_ok=True)
    pkg = '/tmp/vinagpu_linux.tar.gz'
    # Con internet: descarga desde GitHub
    internet = True
    try:
        urllib.request.urlopen('https://raw.githubusercontent.com', timeout=10).close()
    except Exception:
        internet = False
    if internet:
        try:
            urllib.request.urlretrieve(VINA_GPU_URL, pkg)
        except Exception as ex:
            print('Descarga GitHub fallo:', str(ex)[:100], flush=True)
            internet = False
    if not internet and EN_KAGGLE:
        # Sin internet: usar el dataset (si Fredy lo sube)
        src = IN + '/vinagpu_linux.tar.gz'
        if os.path.exists(src):
            shutil.copy(src, pkg)
    if not os.path.exists(pkg):
        print('ERROR: no se pudo obtener el binario Vina-GPU.', flush=True)
        raise SystemExit(1)
    import tarfile
    with tarfile.open(pkg) as t:
        t.extractall(BIN_DIR)
    os.chmod(VINA_GPU_BIN, 0o755)
print('Binario Vina-GPU listo:', VINA_GPU_BIN, flush=True)

# --- 2) Datos (Kaggle: dataset; Colab: GitHub) ---
for sub in ['receptores', 'ligandos', 'resultados']:
    os.makedirs(WORK + '/' + sub, exist_ok=True)

def es_receptor(n):
    return n in ('TDP43.pdbqt', 'SOD1.pdbqt', 'FUS.pdbqt')

if EN_KAGGLE:
    for d in ('receptores', 'ligandos'):
        for f in glob.glob(WORK + '/' + d + '/*.pdbqt'):
            os.remove(f)
    for raiz, _, archivos in os.walk(IN):
        for a in archivos:
            if a.endswith('.pdbqt'):
                destino = WORK + '/receptores/' if es_receptor(a) else WORK + '/ligandos/'
                shutil.copy(os.path.join(raiz, a), destino + a)
    # PDBQT reparados
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
else:
    import tarfile
    URLS = [
        'https://raw.githubusercontent.com/fredy30-Rojas/masive-als-data/main/colab_receptores_plano.tar.gz',
        'https://raw.githubusercontent.com/fredy30-Rojas/masive-als-data/main/colab_ligandos50_plano.tar.gz',
    ]
    extract_dir = '/content/tmp_extract'
    if os.path.exists(extract_dir):
        shutil.rmtree(extract_dir)
    os.makedirs(extract_dir)
    for url in URLS:
        dest = '/content/' + os.path.basename(url)
        if os.path.exists(dest):
            os.remove(dest)
        urllib.request.urlretrieve(url, dest)
    for tar in sorted(glob.glob('/content/*.tar.gz')):
        with tarfile.open(tar) as t:
            try:
                t.extractall(extract_dir, filter='data')
            except TypeError:
                t.extractall(extract_dir)
    for raiz, _, archivos in os.walk(extract_dir):
        for a in archivos:
            if a.endswith('.pdbqt'):
                destino = WORK + '/receptores/' if es_receptor(a) else WORK + '/ligandos/'
                shutil.copy(os.path.join(raiz, a), destino + a)

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
CSV = WORK + '/resultados/resultados_vinagpu.csv'
if not os.path.exists(CSV):
    with open(CSV, 'w', newline='') as f:
        csv.writer(f).writerow(['ligand', 'target', 'energy', 'timestamp'])
SALTADOS = WORK + '/resultados/saltados_vinagpu.txt'
saltados = set()
if os.path.exists(SALTADOS):
    saltados = set(l.strip() for l in open(SALTADOS) if l.strip())
hechos = {}
for r in csv.DictReader(open(CSV)):
    hechos.setdefault(r['ligand'], set()).add(r['target'])

# --- 6) Base Colab (pares ya hechos) ---
BASE = WORK + '/resultados/resultados_colab.csv'
try:
    if EN_KAGGLE:
        shutil.copy(IN + '/resultados_colab.csv', BASE)
    else:
        urllib.request.urlretrieve(
            'https://raw.githubusercontent.com/fredy30-Rojas/masive-als-data/main/resultados_colab.csv', BASE)
    for r in csv.DictReader(open(BASE)):
        hechos.setdefault(r['ligand'], set()).add(r['target'])
    print('Base Colab fusionada: %d pares ya hechos' % sum(len(v) for v in hechos.values()), flush=True)
except Exception as ex:
    print('AVISO sin base Colab:', str(ex)[:80], flush=True)

# --- 7) Seleccion de pares ---
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

if VALIDAR and pendientes:
    completos = {n for n in hechos if len(hechos[n]) >= 3}
    muestras = [p for p in pendientes if p[1] in completos][:30]
    if muestras:
        print('VALIDACION: re-acoplo %d pares ya hechos' % len(muestras), flush=True)
        pendientes = muestras

# --- 8) Acoplar con Vina-GPU (config file por par) ---
# NOTA: Vina-GPU-2.1 NO acepta --exhaustiveness ni --cpu (comentados en su parser).
# El control de esfuerzo es --thread (numero de tareas de computo en la GPU).
# Ejecutamos con cwd=BIN_DIR para que encuentre los kernels ./OpenCL/ (default_work_path=".").
def acoplar(lig, target, thread=32):
    info = RECEPTORES[target]
    cfg = '/tmp/cfg_%s_%s.txt' % (os.path.basename(lig).replace('.pdbqt', '')[:20], target)
    with open(cfg, 'w') as f:
        f.write('receptor = %s\n' % info['archivo'])
        f.write('ligand = %s\n' % lig)
        f.write('center_x = %s\n' % info['centro'][0])
        f.write('center_y = %s\n' % info['centro'][1])
        f.write('center_z = %s\n' % info['centro'][2])
        f.write('size_x = %s\n' % info['tamano'][0])
        f.write('size_y = %s\n' % info['tamano'][1])
        f.write('size_z = %s\n' % info['tamano'][2])
        f.write('num_modes = 3\n')
        f.write('seed = 42\n')
        f.write('thread = %d\n' % thread)
    try:
        r = subprocess.run([VINA_GPU_BIN, '--config', cfg], capture_output=True, text=True,
                           timeout=1800, cwd=BIN_DIR)
        out = (r.stdout or '') + (r.stderr or '')
        if r.returncode != 0:
            return None, 'vinagpu rc=%d %s' % (r.returncode, out[-200:])
        for ln in out.splitlines():
            s = ln.split()
            if len(s) >= 2 and s[0] == '1':
                try:
                    return round(float(s[1]), 4), None
                except ValueError:
                    pass
        return None, 'sin afinidad: ' + out[-150:]
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

# --- 9) Comparacion con Vina 1.2.5 ---
if VALIDAR and os.path.exists(BASE):
    print()
    print('=== COMPARACION Vina-GPU vs Vina 1.2.5 (seed 42) ===', flush=True)
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
            print('%-22s %-5s Vina=%7.2f GPU=%7.2f diff=%5.2f' % (lig, tgt, ref, q, d), flush=True)
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
