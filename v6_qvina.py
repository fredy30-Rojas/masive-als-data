# CELDA 5: Docking con QVina 2.1 (GPU/CPU) + validacion de consistencia vs Vina 1.2.5
# Adaptado de v6_docking.py: cambia SOLO el binario de acople, mantiene el resto.
# - Mismo formato CSV: ligand,target,energy,timestamp
# - Subconjuntos pares/impares, reparacion con Open Babel, reanudacion por CSV
# - VALIDAR=True: acopla una muestra de pares YA hechos (seed 42) y compara con el CSV
import csv, glob, json, os, shutil, subprocess, sys, tarfile, time, urllib.request

WORK = '/content/masive_als'
# QVina 2.1 (binario x86_64). Descarga desde GitHub (repo masive-als-data) o dataset Kaggle.
QVINA_BIN = '/content/qvina_linux'
QVINA_URL = 'https://raw.githubusercontent.com/fredy30-Rojas/masive-als-data/main/qvina_linux'

# --- 0) Verificar GPU (requisito si se usa Vina-GPU; informativo para QVina2 CPU) ---
print('=== VERIFICACION GPU ===', flush=True)
gpu = subprocess.run(['nvidia-smi'], capture_output=True, text=True, timeout=30)
print('nvidia-smi rc=%d' % gpu.returncode, flush=True)
if gpu.returncode == 0:
    print('GPU DETECTADA:', [l for l in gpu.stdout.splitlines() if 'Tesla' in l or 'NVIDIA' in l][:2], flush=True)
    GPU_OK = True
else:
    print('SIN GPU (nvidia-smi fallo). QVina2 CPU funciona; Vina-GPU NO.', flush=True)
    GPU_OK = False

# --- 1) Descargar binario QVina ---
if not os.path.exists(QVINA_BIN):
    print('Descargando binario QVina 2.1...', flush=True)
    urllib.request.urlretrieve(QVINA_URL, QVINA_BIN)
    os.chmod(QVINA_BIN, 0o755)
print('Binario QVina listo', flush=True)

# El binario conda-forge necesita libboost; instalarlas si faltan
try:
    subprocess.run(['ldconfig', '-p'], capture_output=True, timeout=15)
    boost_ok = subprocess.run(
        [QVINA_BIN, '--help'], capture_output=True, text=True, timeout=30).returncode == 0
    if not boost_ok:
        print('Faltan libboost - instalando via apt...', flush=True)
        subprocess.run(['apt-get', 'update', '-qq'], timeout=300, capture_output=True)
        subprocess.run(['apt-get', 'install', '-y', '-qq',
                        'libboost-program-options-dev', 'libboost-filesystem-dev',
                        'libboost-thread-dev', 'libboost-system-dev'],
                       timeout=600, capture_output=True)
except Exception as ex:
    print('Nota boost:', str(ex)[:100], flush=True)
print(subprocess.run([QVINA_BIN, '--help'], capture_output=True, text=True, timeout=30).stdout[:600], flush=True)

# --- 2) Datos frescos (igual que v6) ---
for sub in ['receptores', 'ligandos', 'resultados', 'ligandos_reparados']:
    os.makedirs(WORK + '/' + sub, exist_ok=True)
URLS = [
    'https://raw.githubusercontent.com/fredy30-Rojas/masive-als-data/main/colab_receptores_plano.tar.gz',
    'https://raw.githubusercontent.com/fredy30-Rojas/masive-als-data/main/colab_ligandos50_plano.tar.gz',
]
for url in URLS:
    dest = '/content/' + os.path.basename(url)
    if os.path.exists(dest):
        os.remove(dest)
    urllib.request.urlretrieve(url, dest)
extract_dir = '/content/tmp_extract'
if os.path.exists(extract_dir):
    shutil.rmtree(extract_dir)
os.makedirs(extract_dir)
for tar in sorted(glob.glob('/content/*.tar.gz')):
    with tarfile.open(tar) as t:
        try:
            t.extractall(extract_dir, filter='data')
        except TypeError:
            t.extractall(extract_dir)
def es_receptor(n):
    return n in ('TDP43.pdbqt', 'SOD1.pdbqt', 'FUS.pdbqt')
for d in ('receptores', 'ligandos'):
    for f in glob.glob(WORK + '/' + d + '/*.pdbqt'):
        os.remove(f)
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

# --- 4) FUS plano ---
fus = RECEPTORES['FUS']['archivo']
try:
    lines = open(fus).read().splitlines()
    ini = next((i for i, l in enumerate(lines) if l.startswith('MODEL')), None)
    fin = next((i for i, l in enumerate(lines) if l.strip() == 'ENDMDL'), None)
    keep = lines[ini + 1:fin] if (ini is not None and fin is not None and fin >= ini) else lines
    if not any(l.startswith(('ATOM', 'HETATM')) for l in keep):
        keep = lines
    open(fus, 'w').write('\n'.join(keep) + '\n')
except Exception as ex:
    print('ERROR reparando FUS:', str(ex)[:100], flush=True)

# --- 5) CSV (reanuda) ---
CSV = WORK + '/resultados/resultados_qvina.csv'
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

# --- 6) Validacion: muestra de 8-10 pares YA hechos (compara con resultados_colab.csv) ---
VALIDAR = True   # True = primero valida consistencia; False = tanda completa
N_MUESTRA = 10
if VALIDAR:
    base_csv = WORK + '/resultados/resultados_colab.csv'
    if not os.path.exists(base_csv):
        print('Descargando resultados_colab.csv para comparar...', flush=True)
        urllib.request.urlretrieve(
            'https://raw.githubusercontent.com/fredy30-Rojas/masive-als-data/main/resultados_colab.csv',
            base_csv)
    refs = {}
    for r in csv.DictReader(open(base_csv)):
        refs.setdefault(r['ligand'], {})[r['target']] = float(r['energy'])
    # tomar hasta N_MUESTRA ligandos con las 3 proteinas hechas
    muestra = sorted([k for k, v in refs.items() if len(v) >= 3])[:N_MUESTRA]
    pendientes = []
    for nombre in muestra:
        for target in RECEPTORES:
            if target not in hechos.get(nombre, set()):
                pendientes.append((WORK + '/ligandos/' + nombre + '.pdbqt', nombre, target))
    print('VALIDACION: %d ligandos x 3 = %d pares a re-acoplar con QVina' % (len(muestra), len(pendientes)), flush=True)

# --- 7) Reparacion con Open Babel (igual que v6) ---
def instalar_obabel():
    try:
        from openbabel import openbabel as ob
        return ob
    except ImportError:
        subprocess.run([sys.executable, '-m', 'pip', 'install', '-q', 'openbabel-wheel'], timeout=600)
        from openbabel import openbabel as ob
        return ob
def smiles_chembl(nombre):
    url = 'https://www.ebi.ac.uk/chembl/api/data/molecule/%s.json' % nombre
    req = urllib.request.Request(url, headers={'User-Agent': 'masive-als/1.0'})
    with urllib.request.urlopen(req, timeout=30) as r:
        d = json.load(r)
    return ((d.get('molecule_structures') or {}).get('canonical_smiles') or '').strip()
def reparar_pdbqt(ob, nombre, lig):
    try:
        s = smiles_chembl(nombre)
        if not s:
            return False
        frags = sorted([f.strip() for f in s.split('.') if len(f.strip()) > 1], key=len, reverse=True)
        if frags:
            s = frags[0]
        conv = ob.OBConversion()
        conv.SetInAndOutFormats('smi', 'pdbqt')
        mol = ob.OBMol()
        if not conv.ReadString(mol, s):
            return False
        builder = ob.OBBuilder()
        builder.Build(mol)
        ob.OBChargeModel.FindType('gasteiger').ComputeCharges(mol)
        txt = conv.WriteString(mol)
        if not any(l.startswith(('ATOM', 'HETATM')) for l in txt.splitlines()):
            return False
        open(lig, 'w').write(txt)
        shutil.copy(lig, WORK + '/ligandos_reparados/' + nombre + '.pdbqt')
        return True
    except Exception:
        return False

# --- 8) Acoplar con QVina (mismos flags que Vina; control de hilos via --threads) ---
def acoplar(lig, target, exhaustiveness=2, threads=2):
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
    variantes = [base + ['--threads', str(threads)], base]  # con threads primero, fallback sin el
    try:
        for cmd in variantes:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
            out = (r.stdout or '') + (r.stderr or '')
            if r.returncode == 0:
                break
        if r.returncode != 0:
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

ob = None
reparados = set()
if pendientes:
    ob = instalar_obabel()
    t0 = time.time()
    n_ok = 0
    for lig, nombre, target in pendientes:
        if not tiene_atomos(lig):
            if nombre not in reparados:
                reparados.add(nombre)
                if reparar_pdbqt(ob, nombre, lig):
                    print('REPARADO', nombre, flush=True)
                else:
                    with open(SALTADOS, 'a') as f:
                        f.write(nombre + '\n')
                    continue
        energia, err = acoplar(lig, target)
        if energia is None and nombre not in reparados:
            reparados.add(nombre)
            if reparar_pdbqt(ob, nombre, lig):
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

# --- 9) Comparacion con Vina 1.2.5 (si validacion) ---
if VALIDAR:
    print()
    print('=== COMPARACION QVina vs Vina 1.2.5 (seed 42) ===', flush=True)
    diffs = []
    for r in csv.DictReader(open(CSV)):
        q = float(r['energy'])
        ref = refs.get(r['ligand'], {}).get(r['target'])
        if ref is not None:
            d = abs(q - ref)
            diffs.append((r['ligand'], r['target'], ref, q, d))
    if diffs:
        for lig, tgt, ref, q, d in sorted(diffs, key=lambda x: -x[4])[:15]:
            print('%-22s %-5s Vina=%7.2f QVina=%7.2f diff=%5.2f' % (lig, tgt, ref, q, d), flush=True)
        mean = sum(d[4] for d in diffs) / len(diffs)
        maxd = max(d[4] for d in diffs)
        print('MEDIA diff=%.3f | MAX diff=%.3f | n=%d' % (mean, maxd, len(diffs)), flush=True)
        # RANKING: top de cada motor
        vina_top = sorted([(r['ligand'], r['target'], float(r['energy'])) for r in
                           csv.DictReader(open(WORK + '/resultados/resultados_colab.csv'))
                           if r['ligand'] in {d[0] for d in diffs}], key=lambda x: x[2])[:5]
        qvina_top = sorted([(r['ligand'], r['target'], float(r['energy'])) for r in
                            csv.DictReader(open(CSV))], key=lambda x: x[2])[:5]
        print('Top Vina :', ['%s/%s %.2f' % t for t in vina_top], flush=True)
        print('Top QVina:', ['%s/%s %.2f' % t for t in qvina_top], flush=True)
        if maxd > 1.5:
            print('ALERTA: diferencias grandes - NO mezclar motores sin revisar.', flush=True)
        else:
            print('OK: diferencias pequenas, motores compatibles.', flush=True)

print(flush=True)
print('=== TANDAS COMPLETADAS ===', flush=True)
print('Filas en CSV:', len(list(csv.DictReader(open(CSV)))), flush=True)
