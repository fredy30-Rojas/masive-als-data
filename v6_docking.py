# CELDA 5: Docking impares (binario Vina) + auto-reparacion de PDBQT corruptos
import csv, glob, json, os, shutil, subprocess, sys, tarfile, time, urllib.request
WORK = '/content/masive_als'
VINA_BIN = '/content/vina_bin'

# --- 1) Descargar el binario oficial de Vina ---
if not os.path.exists(VINA_BIN):
    url = 'https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64'
    print('Descargando binario Vina...', flush=True)
    urllib.request.urlretrieve(url, VINA_BIN)
    os.chmod(VINA_BIN, 0o755)
print('Binario Vina listo', flush=True)

# --- 2) Auto-descarga de datos SIEMPRE frescos ---
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
    print('Descargando:', os.path.basename(url), flush=True)
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

# --- 4) FUS.pdbqt plano (SIN MODEL/ENDMDL) ---
fus = RECEPTORES['FUS']['archivo']
try:
    lines = open(fus).read().splitlines()
    ini = next((i for i, l in enumerate(lines) if l.startswith('MODEL')), None)
    fin = next((i for i, l in enumerate(lines) if l.strip() == 'ENDMDL'), None)
    if ini is not None and fin is not None and fin >= ini:
        keep = lines[ini + 1:fin]
    else:
        keep = lines
    if not any(l.startswith(('ATOM', 'HETATM')) for l in keep):
        keep = lines
    open(fus, 'w').write('\n'.join(keep) + '\n')
    print('FUS plano: ATOM=%d' % sum(1 for l in keep if l.startswith(('ATOM', 'HETATM'))), flush=True)
except Exception as ex:
    print('ERROR reparando FUS:', str(ex)[:100], flush=True)

# --- 5) CSV (reanuda) ---
CSV = WORK + '/resultados/resultados_colab.csv'
if not os.path.exists(CSV):
    with open(CSV, 'w', newline='') as f:
        csv.writer(f).writerow(['ligand', 'target', 'energy', 'timestamp'])
SALTADOS = WORK + '/resultados/saltados.txt'
saltados = set()
if os.path.exists(SALTADOS):
    saltados = set(l.strip() for l in open(SALTADOS) if l.strip())

hechos = {}
for r in csv.DictReader(open(CSV)):
    hechos.setdefault(r['ligand'], set()).add(r['target'])

ligs = sorted(glob.glob(WORK + '/ligandos/*.pdbqt'))
# --- subconjunto: 'pares', 'impares' o 'todos' ---
SUBCONJUNTO = 'impares'
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

# --- 6) Reparacion con Open Babel (SMILES de ChEMBL) ---
def instalar_obabel():
    try:
        from openbabel import openbabel as ob
        return ob
    except ImportError:
        print('Instalando openbabel-wheel...', flush=True)
        subprocess.run([sys.executable, '-m', 'pip', 'install', '-q', 'openbabel-wheel'], timeout=600)
        from openbabel import openbabel as ob
        return ob

def smiles_chembl(nombre):
    url = 'https://www.ebi.ac.uk/chembl/api/data/molecule/%s.json' % nombre
    req = urllib.request.Request(url, headers={'User-Agent': 'masive-als/1.0'})
    with urllib.request.urlopen(req, timeout=30) as r:
        d = json.load(r)
    s = (d.get('molecule_structures') or {}).get('canonical_smiles') or ''
    return s.strip()

def reparar_pdbqt(ob, nombre, lig):
    """Regenera el PDBQT desde SMILES. Devuelve True si quedo con atomos."""
    try:
        s = smiles_chembl(nombre)
        if not s:
            print('  SIN_SMILES', nombre, flush=True)
            return False
        # FIX sales/fragmentos: conservar solo el fragmento mas grande
        frags = sorted([f.strip() for f in s.split('.') if len(f.strip()) > 1],
                       key=len, reverse=True)
        if frags:
            s = frags[0]
        conv = ob.OBConversion()
        conv.SetInAndOutFormats('smi', 'pdbqt')
        mol = ob.OBMol()
        if not conv.ReadString(mol, s):
            print('  SMILES_ILEGIBLE', nombre, flush=True)
            return False
        builder = ob.OBBuilder()
        builder.Build(mol)
        ob.OBChargeModel.FindType('gasteiger').ComputeCharges(mol)
        txt = conv.WriteString(mol)
        if not any(l.startswith(('ATOM', 'HETATM')) for l in txt.splitlines()):
            print('  SIN_ATOMOS', nombre, flush=True)
            return False
        open(lig, 'w').write(txt)
        # respaldo para descargar
        shutil.copy(lig, WORK + '/ligandos_reparados/' + nombre + '.pdbqt')
        return True
    except Exception as ex:
        print('  REPARAR_ERR', nombre, str(ex)[:80], flush=True)
        return False

# --- 7) Acoplar cada par en un SUBPROCESO ---
def acoplar(lig, target):
    info = RECEPTORES[target]
    cmd = [VINA_BIN,
           '--receptor', info['archivo'],
           '--ligand', lig,
           '--center_x', str(info['centro'][0]),
           '--center_y', str(info['centro'][1]),
           '--center_z', str(info['centro'][2]),
           '--size_x', str(info['tamano'][0]),
           '--size_y', str(info['tamano'][1]),
           '--size_z', str(info['tamano'][2]),
           '--exhaustiveness', '2', '--num_modes', '3', '--seed', '42']
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        out = (r.stdout or '') + (r.stderr or '')
        if r.returncode != 0:
            return None, 'vina rc=%d %s' % (r.returncode, out[-150:])
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
                print('REPARANDO_VACIO', nombre, flush=True)
                if reparar_pdbqt(ob, nombre, lig):
                    print('REPARADO', nombre, flush=True)
                else:
                    print('REPARAR_FALLO', nombre, '- se salta', flush=True)
                    with open(SALTADOS, 'a') as f:
                        f.write(nombre + '\n')
                    continue
        energia, err = acoplar(lig, target)
        if energia is None and nombre not in reparados:
            reparados.add(nombre)
            print('REPARANDO_TRAS_ERROR', nombre, target, flush=True)
            if reparar_pdbqt(ob, nombre, lig):
                print('REPARADO', nombre, flush=True)
                energia, err = acoplar(lig, target)
        if energia is not None:
            with open(CSV, 'a', newline='') as f:
                csv.writer(f).writerow([nombre, target, energia, time.strftime('%Y-%m-%d %H:%M:%S')])
            n_ok += 1
        else:
            print('ERROR', nombre, target, err, flush=True)
        if n_ok % 15 == 0 and n_ok > 0:
            print('[%d acoplados] %.1f min' % (n_ok, (time.time() - t0) / 60), flush=True)
    print('Acoplados ahora:', n_ok, flush=True)

print(flush=True)
print('=== TANDAS COMPLETADAS ===', flush=True)
print('Filas en CSV:', len(list(csv.DictReader(open(CSV)))), flush=True)
print('PDBQT reparados:', sorted(os.path.basename(p) for p in glob.glob(WORK + '/ligandos_reparados/*.pdbqt')), flush=True)