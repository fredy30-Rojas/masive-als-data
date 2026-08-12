
import base64, glob, os
from google.colab import files
WORK = '/content/masive_als'
archivos = sorted(glob.glob(WORK + '/ligandos_reparados/*.pdbqt'))
print('@@TOTAL@@', len(archivos), flush=True)
for a in archivos:
    txt = open(a, errors='replace').read()
    b = base64.b64encode(txt.encode('utf-8')).decode('ascii')
    print('@@ARCHIVO@@ ' + os.path.basename(a), flush=True)
    for i in range(0, len(b), 2500):
        print('@@P@@' + b[i:i+2500], flush=True)
    print('@@FIN_A@@', flush=True)
print('@@CSV@@', flush=True)
print(open(WORK + '/resultados/resultados_colab.csv').read(), flush=True)
print('@@FIN@@', flush=True)
try:
    files.download(WORK + '/resultados/resultados_colab.csv')
except Exception as e:
    print('DL_ERR', str(e)[:60], flush=True)
