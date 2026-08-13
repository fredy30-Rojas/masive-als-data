# Prompt para Buffy: SATURAR las GPUs diarias (Colab T4 + Kaggle)

Fredy quiere aprovechar al MÁXIMO las horas de GPU gratis cada día. Vina-GPU ya está validado: los motores son compatibles con Vina 1.2.5 (diferencia media 0,60, máxima 1,31, umbral 1,5 → OK). Ahora toca trabajar al 100%. Léelo todo antes de empezar.

## Estado actual (13 agosto 2026, tarde)

- **135 ligandos listos en el repo** `fredy30-Rojas/masive-als-data`:
  - 50 ChEMBL (colab_ligandos50_plano.tar.gz)
  - 85 fármacos conocidos de ELA (ligandos_ela_plano.tar.gz, NUEVO)
- 3 receptores: TDP43 (4IUF), SOD1, FUS — cajas centradas de la literatura.
- Notebooks actualizados para **tanda completa** (`VALIDAR = False`), la validación ya pasó:
  - `MASIVE-ALS_Colab_GPU_VinaGPU.ipynb` → subconjunto `impares` (68 ligandos × 3 = ~204 pares)
  - `MASIVE-ALS_Kaggle_GPU_VinaGPU.ipynb` → subconjunto `pares` (67 × 3 = ~201 pares) + descarga de los 85 ELA desde GitHub (bloque 2c)
- Velocidad real medida: **12 pares en 1,2 min en T4** (~10 pares/min) → tanda completa ≈ 20-25 min por GPU.

## Misión: GPUs trabajando todo el tiempo que dure la sesión

1. **Colab T4 (cuenta Yograbotodo@gmail.com):**
   - Abrir `MASIVE-ALS_Colab_GPU_VinaGPU.ipynb` desde Drive, entorno T4 GPU.
   - Ejecutar celdas en orden. Verificar: `nvidia-smi` muestra la T4, "Ligandos (subconjunto=impares): 68", "Pares pendientes: ~204".
   - Al terminar: copiar `resultados_vinagpu_colab.csv` al PC local (`Desktop\MASIVE-ALS-Colab\resultados\`).

2. **Kaggle (cuenta yograbotodo, kernel masive-als-vinagpu):**
   - RECORDATORIO: la última vez `nvidia-smi` falló (posible verificación de teléfono). Si la GPU sigue sin aparecer: avisar a Fredy y NO perder tiempo — toda la tanda va a Colab (impares + pares, ajustar SUBCONJUNTO='todos').
   - Si hay GPU: ejecutar, verificar "Ligandos totales con ELA: 135" y pares pendientes ~201.
   - Al terminar: copiar `resultados_vinagpu_kaggle.csv` al PC local.

3. **Mientras las GPUs acoplan, preparar el siguiente lote** (para que no paren):
   - Generar más PDBQT: compuestos ChEMBL con bioactividad en TDP-43/SOD1/FUS (o la lista que tenga el pipeline) → Open Babel a PDBQT plano (una sola conformación, sin MODEL/ENDMDL).
   - Empaquetar como `ligandos_batch2_plano.tar.gz` y subirlo al repo.
   - NOTA IMPORTANTE: los PDBQT deben tener átomos (evitar los que den archivos vacíos) y nombres únicos.

4. **Cuando termine la tanda y exista batch2:** volver a ejecutar los notebooks. El CSV de checkpoint reanuda solo: los pares ya hechos se saltan. Cada re-ejecución consume la sesión de GPU disponible.

5. **Fusionar resultados** (cuando haya resultados de ambas plataformas):
   - Mezclar `resultados_vinagpu_colab.csv` + `resultados_vinagpu_kaggle.csv` en `resultados_vinagpu_total.csv` (sin duplicados: mismo ligand+target).
   - Subir a GitHub masive-als-data y avisar a Fredy con el Top 10.

## Reglas

- NUNCA incrustar tokens ni claves en notebooks, kernels o repos públicos. Los datos y binarios viajan por GitHub público o dataset privado.
- Preguntar a Fredy ANTES de re-pushear kernels de Kaggle o crear versiones nuevas del dataset.
- Si una sesión de Colab muere a mitad: relanzar y re-ejecutar — el checkpoint salta lo hecho. No volver a acoplar pares ya hechos.
- No mezclar resultados de motores sin revisar (Vina-GPU ya validado; si se añade QVina, validar igual que antes).
- Reportar a Fredy: pares acoplados por plataforma, tiempo, Top 10 de afinidades, y qué falta para el siguiente lote.
