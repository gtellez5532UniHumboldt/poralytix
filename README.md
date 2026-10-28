# 🧬 Poralytix — Estación de Trabajo de Cinética de Citometría de Flujo (.fcs)

> **Poralytix** (*pore* + *analytics*) es la versión refactorizada y modular del analizador
> de cinética de permeabilización de membrana. La versión original (script único) se
> conserva en el repositorio [`analizador-citometria`](https://github.com/gtellez5532UniHumboldt/analizador-citometria).
>
> Los datos experimentales (`.fcs`) y las carpetas de resultados (`Resultados_*`) no se
> versionan (ver `.gitignore`).

**Versión:** 2.0.0 (Arquitectura Modular Refactorizada)  
**Entorno Validado:** Python 3.11 en Windows x64  
**Dominio:** Análisis biofísico de permeabilización de membrana celular por péptidos citolíticos (Calceína-AM)

---

## 1. Estructura del Proyecto

```
poralytix/
├── main.py                     # Punto de entrada de la Interfaz Gráfica (GUI)
├── cli.py                      # Herramienta de procesamiento batch por línea de comandos (Headless)
├── requirements.txt            # Dependencias científicas
├── README.md                   # Este manual
├── cytoworkstation/
│   ├── __init__.py
│   ├── config.py               # Paletas cromáticas, PEPTIDE_ALIASES (sat1, 813, 820, CZS5), constantes
│   ├── io/
│   │   ├── fcs_parser.py       # Parser binario FCS 2.0 / 3.0 / 3.1 (sin librerías externas)
│   │   ├── metadata.py         # Extractor robusto de metadatos (soporte para subcarpetas y CZS5)
│   │   └── session.py          # Serialización y deserialización de sesiones .json
│   ├── core/
│   │   ├── statistics.py       # Intervalos de Wilson 95%, Cohen's h, Cochran-Armitage, Fisher/Chi2
│   │   ├── models.py           # Modelos no lineales: Hill, Michaelis-Menten, Gompertz, Weibull, Bi-exp
│   │   ├── kinetics.py         # Gating morfológico, doble cut-off V1L/V2L, T_lag, T_max, V_active
│   │   └── synthetic_demo.py   # Generador determinista de archivos FCS sintéticos
│   ├── export/
│   │   ├── tables_export.py    # Exportación tabular a CSV / TSV
│   │   ├── plots_600dpi.py     # Renderizado de figuras a 600 DPI (TIFF / PNG)
│   │   └── pdf_report.py       # Motor de reportes PDF vectoriales auditables (PDF 1.4)
│   └── gui/
│       ├── widgets/            # Canvases interactivos Tkinter (histograma, cinéticas, saturación)
│       ├── dialogs/            # Modales de metadatos, fórmulas matemáticas y guía F1
│       └── app.py              # Controlador principal de la Workstation (5 pestañas)
└── tests/
    ├── test_metadata.py        # Pruebas de deducción de metadatos sobre datasets de prueba
    ├── test_fcs_parser.py      # Pruebas de lectura binaria de archivos reales
    ├── test_statistics.py      # Validación de métodos bioestadísticos
    └── test_models.py          # Validación de convergencia de modelos no lineales
```

---

## 2. Instalación y Requisitos

Requiere Python 3.11 con las siguientes librerías:
```bash
pip install -r requirements.txt
```
Librerías principales:
- `numpy >= 1.24.0`
- `scipy >= 1.10.0`
- `pandas >= 2.0.0`
- `matplotlib >= 3.7.0`
- `reportlab >= 4.0.0`

---

## 3. Modo de Uso

### Interfaz Gráfica (GUI)
Para iniciar la estación de trabajo con todas sus herramientas visuales e interactivas:
```powershell
py -3.11 main.py
```
**Flujo de Trabajo en las 5 Pestañas:**
1. **Pestaña 1 (Carga):** Sube archivos `.fcs`, edita metadatos en lote o genera el set demo sintético.
2. **Pestaña 2 (Pre-Análisis):** Ajusta los Cut-Offs V1L (exclusión) y V2L (lisis) arrastrando directamente con el ratón sobre el histograma. Ajusta las ventanas $T_{lag}$ y $T_{max}$ con deslizadores o sobre las curvas.
3. **Pestaña 3 (Panel Categórico):** Evalúa el balance temporal de eventos Afectados vs No Afectados, mapa de calor espaciotemporal, pruebas de Cochran-Armitage e isocronas $EC_{50}(t)$.
4. **Pestaña 4 (Modelos):** Ajusta las curvas cinéticas con el modelo de Hill regularizado ($n_H \le 6.0$), Michaelis-Menten, Gompertz o Weibull. Guarda la sesión en `.json` (Ctrl+S).
5. **Pestaña 5 (Exportación):** Genera el reporte formal en PDF y figuras de alta resolución a 600 DPI.

### Modo Batch por Consola (CLI Headless)
Para procesar carpetas completas sin interfaz gráfica:
```powershell
py -3.11 cli.py --input-dir ../CZS5c- --output-dir Resultados_CZS5 --target-cell THP-1 --control-cell PBMCs --model "Allosteric Hill"
```

---

## 4. Ejecución de la Suite de Pruebas Automatizadas

Para verificar la integridad del parser, la heurística de metadatos y los modelos matemáticos:
```powershell
# Desde la raíz del repositorio
py -3.11 -m unittest discover -s tests -t . -v
```

---

## 5. Modelos Matemáticos y Ecuaciones Biofísicas

### Modelo Alostérico de Hill Regularizado (Por defecto)
$$V_{active}(C) = \frac{V_{max} \cdot C^{n_H}}{K_d^{n_H} + C^{n_H}}$$
- $V_{max}$ (%/s): Velocidad máxima de permeabilización a saturación.
- $K_d$ (µM): Concentración molar de semi-saturación ($EC_{50}$ cinética).
- $n_H$: Coeficiente de Hill en $[0.5, 6.0]$ que refleja la cooperatividad del ensamblado oligomérico del poro.

### Intervalo de Confianza Binomial Exacto de Wilson (95%)
$$Center = \frac{p + \frac{z^2}{2N}}{1 + \frac{z^2}{N}}, \quad Half = \frac{z}{1 + \frac{z^2}{N}} \sqrt{\frac{p(1-p)}{N} + \frac{z^2}{4N^2}}$$
Garantiza propiedades asintóticas y de cobertura exactas aun cuando las proporciones se aproximan a 0% o 100%.
