"""
Modal con la Guía de Usuario Paso a Paso del Software.
Soporte bilingüe Español e Inglés.
"""
import tkinter as tk
from tkinter import ttk

from ...i18n import get_current_language


def open_user_guide_dialog(parent, lang=None):
    """Abre la ventana interactiva con la Guía Paso a Paso de Uso del Software."""
    active_lang = lang or get_current_language()
    dlg = tk.Toplevel(parent)
    
    if active_lang == 'en':
        dlg.title("📖 Step-by-Step Flow Cytometry Kinetics Software User Guide")
        hdr_txt = "📖 SOFTWARE PRACTICAL USER GUIDE - WORKFLOW (TABS 1 TO 5)"
        steps = [
            (" 📁 Step 1: Ingestion & Metadata ",
             "STEP 1: DATA INGESTION (.FCS) AND METADATA (TAB 1)\n\n"
             "1. File Ingestion:\n"
             "   • Click 'Upload .fcs Files' to select your experimental data.\n"
             "   • Or click 'Generate & Load Demo Set' to explore the workstation with simulated data.\n\n"
             "2. Recommended Nomenclature:\n"
             "   • Format: [Compound]_[CellLine]_[Concentration]_[CC_optional]_[Replicate].fcs\n"
             "   • Baseline Control Example (-60 to 0 s): THP1_CC_sat1_0.7_1.fcs\n"
             "   • Experimental Sample Example (0 to 600 s): THP1_sat1_0.7_1.fcs\n\n"
             "3. Metadata Editing:\n"
             "   • Double-click any table cell to edit Compound, Dose, Cell or Role.\n"
             "   • Select multiple rows and use 'Batch Edit Selected' for bulk changes.\n\n"
             "4. Execution:\n"
             "   • Select permeability channel (FITC/FL1) and click 'RUN PRE-ANALYSIS'."),
            (" 🔍 Step 2: Dynamic Pre-Analysis ",
             "STEP 2: DYNAMIC PRE-ANALYSIS AND CUT-OFF GATING (TAB 2)\n\n"
             "1. Fluorescence Histogram (Box 1):\n"
             "   • Dotted Purple Line (Cut-Off 1 V1L): Events < threshold are debris or unloaded cells.\n"
             "   • Dashed Red Line (Cut-Off 2 V2L): Classifies Affected (< threshold) vs Non-Affected (>= threshold).\n"
             "   • Drag both lines directly on the canvas to adjust thresholds in real time.\n\n"
             "2. Time Windows T_lag and T_max (Boxes 2 & 3):\n"
             "   • Drag the blue line (T_lag) and red line (T_max) on the Permeabilization Plot 1.\n"
             "   • The active rate V_active (%/s) recalculates immediately.\n\n"
             "3. Sample Cut-Offs Table (Box 4):\n"
             "   • Double-click any numeric cell to assign sample-specific gating thresholds."),
            (" 🧬 Step 3: Categorical Panel ",
             "STEP 3: CATEGORICAL PANEL AFFECTED vs NON-AFFECTED (TAB 3)\n\n"
             "1. Visualization Modes:\n"
             "   • Complementary Time Balance: Evolution of % V2L and % V2R with milestone markers.\n"
             "   • Milestones Comparison: Bar charts for T_onset, effect T50 and T_effect_max by dose.\n"
             "   • Heatmap: Matrix of % V2L by Compound × Dose × Time Window.\n"
             "   • Instantaneous Rate: Derivative d(%V2L)/dt of cellular conversion.\n\n"
             "2. Table Sub-tabs:\n"
             "   • Effect Milestones, Proportions by Window, Trend Tests and Biophysical Interpretation."),
            (" 📊 Step 4: Models & Statistics ",
             "STEP 4: KINETIC MODELS AND PRE-REPORT (TAB 4)\n\n"
             "1. Mathematical Model Selection:\n"
             "   • Choose between Allosteric Hill, Michaelis-Menten, Gompertz, Weibull, or Bi-Exponential.\n"
             "   • Inspect evaluated equations, R², and SEM bars in real time.\n\n"
             "2. Save & Load Sessions:\n"
             "   • Use 'Save Session' (Ctrl+S) to export complete state into an auditable .json file."),
            (" 📄 Step 5: Final Report ",
             "STEP 5: FINAL AUDIT REPORT AND EXPORT (TAB 5)\n\n"
             "1. Report Generation:\n"
             "   • Set report title, author and target output directory.\n"
             "   • Click 'EXPORT ALL (600 DPI PLOTS + CSV TABLES + PDF REPORT)'.\n"
             "   • Automatically produces 600 DPI figures, CSV tables and an executive PDF document.")
        ]
    else:
        dlg.title("📖 Guía Paso a Paso de Uso del Software de Cinética de Citometría")
        hdr_txt = "📖 GUÍA PRÁCTICA DE USO DEL SOFTWARE - FLUJO DE TRABAJO (PESTAÑAS 1 A 5)"
        steps = [
            (" 📁 Paso 1: Carga y Metadatos ",
             "PASO 1: CARGA DE DATOS (.FCS) Y METADATOS (PESTAÑA 1)\n\n"
             "1. Carga de Archivos:\n"
             "   • Haz clic en 'Subir Archivos .fcs' para seleccionar tus datos experimentales.\n"
             "   • O pulsa 'Generar y Cargar Set Demo' para explorar el sistema con datos simulados.\n\n"
             "2. Nomenclatura Recomendada:\n"
             "   • Formato: [Compuesto]_[LineaCelular]_[Concentracion]_[CC_opcional]_[Replica].fcs\n"
             "   • Ejemplo Control Basal (-60 a 0 s): THP1_CC_sat1_0.7_1.fcs\n"
             "   • Ejemplo Muestra Experimental (0 a 600 s): THP1_sat1_0.7_1.fcs\n\n"
             "3. Edición de Metadatos:\n"
             "   • Doble clic sobre cualquier celda de la tabla para editar Compuesto, Dosis, Célula o Rol.\n"
             "   • Selecciona varias filas y usa 'Editar Metadatos de Seleccionados' para edición en lote.\n\n"
             "4. Ejecución:\n"
             "   • Selecciona el canal de permeabilidad (FITC/FL1) y haz clic en 'EJECUTAR PRE-ANÁLISIS'."),
            (" 🔍 Paso 2: Pre-Análisis Dinámico ",
             "PASO 2: PRE-ANÁLISIS DINÁMICO Y CORTE DE CUT-OFFS (PESTAÑA 2)\n\n"
             "1. Histograma de Fluorescencia (Caja 1):\n"
             "   • Línea Morada Punteada (Cut-Off 1 V1L): Eventos < umbral son detritos o células no cargadas.\n"
             "   • Línea Roja Discontinua (Cut-Off 2 V2L): Clasifica en Afectados (< umbral) vs No Afectados (>= umbral).\n"
             "   • Arrastra ambas líneas directamente con el ratón sobre el lienzo para ajustar los cortes en tiempo real.\n\n"
             "2. Ventanas Temporales T_lag y T_max (Caja 2 y 3):\n"
             "   • Arrastra la línea azul (T_lag) y la línea roja (T_max) en la Gráfica 1 de Permeabilización.\n"
             "   • La tasa activa V_active (%/s) se recalcula instantáneamente.\n\n"
             "3. Tabla de Cut-Offs por Muestra (Caja 4):\n"
             "   • Doble clic sobre cualquier celda numérica para asignar umbrales individuales a una muestra."),
            (" 🧬 Paso 3: Panel Categórico ",
             "PASO 3: PANEL CATEGÓRICO AFECTADOS vs NO AFECTADOS (PESTAÑA 3)\n\n"
             "1. Modos de Visualización:\n"
             "   • Balance Temporal Complementario: Curvas de evolución de % V2L y % V2R con marcadores de hitos.\n"
             "   • Comparativa de Hitos: Barras de T_onset, T50 del efecto y T_efecto_máx por dosis.\n"
             "   • Mapa de Calor: Matriz de % V2L por Compuesto × Dosis × Ventana de Tiempo.\n"
             "   • Velocidad Instantánea: Derivada d(%V2L)/dt de conversión celular.\n\n"
             "2. Sub-pestañas de Tablas:\n"
             "   • Hitos del Efecto, Proporciones por Ventana, Tests de Tendencia e Interpretación Biofísica."),
            (" 📊 Paso 4: Modelos y Estadística ",
             "PASO 4: MODELOS CINÉTICOS Y PRE-REPORTE (PESTAÑA 4)\n\n"
             "1. Selección de Modelos Matemáticos:\n"
             "   • Elige entre Allosteric Hill, Michaelis-Menten, Gompertz, Weibull o Bi-Exponencial.\n"
             "   • Observa en tiempo real la ecuación evaluada, el R² y las barras SEM.\n\n"
             "2. Guardar y Cargar Sesión:\n"
             "   • Usa 'Guardar Sesión' (Ctrl+S) para exportar todo el estado a un archivo .json auditable."),
            (" 📄 Paso 5: Reporte Final ",
             "PASO 5: REPORTE FINAL Y EXPORTACIÓN (PESTAÑA 5)\n\n"
             "1. Generación del Reporte:\n"
             "   • Define el título del informe, autor y carpeta de salida.\n"
             "   • Pulsa 'GENERAR Y EXPORTAR REPORTE FINAL EN PDF E IMÁGENES 600 DPI'.\n"
             "   • Se generan automáticamente todas las figuras TIFF a 600 DPI, tablas CSV/TSV y el PDF.")
        ]

    dlg.geometry("960x660")
    dlg.transient(parent)

    hdr = tk.Frame(dlg, bg="#002B49", pady=6, padx=12)
    hdr.pack(fill="x")
    tk.Label(
        hdr,
        text=hdr_txt,
        bg="#002B49", fg="white", font=("DejaVu Sans", 10, "bold")
    ).pack(side="left")

    nb = ttk.Notebook(dlg)
    nb.pack(fill="both", expand=True, padx=10, pady=8)

    for title, body in steps:
        frm = ttk.Frame(nb, padding=12)
        nb.add(frm, text=title)
        txt = tk.Text(frm, wrap="word", font=("DejaVu Sans", 9), padx=10, pady=10, bg="#FDFEFE", relief="flat")
        txt.insert("1.0", body)
        txt.config(state="disabled")
        txt.pack(fill="both", expand=True)

    return dlg
