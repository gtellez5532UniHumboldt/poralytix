"""
Modal con la Guía de Usuario Paso a Paso del Software.
"""
import tkinter as tk
from tkinter import ttk


def open_user_guide_dialog(parent):
    """Abre la ventana interactiva con la Guía Paso a Paso de Uso del Software."""
    dlg = tk.Toplevel(parent)
    dlg.title("📖 Guía Paso a Paso de Uso del Software de Cinética de Citometría")
    dlg.geometry("960x660")
    dlg.transient(parent)

    hdr = tk.Frame(dlg, bg="#002B49", pady=6, padx=12)
    hdr.pack(fill="x")
    tk.Label(
        hdr,
        text="📖 GUÍA PRÁCTICA DE USO DEL SOFTWARE - FLUJO DE TRABAJO (PESTAÑAS 1 A 5)",
        bg="#002B49", fg="white", font=("DejaVu Sans", 10, "bold")
    ).pack(side="left")

    nb = ttk.Notebook(dlg)
    nb.pack(fill="both", expand=True, padx=10, pady=8)

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

    for title, body in steps:
        frm = ttk.Frame(nb, padding=12)
        nb.add(frm, text=title)
        txt = tk.Text(frm, wrap="word", font=("DejaVu Sans", 9), padx=10, pady=10, bg="#FDFEFE", relief="flat")
        txt.insert("1.0", body)
        txt.config(state="disabled")
        txt.pack(fill="both", expand=True)

    return dlg
