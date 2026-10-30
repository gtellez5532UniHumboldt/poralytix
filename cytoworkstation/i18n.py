"""
Módulo de Internacionalización (i18n) para la Estación de Trabajo de Citometría de Flujo.
Soporta Español (es) e Inglés (en) de forma interactiva y persistente.
"""

SUPPORTED_LANGUAGES = {
    'es': 'Español',
    'en': 'English'
}

DEFAULT_LANGUAGE = 'es'
_CURRENT_LANGUAGE = 'es'

TRANSLATIONS = {
    'es': {
        # App & Ventana Principal
        'app_window_title': "Estación de Trabajo de Cinética de Permeabilización Celular (Calceína-AM | .fcs)",
        'app_top_title': "🧬 WORKSTATION DE CINÉTICA DE PERMEABILIZACIÓN (.fcs) — Preanálisis, Cut-Offs por Muestra, Modelos y Auditoría",
        'status_ready': "Estado: Listo (Motor Nativo Activo)",
        'status_processing': "Estado: Procesando análisis...",
        'status_completed': "Estado: Análisis completado con éxito",
        'lang_label': "🌐 Idioma:",
        'lang_es': "🇪🇸 Español",
        'lang_en': "🇬🇧 English",

        # Botones de barra superior
        'btn_save_session': "💾 Guardar Sesión",
        'btn_load_session': "📂 Cargar Sesión",
        'btn_math_models': "📐 Modelos y Cálculos",
        'btn_user_guide': "📖 Guía de Uso",

        # Menú Principal
        'menu_session': "💾 Sesión y Configuración",
        'menu_session_save': "💾 Guardar / Exportar Sesión de Trabajo (.json)...    [Ctrl+S]",
        'menu_session_load': "📂 Cargar / Importar Sesión de Trabajo (.json)...     [Ctrl+O]",
        'menu_session_quicksave': "⚡ Guardado Rápido en Workspace (last_session_config.json)",
        'menu_session_quickload': "🔄 Restaurar Última Sesión Guardada (last_session_config.json)",
        'menu_session_export_cutoffs': "📋 Exportar Tabla de Cut-Offs por Muestra (.csv)...",
        'menu_session_exit': "🚪 Salir",

        'menu_models': "❓ Ayuda: Modelos Matemáticos y Cálculos",
        'menu_models_1': "📐 1. Modelos Farmacodinámicos y Cinéticos (Hill, Michaelis-Menten, Gompertz, Weibull, Bi-Exp)",
        'menu_models_2': "🔬 2. Gating (V1L/V1R), Cut-Off de Histograma (V2L/V2R) y Extracción Cinética Continua",
        'menu_models_3': "⏱️ 3. Hitos de Tiempo (T50 Efecto vs LT50, ΔT10-90), Entropía de Shannon, Hazard y Bimodalidad",
        'menu_models_4': "📊 4. Pruebas Estadísticas (Fisher, χ², h de Cohen, Wilson IC95%, Cochran-Armitage, Welch + FDR)",

        'menu_guide': "📖 Guía de Uso del Software",
        'menu_guide_steps': "📖 Guía Paso a Paso de Uso del Software    [F1]",
        'menu_guide_export': "📄 Exportar Guía de Uso y Manual Matemático Completo (.md)...",

        'menu_language': "🌐 Idioma / Language",

        # Pestañas del Notebook
        'tab_1_title': " 📁 1. Carga de Datos y Guía ",
        'tab_2_title': " 🔍 2. Pre-Análisis Dinámico (FITC / FSC / SSC) ",
        'tab_3_title': " 🧬 3. Panel Categórico (Afectados vs No Afectados y Ventanas) ",
        'tab_4_title': " 📊 4. Pre-Reporte, Modelos Cinéticos y Estadística ",
        'tab_5_title': " 📄 5. Reporte Final y Exportación (600 DPI / PDF) ",

        # Pestaña 1
        'tab1_guide_title': " 💡 Guía de Nomenclatura de Archivos .fcs para Autodetección Inmediata ",
        'tab1_guide_body': (
            "Formato recomendado:   [Compuesto]_[LineaCelular]_[Concentracion_uM]_[CC_opcional]_[Replica].fcs\n"
            "  • Control Basal Previo (-60 a 0 s): Incluir 'CC' o 'BASE'  ➔  Ej:  820_THP-1_CC_4.35_1.fcs   o   Sat1_PBMC_CC_5.70_R1.fcs\n"
            "  • Muestra Experimental (0 a 600 s): Sin 'CC' en el nombre  ➔  Ej:  820_THP-1_4.35_1.fcs      o   Sat1_THP-1_2.80_R2.fcs\n"
            "Nota: Si tus archivos tienen otros nombres, puedes editar cualquier celda con Doble Clic en la tabla inferior."
        ),
        'tab1_params_title': " Parámetros Celulares Globales y Selección del Canal de Fluorescencia ",
        'tab1_target_cell': "Célula Diana (Target):",
        'tab1_control_cell': "Célula Control Sano:",
        'tab1_density': "Densidad Celular (cél/mL):",
        'tab1_perm_channel': "Canal Permeabilidad:",
        'tab1_files_title': " Gestión de Archivos Binarios (.fcs) ",
        'tab1_btn_upload_files': "➕ Subir Archivos .fcs",
        'tab1_btn_upload_folder': "📁 Subir Carpeta Completa",
        'tab1_btn_demo': "🧪 Generar y Cargar Set Demo (.fcs)",
        'tab1_btn_batch_edit': "✏️ Editar en Lote Seleccionados",
        'tab1_btn_remove_selected': "🗑️ Quitar Seleccionados",
        'tab1_btn_clear_all': "🧹 Limpiar Tabla",
        'tab1_btn_run_preanalysis': "🔍 EJECUTAR PRE-ANÁLISIS Y VER GRÁFICAS ➔",
        'tab1_tbl_title': " Tabla de Muestras (Doble clic en cualquier celda para editar Compuesto, Célula, Rol CC, Concentración µM, Réplica o Exclusión t_start) ",

        # Columnas de tabla 1
        'col_filename': "Archivo .fcs",
        'col_compound': "Compuesto / Péptido",
        'col_cell_line': "Línea Celular",
        'col_sample_role': "Rol (Control CC vs Experimental)",
        'col_concentration': "Concentración (µM)",
        'col_replicate': "Réplica (R1..Rn)",
        'col_exclusion': "Exclusión Inicial t_start (s)",

        # Pestaña 2
        'tab2_ctrl_title': " 🎛️ Controles Globales de Gating Morfológico, Cut-Offs y Modos de Normalización ",
        'tab2_btn_recalc': "🔄 Recalcular y Actualizar Todo",
        'tab2_btn_sample_editor': "🔬 Editor Individual de Cut-Offs",
        'tab2_btn_apply_global': "🌐 Aplicar Valores Globales a Todas las Muestras",
        'tab2_btn_next': "Ir a Panel Categórico ➔",

        # Pestaña 3
        'tab3_btn_run_cat': "🧬 Ejecutar Análisis Categórico Completo",
        'tab3_btn_next': "Ir a Modelos Cinéticos ➔",

        # Pestaña 4
        'tab4_btn_fit_models': "📈 Ajustar Modelos No Lineales y Pruebas",
        'tab4_btn_next': "Ir a Reporte Final ➔",

        # Pestaña 5
        'tab5_btn_export_all': "🚀 EXPORTAR TODO (GRÁFICOS 600 DPI + TABLAS CSV + REPORTE PDF)",
        'tab5_btn_export_plots': "🖼️ Exportar Sólo Gráficos de Alta Resolución (600 DPI PNG)",
        'tab5_btn_export_tables': "📊 Exportar Sólo Tablas de Datos (.csv)",
        'tab5_btn_export_pdf': "📄 Generar y Exportar Reporte Ejecutivo PDF",

        # Diálogos y Alertas
        'dialog_confirm': "Confirmación",
        'dialog_info': "Información",
        'dialog_warning': "Advertencia",
        'dialog_error': "Error",
        'msg_no_files': "No se han cargado archivos .fcs en la tabla.",
        'msg_analysis_complete': "Análisis ejecutado exitosamente.",
        'msg_session_saved': "Sesión de trabajo guardada exitosamente en:\n{}",
        'msg_session_loaded': "Sesión de trabajo restaurada exitosamente.",
        'msg_lang_changed': "Idioma cambiado a Español.",

        # Editor de muestra
        'editor_title': "Editar Cut-Offs de Muestra: {}",
        'editor_group': " Parámetros Individuales: {} ",
        'editor_v1l': "Exclusión Detritos V1L (< log10):",
        'editor_v2l': "Cut-Off Afectados V2L/V2R (log10):",
        'editor_fsc': "Umbral Morfológico FSC-A Mín (>):",
        'editor_ssc': "Umbral Morfológico SSC-A Mín (>):",
        'editor_tlag': "Inicio Fase Activa T_lag (s):",
        'editor_tmax': "Fin Fase Activa / Meseta T_max (s):",
        'editor_btn_save': "💾 Guardar y Aplicar a esta Muestra",
        'editor_btn_cancel': "Cancelar",
    },
    'en': {
        # App & Main Window
        'app_window_title': "Cell Permeabilization Kinetics Workstation (Calcein-AM | .fcs)",
        'app_top_title': "🧬 CELL PERMEABILIZATION KINETICS WORKSTATION (.fcs) — Pre-Analysis, Sample Cut-Offs, Models & Audit",
        'status_ready': "Status: Ready (Native Engine Active)",
        'status_processing': "Status: Processing analysis...",
        'status_completed': "Status: Analysis completed successfully",
        'lang_label': "🌐 Language:",
        'lang_es': "🇪🇸 Español",
        'lang_en': "🇬🇧 English",

        # Top Bar Buttons
        'btn_save_session': "💾 Save Session",
        'btn_load_session': "📂 Load Session",
        'btn_math_models': "📐 Models & Calculations",
        'btn_user_guide': "📖 User Guide",

        # Main Menu
        'menu_session': "💾 Session & Settings",
        'menu_session_save': "💾 Save / Export Work Session (.json)...    [Ctrl+S]",
        'menu_session_load': "📂 Load / Import Work Session (.json)...     [Ctrl+O]",
        'menu_session_quicksave': "⚡ Quick Save in Workspace (last_session_config.json)",
        'menu_session_quickload': "🔄 Restore Last Saved Session (last_session_config.json)",
        'menu_session_export_cutoffs': "📋 Export Sample Cut-Offs Table (.csv)...",
        'menu_session_exit': "🚪 Exit",

        'menu_models': "❓ Help: Mathematical Models & Calculations",
        'menu_models_1': "📐 1. Pharmacodynamic & Kinetic Models (Hill, Michaelis-Menten, Gompertz, Weibull, Bi-Exp)",
        'menu_models_2': "🔬 2. Gating (V1L/V1R), Histogram Cut-Off (V2L/V2R) & Continuous Kinetic Extraction",
        'menu_models_3': "⏱️ 3. Time Milestones (T50 Effect vs LT50, ΔT10-90), Shannon Entropy, Hazard & Bimodality",
        'menu_models_4': "📊 4. Statistical Tests (Fisher, χ², Cohen's h, Wilson 95% CI, Cochran-Armitage, Welch + FDR)",

        'menu_guide': "📖 Software User Guide",
        'menu_guide_steps': "📖 Step-by-Step Software User Guide    [F1]",
        'menu_guide_export': "📄 Export Complete User Guide & Math Manual (.md)...",

        'menu_language': "🌐 Language / Idioma",

        # Notebook Tabs
        'tab_1_title': " 📁 1. Data Ingestion & Guide ",
        'tab_2_title': " 🔍 2. Dynamic Pre-Analysis (FITC / FSC / SSC) ",
        'tab_3_title': " 🧬 3. Categorical Panel (Affected vs Non-Affected & Windows) ",
        'tab_4_title': " 📊 4. Pre-Report, Kinetic Models & Statistics ",
        'tab_5_title': " 📄 5. Final Report & Export (600 DPI / PDF) ",

        # Tab 1
        'tab1_guide_title': " 💡 .fcs File Naming Guide for Immediate Auto-Detection ",
        'tab1_guide_body': (
            "Recommended format:   [Compound]_[CellLine]_[Concentration_uM]_[CC_optional]_[Replicate].fcs\n"
            "  • Baseline Control (-60 to 0 s): Include 'CC' or 'BASE'  ➔  Ex:  820_THP-1_CC_4.35_1.fcs   or   Sat1_PBMC_CC_5.70_R1.fcs\n"
            "  • Experimental Sample (0 to 600 s): Without 'CC' in name  ➔  Ex:  820_THP-1_4.35_1.fcs      or   Sat1_THP-1_2.80_R2.fcs\n"
            "Note: If your files use different names, you can edit any cell by Double-Clicking in the table below."
        ),
        'tab1_params_title': " Global Cellular Parameters and Fluorescence Channel Selection ",
        'tab1_target_cell': "Target Cell Line:",
        'tab1_control_cell': "Healthy Control Cell:",
        'tab1_density': "Cell Density (cells/mL):",
        'tab1_perm_channel': "Permeability Channel:",
        'tab1_files_title': " Binary File Management (.fcs) ",
        'tab1_btn_upload_files': "➕ Upload .fcs Files",
        'tab1_btn_upload_folder': "📁 Upload Full Folder",
        'tab1_btn_demo': "🧪 Generate & Load Demo Set (.fcs)",
        'tab1_btn_batch_edit': "✏️ Batch Edit Selected",
        'tab1_btn_remove_selected': "🗑️ Remove Selected",
        'tab1_btn_clear_all': "🧹 Clear Table",
        'tab1_btn_run_preanalysis': "🔍 RUN PRE-ANALYSIS & VIEW PLOTS ➔",
        'tab1_tbl_title': " Sample Table (Double-click any cell to edit Compound, Cell, Role CC, Concentration µM, Replicate or Exclusion t_start) ",

        # Tab 1 table columns
        'col_filename': ".fcs File",
        'col_compound': "Compound / Peptide",
        'col_cell_line': "Cell Line",
        'col_sample_role': "Role (Control CC vs Experimental)",
        'col_concentration': "Concentration (µM)",
        'col_replicate': "Replicate (R1..Rn)",
        'col_exclusion': "Initial Exclusion t_start (s)",

        # Tab 2
        'tab2_ctrl_title': " 🎛️ Global Morphological Gating Controls, Cut-Offs and Normalization Modes ",
        'tab2_btn_recalc': "🔄 Recalculate and Update All",
        'tab2_btn_sample_editor': "🔬 Individual Sample Cut-Offs Editor",
        'tab2_btn_apply_global': "🌐 Apply Global Values to All Samples",
        'tab2_btn_next': "Go to Categorical Panel ➔",

        # Tab 3
        'tab3_btn_run_cat': "🧬 Run Full Categorical Analysis",
        'tab3_btn_next': "Go to Kinetic Models ➔",

        # Tab 4
        'tab4_btn_fit_models': "📈 Fit Non-Linear Models & Hypotheses",
        'tab4_btn_next': "Go to Final Report ➔",

        # Tab 5
        'tab5_btn_export_all': "🚀 EXPORT ALL (600 DPI PLOTS + CSV TABLES + PDF REPORT)",
        'tab5_btn_export_plots': "🖼️ Export Only High-Resolution Plots (600 DPI PNG)",
        'tab5_btn_export_tables': "📊 Export Only Data Tables (.csv)",
        'tab5_btn_export_pdf': "📄 Generate & Export Executive PDF Report",

        # Dialogs and Alerts
        'dialog_confirm': "Confirmation",
        'dialog_info': "Information",
        'dialog_warning': "Warning",
        'dialog_error': "Error",
        'msg_no_files': "No .fcs files loaded in the table.",
        'msg_analysis_complete': "Analysis executed successfully.",
        'msg_session_saved': "Work session successfully saved to:\n{}",
        'msg_session_loaded': "Work session successfully restored.",
        'msg_lang_changed': "Language switched to English.",

        # Sample editor
        'editor_title': "Edit Sample Cut-Offs: {}",
        'editor_group': " Individual Parameters: {} ",
        'editor_v1l': "Debris Exclusion V1L (< log10):",
        'editor_v2l': "Affected Cut-Off V2L/V2R (log10):",
        'editor_fsc': "Morphological Threshold FSC-A Min (>):",
        'editor_ssc': "Morphological Threshold SSC-A Min (>):",
        'editor_tlag': "Active Phase Start T_lag (s):",
        'editor_tmax': "Active Phase End / Plateau T_max (s):",
        'editor_btn_save': "💾 Save and Apply to this Sample",
        'editor_btn_cancel': "Cancel",
    }
}


def get_current_language():
    """Devuelve el código del idioma activo ('es' o 'en')."""
    return _CURRENT_LANGUAGE


def set_current_language(lang):
    """Establece el idioma activo ('es' o 'en')."""
    global _CURRENT_LANGUAGE
    if lang in SUPPORTED_LANGUAGES:
        _CURRENT_LANGUAGE = lang
    return _CURRENT_LANGUAGE


def t(key, lang=None, **kwargs):
    """
    Traduce una clave al idioma especificado o al actual.
    
    Parameters
    ----------
    key : str
        Clave de traducción en TRANSLATIONS.
    lang : str, optional
        Código de idioma ('es' o 'en'). Si es None, usa el idioma global actual.
    **kwargs
        Argumentos para formatear strings dinámicos con .format(**kwargs).
        
    Returns
    -------
    str
        Cadena traducida o la clave original si no se encuentra.
    """
    target_lang = lang or _CURRENT_LANGUAGE
    dict_lang = TRANSLATIONS.get(target_lang, TRANSLATIONS['es'])
    text = dict_lang.get(key)
    if text is None:
        # Fallback a español
        text = TRANSLATIONS['es'].get(key, key)
    
    if kwargs and isinstance(text, str):
        try:
            return text.format(**kwargs)
        except Exception:
            return text
    return text
