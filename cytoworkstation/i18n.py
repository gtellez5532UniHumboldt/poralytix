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

        # Pestañas Principales del Notebook
        'tab_1_title': " 📁 1. Carga de Datos y Guía ",
        'tab_2_title': " 🔍 2. Pre-Análisis Dinámico (FITC / FSC / SSC) ",
        'tab_3_title': " 🧬 3. Panel Categórico (Afectados vs No Afectados y Ventanas) ",
        'tab_4_title': " 📊 4. Pre-Reporte, Modelos Cinéticos y Estadística ",
        'tab_5_title': " 📄 5. Reporte Final y Exportación (600 DPI / PDF) ",

        # Pestaña 1: Carga y Metadatos
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

        # Columnas de Tabla 1
        'col_filename': "Archivo .fcs",
        'col_compound': "Compuesto / Péptido",
        'col_cell_line': "Línea Celular",
        'col_sample_role': "Rol (Control CC vs Experimental)",
        'col_concentration': "Concentración (µM)",
        'col_replicate': "Réplica (R1..Rn)",
        'col_exclusion': "Exclusión Inicial t_start (s)",

        # Pestaña 2: Pre-Análisis Dinámico
        'tab2_ctrl_title': " 🎛️ Controles Globales de Gating Morfológico, Cut-Offs y Modos de Normalización ",
        'tab2_btn_recalc': "🔄 Recalcular y Actualizar Todo",
        'tab2_btn_sample_editor': "🔬 Editor Individual de Cut-Offs",
        'tab2_btn_apply_global': "🌐 Aplicar Valores Globales a Todas las Muestras",
        'tab2_btn_next': "Ir a Panel Categórico ➔",
        'tab2_gate_title': " 1. Canal, Prefiltros (FSC/SSC/V1L) y Cut-Off V2L/V2R ",
        'tab2_perm_channel': "Canal Permeabilidad:",
        'tab2_fsc_min': "FSC-A Mín (>):",
        'tab2_fsc_max': "FSC-A Máx:",
        'tab2_ssc_min': "SSC-A Mín (>):",
        'tab2_ssc_max': "SSC-A Máx:",
        'tab2_no_limit': "Sin límite",
        'tab2_v1l_excl': "Exclusión V1L (<log10):",
        'tab2_v2l_cutoff': "Cut-Off V2L/V2R (log10):",
        'tab2_bin_size': "Bin Temporal (s):",
        'tab2_hist_ymode': "Eje Y Hist.:",
        'tab2_cutoff_scope': "Alcance Filtro/Cut-Off:",
        'tab2_scope_global': "Global (Todas las muestras)",
        'tab2_scope_selected': "Discriminado (Muestras seleccionadas)",
        'tab2_btn_apply_all': "🌐 Aplicar a TODAS",
        'tab2_btn_apply_sel': "🎯 Aplicar a SELECCIONADAS",
        'tab2_btn_auto_cc': "Auto CC",
        'tab2_plot1_mode': "Modo Gráfica 1:",
        'tab2_mode_cont': "Continua (Mediana FITC %)",
        'tab2_mode_cat': "Categórica (Afectados V2L % por Cut-Off)",
        'tab2_mode_dual': "Dual (Continua + Afectados V2L %)",
        'tab2_chk_scatter_norm': "Mostrar FSC/SSC Normalizados (% vs Basal)",
        'tab2_gate_stats': "Canal: {} | V1R Retenidos: {} | V1L Omitidos: {}",
        'tab2_sel_title': " 2. Selección Dinámica de Muestras a Graficar ",
        'tab2_compounds': "Compuestos:",
        'tab2_cell_lines': "Líneas:",
        'tab2_doses': "Dosis (µM):",
        'tab2_btn_sel_all': "Seleccionar Todas",
        'tab2_btn_refresh': "Refrescar Gráficas",
        'tab2_thresh_title': " 3. Ajuste de T_lag y T_max por Curva (Arrastrable) ",
        'tab2_active_curve': "Curva Activa:",
        'tab2_radio_auto': "Automático",
        'tab2_radio_manual': "Manual / Arrastrar",
        'tab2_btn_auto_all': "⚡ Auto (Todas)",
        'tab2_btn_auto_active': "⚡ Auto (Activa)",
        'tab2_tlag_lbl': "T_lag (s) [Azul :]:",
        'tab2_tmax_lbl': "T_max (s) [Rojo -.]:",
        'tab2_active_metrics': "T_lag: {} s | T_max: {} s\nV_active: {} %/s | V2L Máx: {} % | LT50: {} s",
        'tab2_cut_tbl_title': " 4. Tabla Interactiva de Cut-Offs por Muestra (Doble Clic en Celda = Editar) ",
        'tab2_btn_edit_sample': "✏️ Editar Muestra Sel.",
        'tab2_btn_export_csv': "📋 Exportar Cut-Offs (.csv)",
        'tab2_btn_go_cat': "🧬 Panel Categórico ➔",
        'tab2_col_sample': "Muestra (Comp | Cél | µM)",
        'tab2_col_scope': "Alcance",
        'tab2_col_v1l': "V1L <log",
        'tab2_col_v2l': "CutOff V2L",
        'tab2_col_fsc': "FSC Mín",
        'tab2_col_ssc': "SSC Mín",
        'tab2_col_tlag': "T_lag(s)",
        'tab2_col_tmax': "T_max(s)",
        'tab2_col_mode': "Modo T",
        'tab2_col_v2l_max': "%V2L Máx",

        # Pestaña 3: Panel Categórico
        'tab3_btn_run_cat': "🧬 Ejecutar Análisis Categórico Completo",
        'tab3_btn_next': "Ir a Modelos Cinéticos ➔",
        'tab3_ctrl_title': " Configuración del Panel Categórico: Proporciones Afectados (V2L) vs. No Afectados (V2R), Ventanas de Tiempo, T50, T_efecto_máx y Pruebas de Afectación ",
        'tab3_compound': "Compuesto:",
        'tab3_cell_line': "Línea Celular:",
        'tab3_curve': "Curva Destacada:",
        'tab3_vis_mode': "Visualización Gráfica:",
        'tab3_opt_all_compounds': "Todos",
        'tab3_opt_all_cells': "Todas",
        'tab3_vis_opt1': "Balance Temporal: % Afectados (V2L) vs % No Afectados (V2R) y Hitos",
        'tab3_vis_opt2': "Comparativa de Hitos: T_onset, T50 del Efecto y T_efecto_máx + % Máx/Mín por Dosis",
        'tab3_vis_opt3': "Mapa de Calor Espaciotemporal: % Eventos Afectados (V2L) por Compuesto, Dosis y Tiempo",
        'tab3_vis_opt4': "Velocidad Instantánea de Conversión d(%V2L)/dt y Entropía de Shannon H(t)",
        'tab3_btn_recalc': "🔄 Recalcular Panel Categórico",
        'tab3_btn_export_report': "📄 Exportar en Reporte Final ➔",
        'tab3_subtab_milestones': " 1. Hitos del Efecto (% Afectados Máx/Mín, T_efecto_máx, T50 Efecto, LT50 y Definición de Afectación) ",
        'tab3_subtab_windows': " 2. Proporciones por Ventana de Tiempo (Afectados V2L vs No Afectados V2R, IC95%, Odds y Fisher/Chi²) ",
        'tab3_subtab_tests': " 3. Pruebas Estadísticas (Tendencia Cochran-Armitage, Log-Rank, Isocronas EC50(t) y Contrastes) ",
        'tab3_subtab_insights': " 4. Interpretación Biofísica Automática y Análisis Sugeridos del Fenómeno ",

        # Pestaña 4: Pre-Reporte y Modelos
        'tab4_btn_fit_models': "📈 Ajustar Modelos No Lineales y Pruebas",
        'tab4_ctrl_title': " Configuración de Analizadores Cinéticos / Farmacodinámicos y Diseño de Comparación Estadística ",
        'tab4_analyzer': "Analizador Cinético / Modelo:",
        'tab4_stat_mode': "Modo de Comparación Estadística:",
        'tab4_fdr': "Corrección Múltiple:",
        'tab4_bootstrap': "Réplicas Bootstrap (IC 95%):",
        'tab4_btn_recalc': "🔄 Recalcular Modelos y Tablas",
        'tab4_btn_explain': "❓ Explicación del Modelo y Ecuaciones",
        'tab4_subtab_t0': "Tabla 0: Prerregistro Pre-Análisis",
        'tab4_subtab_t1': "Tabla 1: Tasa 0-120s",
        'tab4_subtab_t2': "Tabla 2: Tasa 0-600s",
        'tab4_subtab_t3': "Tabla 3: V_active [T_lag-T_max]",
        'tab4_subtab_tcat': "Tabla Cat: Cut-Off V2L%, LT50, Entropía y Modo",
        'tab4_subtab_t4': "Tabla 4: Selectividad y Odds Ratio",
        'tab4_subtab_t5': "Tabla 5: Contrastes Dinámicos / Dosis",
        'tab4_subtab_t6': "Tabla 6: Parámetros del Modelo (Hill / Otros)",
        'tab4_subtab_ttimebin': "Tabla Bin: Datos por Bin de Tiempo (Diana vs Control)",
        'tab4_sat_title': " 📐 Panel de Previsualización: Gráfica de la Función, Ecuación de los Modelos y R² del Grupo de Muestras ",
        'tab4_sample_group': "Grupo de Muestras:",
        'tab4_all_groups': "Todos los Grupos",
        'tab4_view': "Vista:",
        'tab4_view_selected': "Modelo Seleccionado (Por Grupo)",
        'tab4_view_compare': "Comparar Todos los Modelos en el Grupo",
        'tab4_eq_title': " Ecuación Parametrizada de la Función y Bondad de Ajuste (R² del Grupo de Muestras) ",
        'tab4_col_group': "Grupo de Muestras",
        'tab4_col_model': "Modelo",
        'tab4_col_eq': "Ecuación Evaluada de la Función f(x)",
        'tab4_col_r2_group': "R² Grupo (Medias)",
        'tab4_col_r2_replicas': "R² Réplicas (Todas)",
        'tab4_col_r2_adj': "R² Ajustado",
        'tab4_col_rmse': "RMSE",

        # Pestaña 5: Reporte Final y Exportación
        'tab5_btn_export_all': "🚀 EXPORTAR TODO (GRÁFICOS 600 DPI + TABLAS CSV + REPORTE PDF)",
        'tab5_btn_export_plots': "🖼️ Exportar Sólo Gráficos de Alta Resolución (600 DPI PNG)",
        'tab5_btn_export_tables': "📊 Exportar Sólo Tablas de Datos (.csv)",
        'tab5_btn_export_pdf': "📄 Generar y Exportar Reporte Ejecutivo PDF",
        'tab5_top_title': " Configuración de Exportación del Reporte Auditable por Categorías (600 DPI PNG + CSV/TSV + PDF) ",
        'tab5_outdir': "Carpeta de Destino:",
        'tab5_btn_browse': "Examinar...",
        'tab5_rep_title': "Título del Reporte PDF:",
        'tab5_rep_author': "Autor / Laboratorio:",
        'tab5_btn_export_bundle': "🚀 GENERAR Y EXPORTAR REPORTE FINAL AUDITABLE (PDF CON IMÁGENES + 600 DPI + CSV/TSV)",
        'tab5_preview_title': " Bitácora Ejecutiva del Reporte Final (El programa permanece abierto tras exportar) ",

        # Diálogos, Edición en Lote y Alertas
        'dialog_confirm': "Confirmación",
        'dialog_info': "Información",
        'dialog_warning': "Advertencia",
        'dialog_error': "Error",
        'dialog_select_row': "Selecciona primero una o más filas en la tabla.",
        'dialog_no_fcs_folder': "No se encontraron archivos .fcs en la carpeta seleccionada.",
        'dialog_select_outdir': "Seleccionar Carpeta de Salida",
        'dialog_select_files': "Seleccionar archivos .fcs",
        'dialog_select_folder': "Seleccionar carpeta con archivos .fcs",
        'msg_no_files': "No se han cargado archivos .fcs en la tabla.",
        'msg_analysis_complete': "Análisis ejecutado exitosamente.",
        'msg_session_saved': "Sesión de trabajo guardada exitosamente en:\n{}",
        'msg_session_loaded': "Sesión de trabajo restaurada exitosamente.",
        'msg_lang_changed': "Idioma cambiado a Español.",
        'batch_edit_title': "Edición en Lote ({} archivos)",
        'batch_edit_compound': "Compuesto:",
        'batch_edit_cell': "Línea Celular:",
        'batch_edit_role': "Rol:",
        'batch_edit_dose': "Dosis (µM):",
        'batch_edit_replicate': "Réplica:",
        'batch_edit_exclusion': "Exclusión t_start (s):",
        'batch_edit_btn_apply': "Aplicar Cambios",

        # Editor individual de muestra
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

        # Textos de Lienzos (Canvas)
        'canvas_hist_title': "Histograma Dinámico (Muestras Seleccionadas) — Exclusión V1L < {:.2f} | Cut-Off V2L/V2R = {:.2f} [{}]",
        'canvas_hist_v1l': "V1L: Excluidos",
        'canvas_hist_v2l': "◄ V2L: Afectados",
        'canvas_hist_v2r': "V2R: Sanos (Intactos) ►",
        'canvas_hist_xlabel': "log10(Fluorescencia Canal Seleccionado)",
        'canvas_hist_ylabel_pct': "Eventos (%)",
        'canvas_hist_ylabel_cnt': "N° Eventos",
        'canvas_hist_legend': "Muestras en Histograma ({}):",
        'canvas_hist_more': "+ {} curvas adicionales...",
        'canvas_kin_time': "Tiempo (segundos)",
        'canvas_kin_legend': "Leyenda de Curvas:",
        'canvas_kin_p1_cat': "1. Proporción de Eventos Afectados (V2L % < Cut-Off) vs. Tiempo — {}",
        'canvas_kin_p1_dual': "1. Superposición Dual: Continua (Sólida) vs. Afectados V2L % (Punteada) — {}",
        'canvas_kin_p1_cont': "1. Cinética de Permeabilización (Mediana Continua) vs. Tiempo — {}",
        'canvas_kin_p2_fsc': "2. Dispersión Frontal (FSC vs. Tiempo — Tamaño Celular / Tumefacción)",
        'canvas_kin_p3_ssc': "3. Dispersión Lateral (SSC vs. Tiempo — Granularidad / Complejidad Interna)",
        'canvas_cat_ylabel': "Proporción (%)",
        'canvas_cat_xlabel': "Ventana de Tiempo (segundos)",
        'canvas_cat_curve': "Curva: {}",
        'canvas_cat_v2l': "% Afectados (V2L < Cut-Off)",
        'canvas_cat_v2r': "% No Afectados (V2R >= Cut-Off)",
        'canvas_cat_milestones_a': "A. Hitos Temporales: T_onset (Azul), T50 Efecto (Naranja), T_máx (Rojo)",
        'canvas_cat_milestones_b': "B. Proporción en Efecto Máximo: % Afectados V2L (Rojo) vs % No Afectados V2R (Verde)",
        'canvas_cat_heatmap_title': "Mapa de Calor: % de Eventos Afectados (V2L) por Ventana de Tiempo",
        'canvas_cat_rate_title': "Velocidad Instantánea de Conversión d(%V2L)/dt (%/s)",
        'canvas_cat_rate_legend': "Velocidad d(%V2L)/dt:",
        'canvas_sat_fn': "Función del Modelo: {}",
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

        # Notebook Main Tabs
        'tab_1_title': " 📁 1. Data Ingestion & Guide ",
        'tab_2_title': " 🔍 2. Dynamic Pre-Analysis (FITC / FSC / SSC) ",
        'tab_3_title': " 🧬 3. Categorical Panel (Affected vs Non-Affected & Windows) ",
        'tab_4_title': " 📊 4. Pre-Report, Kinetic Models & Statistics ",
        'tab_5_title': " 📄 5. Final Report & Export (600 DPI / PDF) ",

        # Tab 1: Ingestion & Metadata
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

        # Tab 2: Dynamic Pre-Analysis
        'tab2_ctrl_title': " 🎛️ Global Morphological Gating Controls, Cut-Offs and Normalization Modes ",
        'tab2_btn_recalc': "🔄 Recalculate and Update All",
        'tab2_btn_sample_editor': "🔬 Individual Sample Cut-Offs Editor",
        'tab2_btn_apply_global': "🌐 Apply Global Values to All Samples",
        'tab2_btn_next': "Go to Categorical Panel ➔",
        'tab2_gate_title': " 1. Channel, Pre-Filters (FSC/SSC/V1L) & Cut-Off V2L/V2R ",
        'tab2_perm_channel': "Permeability Channel:",
        'tab2_fsc_min': "FSC-A Min (>):",
        'tab2_fsc_max': "FSC-A Max:",
        'tab2_ssc_min': "SSC-A Min (>):",
        'tab2_ssc_max': "SSC-A Max:",
        'tab2_no_limit': "No limit",
        'tab2_v1l_excl': "V1L Exclusion (<log10):",
        'tab2_v2l_cutoff': "Cut-Off V2L/V2R (log10):",
        'tab2_bin_size': "Time Bin (s):",
        'tab2_hist_ymode': "Hist Y-Axis:",
        'tab2_cutoff_scope': "Filter/Cut-Off Scope:",
        'tab2_scope_global': "Global (All samples)",
        'tab2_scope_selected': "Discriminated (Selected samples)",
        'tab2_btn_apply_all': "🌐 Apply to ALL",
        'tab2_btn_apply_sel': "🎯 Apply to SELECTED",
        'tab2_btn_auto_cc': "Auto CC",
        'tab2_plot1_mode': "Plot 1 Mode:",
        'tab2_mode_cont': "Continuous (FITC Median %)",
        'tab2_mode_cat': "Categorical (Affected V2L % by Cut-Off)",
        'tab2_mode_dual': "Dual (Continuous + Affected V2L %)",
        'tab2_chk_scatter_norm': "Show Normalized FSC/SSC (% vs Baseline)",
        'tab2_gate_stats': "Channel: {} | V1R Retained: {} | V1L Omitted: {}",
        'tab2_sel_title': " 2. Dynamic Selection of Samples to Plot ",
        'tab2_compounds': "Compounds:",
        'tab2_cell_lines': "Cell Lines:",
        'tab2_doses': "Doses (µM):",
        'tab2_btn_sel_all': "Select All",
        'tab2_btn_refresh': "Refresh Plots",
        'tab2_thresh_title': " 3. Adjust T_lag and T_max per Curve (Draggable) ",
        'tab2_active_curve': "Active Curve:",
        'tab2_radio_auto': "Automatic",
        'tab2_radio_manual': "Manual / Drag",
        'tab2_btn_auto_all': "⚡ Auto (All)",
        'tab2_btn_auto_active': "⚡ Auto (Active)",
        'tab2_tlag_lbl': "T_lag (s) [Blue :]:",
        'tab2_tmax_lbl': "T_max (s) [Red -.]:",
        'tab2_active_metrics': "T_lag: {} s | T_max: {} s\nV_active: {} %/s | V2L Max: {} % | LT50: {} s",
        'tab2_cut_tbl_title': " 4. Interactive Sample Cut-Offs Table (Double-Click Cell = Edit) ",
        'tab2_btn_edit_sample': "✏️ Edit Sel. Sample",
        'tab2_btn_export_csv': "📋 Export Cut-Offs (.csv)",
        'tab2_btn_go_cat': "🧬 Categorical Panel ➔",
        'tab2_col_sample': "Sample (Comp | Cell | µM)",
        'tab2_col_scope': "Scope",
        'tab2_col_v1l': "V1L <log",
        'tab2_col_v2l': "CutOff V2L",
        'tab2_col_fsc': "FSC Min",
        'tab2_col_ssc': "SSC Min",
        'tab2_col_tlag': "T_lag(s)",
        'tab2_col_tmax': "T_max(s)",
        'tab2_col_mode': "T Mode",
        'tab2_col_v2l_max': "%V2L Max",

        # Tab 3: Categorical Panel
        'tab3_btn_run_cat': "🧬 Run Full Categorical Analysis",
        'tab3_btn_next': "Go to Kinetic Models ➔",
        'tab3_ctrl_title': " Categorical Panel Configuration: Proportions Affected (V2L) vs Non-Affected (V2R), Time Windows, T50, T_effect_max & Affectation Tests ",
        'tab3_compound': "Compound:",
        'tab3_cell_line': "Cell Line:",
        'tab3_curve': "Featured Curve:",
        'tab3_vis_mode': "Graphical Visualization:",
        'tab3_opt_all_compounds': "All",
        'tab3_opt_all_cells': "All",
        'tab3_vis_opt1': "Temporal Balance: % Affected (V2L) vs % Non-Affected (V2R) & Milestones",
        'tab3_vis_opt2': "Milestones Comparison: T_onset, Effect T50 & T_effect_max + Max/Min % by Dose",
        'tab3_vis_opt3': "Spatiotemporal Heatmap: % Affected Events (V2L) by Compound, Dose & Time",
        'tab3_vis_opt4': "Instantaneous Conversion Rate d(%V2L)/dt & Shannon Entropy H(t)",
        'tab3_btn_recalc': "🔄 Recalculate Categorical Panel",
        'tab3_btn_export_report': "📄 Export to Final Report ➔",
        'tab3_subtab_milestones': " 1. Effect Milestones (Max/Min % Affected, T_effect_max, Effect T50, LT50 & Classification) ",
        'tab3_subtab_windows': " 2. Proportions by Time Window (Affected V2L vs Non-Affected V2R, 95% CI, Odds & Fisher/Chi²) ",
        'tab3_subtab_tests': " 3. Statistical Tests (Cochran-Armitage Trend, Log-Rank, EC50(t) Isochrones & Contrasts) ",
        'tab3_subtab_insights': " 4. Automated Biophysical Interpretation & Suggested Phenomenon Analysis ",

        # Tab 4: Pre-Report & Models
        'tab4_btn_fit_models': "📈 Fit Non-Linear Models & Hypotheses",
        'tab4_ctrl_title': " Configuration of Kinetic / Pharmacodynamic Analyzers and Statistical Comparison Design ",
        'tab4_analyzer': "Kinetic Analyzer / Model:",
        'tab4_stat_mode': "Statistical Comparison Mode:",
        'tab4_fdr': "Multiple Testing Correction:",
        'tab4_bootstrap': "Bootstrap Replicates (95% CI):",
        'tab4_btn_recalc': "🔄 Recalculate Models and Tables",
        'tab4_btn_explain': "❓ Model & Equations Explanation",
        'tab4_subtab_t0': "Table 0: Pre-Analysis Preregistration",
        'tab4_subtab_t1': "Table 1: Rate 0-120s",
        'tab4_subtab_t2': "Table 2: Rate 0-600s",
        'tab4_subtab_t3': "Table 3: V_active [T_lag-T_max]",
        'tab4_subtab_tcat': "Table Cat: Cut-Off V2L%, LT50, Entropy & Mode",
        'tab4_subtab_t4': "Table 4: Selectivity & Odds Ratio",
        'tab4_subtab_t5': "Table 5: Dynamic Contrasts / Dose",
        'tab4_subtab_t6': "Table 6: Model Parameters (Hill / Others)",
        'tab4_subtab_ttimebin': "Table Bin: Time Bin Data (Target vs Control)",
        'tab4_sat_title': " 📐 Preview Panel: Function Plot, Model Equation & Sample Group R² ",
        'tab4_sample_group': "Sample Group:",
        'tab4_all_groups': "All Groups",
        'tab4_view': "View:",
        'tab4_view_selected': "Selected Model (By Group)",
        'tab4_view_compare': "Compare All Models in Group",
        'tab4_eq_title': " Parameterized Function Equation & Goodness of Fit (Sample Group R²) ",
        'tab4_col_group': "Sample Group",
        'tab4_col_model': "Model",
        'tab4_col_eq': "Evaluated Function Equation f(x)",
        'tab4_col_r2_group': "R² Group (Means)",
        'tab4_col_r2_replicas': "R² Replicates (All)",
        'tab4_col_r2_adj': "Adjusted R²",
        'tab4_col_rmse': "RMSE",

        # Tab 5: Final Report & Export
        'tab5_btn_export_all': "🚀 EXPORT ALL (600 DPI PLOTS + CSV TABLES + PDF REPORT)",
        'tab5_btn_export_plots': "🖼️ Export Only High-Resolution Plots (600 DPI PNG)",
        'tab5_btn_export_tables': "📊 Export Only Data Tables (.csv)",
        'tab5_btn_export_pdf': "📄 Generate & Export Executive PDF Report",
        'tab5_top_title': " Audit Report Export Configuration by Categories (600 DPI PNG + CSV/TSV + PDF) ",
        'tab5_outdir': "Output Directory:",
        'tab5_btn_browse': "Browse...",
        'tab5_rep_title': "PDF Report Title:",
        'tab5_rep_author': "Author / Laboratory:",
        'tab5_btn_export_bundle': "🚀 GENERATE & EXPORT AUDIT FINAL REPORT (PDF WITH PLOTS + 600 DPI + CSV/TSV)",
        'tab5_preview_title': " Executive Audit Log of Final Report (The program remains open after export) ",

        # Dialogs, Batch Edit & Alerts
        'dialog_confirm': "Confirmation",
        'dialog_info': "Information",
        'dialog_warning': "Warning",
        'dialog_error': "Error",
        'dialog_select_row': "Please select one or more rows in the table first.",
        'dialog_no_fcs_folder': "No .fcs files found in selected folder.",
        'dialog_select_outdir': "Select Output Directory",
        'dialog_select_files': "Select .fcs files",
        'dialog_select_folder': "Select folder with .fcs files",
        'msg_no_files': "No .fcs files loaded in the table.",
        'msg_analysis_complete': "Analysis executed successfully.",
        'msg_session_saved': "Work session successfully saved to:\n{}",
        'msg_session_loaded': "Work session successfully restored.",
        'msg_lang_changed': "Language switched to English.",
        'batch_edit_title': "Batch Edit ({} files)",
        'batch_edit_compound': "Compound:",
        'batch_edit_cell': "Cell Line:",
        'batch_edit_role': "Role:",
        'batch_edit_dose': "Dose (µM):",
        'batch_edit_replicate': "Replicate:",
        'batch_edit_exclusion': "Exclusion t_start (s):",
        'batch_edit_btn_apply': "Apply Changes",

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

        # Canvas Strings
        'canvas_hist_title': "Dynamic Histogram (Selected Samples) — V1L Exclusion < {:.2f} | Cut-Off V2L/V2R = {:.2f} [{}]",
        'canvas_hist_v1l': "V1L: Excluded",
        'canvas_hist_v2l': "◄ V2L: Affected",
        'canvas_hist_v2r': "V2R: Healthy (Intact) ►",
        'canvas_hist_xlabel': "log10(Selected Channel Fluorescence)",
        'canvas_hist_ylabel_pct': "Events (%)",
        'canvas_hist_ylabel_cnt': "Event Count",
        'canvas_hist_legend': "Samples in Histogram ({}):",
        'canvas_hist_more': "+ {} additional curves...",
        'canvas_kin_time': "Time (seconds)",
        'canvas_kin_legend': "Curve Legend:",
        'canvas_kin_p1_cat': "1. Proportion of Affected Events (V2L % < Cut-Off) vs. Time — {}",
        'canvas_kin_p1_dual': "1. Dual Overlay: Continuous (Solid) vs. Affected V2L % (Dotted) — {}",
        'canvas_kin_p1_cont': "1. Permeabilization Kinetics (Continuous Median) vs. Time — {}",
        'canvas_kin_p2_fsc': "2. Forward Scatter (FSC vs. Time — Cell Size / Swelling)",
        'canvas_kin_p3_ssc': "3. Side Scatter (SSC vs. Time — Granularity / Internal Complexity)",
        'canvas_cat_ylabel': "Proportion (%)",
        'canvas_cat_xlabel': "Time Window (seconds)",
        'canvas_cat_curve': "Curve: {}",
        'canvas_cat_v2l': "% Affected (V2L < Cut-Off)",
        'canvas_cat_v2r': "% Non-Affected (V2R >= Cut-Off)",
        'canvas_cat_milestones_a': "A. Temporal Milestones: T_onset (Blue), Effect T50 (Orange), T_max (Red)",
        'canvas_cat_milestones_b': "B. Proportion at Maximum Effect: % Affected V2L (Red) vs % Non-Affected V2R (Green)",
        'canvas_cat_heatmap_title': "Heatmap: % Affected Events (V2L) by Time Window",
        'canvas_cat_rate_title': "Instantaneous Conversion Rate d(%V2L)/dt (%/s)",
        'canvas_cat_rate_legend': "Rate d(%V2L)/dt:",
        'canvas_sat_fn': "Model Function: {}",
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
        text = TRANSLATIONS['es'].get(key, key)
    
    if kwargs and isinstance(text, str):
        try:
            return text.format(**kwargs)
        except Exception:
            return text
    return text
