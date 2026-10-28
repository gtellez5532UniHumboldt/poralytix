"""
Controlador principal de la Estación de Trabajo de Cinética de Citometría de Flujo (.fcs).
Organizado en 5 pestañas persistentes con preanálisis interactivo, modelos no lineales y reporte auditable.
"""
import os
import re
import sys
import json
import math
import zlib
import struct
import traceback
import warnings
from pathlib import Path
from datetime import datetime

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.optimize import curve_fit
try:
    from scipy.integrate import trapezoid
except ImportError:
    from scipy.integrate import trapz as trapezoid

from ..config import (
    PEPTIDE_ALIASES, DOSE_RAMPS, GENERIC_RAMPS, CURVE_PALETTE, MAIN_COLORS,
    AVOGADRO, HAS_MPL, HAS_REPORTLAB
)
from ..io.fcs_parser import parse_fcs_native, get_fcs_channels, resolve_perm_channel
from ..io.metadata import guess_metadata_from_filename
from ..io.session import save_session_dict, load_session_dict
from ..core.statistics import (
    wilson_binomial_ci, cohen_h_proportions, contingency_2x2_p_and_or,
    cochran_armitage_trend_test, interpolate_crossing_time,
    classify_affectation_status, adjust_pvalues
)
from ..core.models import (
    hill_model, michaelis_menten_model, gompertz_dose_model,
    weibull_dose_model, biexponential_time_model
)
from ..core.kinetics import (
    extract_tube_kinetics, apply_gate_and_cutoffs_to_cache,
    extract_slope_based_parameters, compute_categorical_kinetics_metrics
)
from ..core.synthetic_demo import generate_demo_fcs_dataset
from ..export.plots_600dpi import (
    render_highres_plot_png, render_highres_histogram_png,
    render_highres_categorical_milestones_png,
    render_highres_categorical_heatmap_png,
    render_highres_target_vs_control_overlap_png
)
from ..export.pdf_report import generate_native_pdf_report
from .widgets import (
    NativeHistogramCutoffCanvas, NativeInteractiveMultiPlot,
    NativeSaturationPreview, NativeCategoricalProportionsCanvas
)
from .dialogs import (
    open_mathematical_models_help_dialog, open_user_guide_dialog,
    open_sample_cutoff_editor_dialog
)

DEFAULT_OUTPUT_DIR = Path.cwd() / "Resultados_Cinetica_Calceina"

class FlowCytometryKineticWorkstation:
    def __init__(self, root):
        self.root = root
        self.root.title("Estación de Trabajo de Cinética de Permeabilización Celular (Calceína-AM | .fcs)")

        sw = self.root.winfo_screenwidth()
        sh = self.root.winfo_screenheight()
        w = min(1380, int(sw * 0.95))
        h = min(860, int(sh * 0.91))
        x = max(5, (sw - w) // 2)
        y = max(5, (sh - h) // 2)
        self.root.geometry(f"{w}x{h}+{x}+{y}")
        self.root.minsize(1024, 680)

        self.file_records = []
        self.tube_event_cache = {}  # {filepath: dict con arrays de eventos}
        self.sample_filter_overrides = {}  # {filepath: dict con umbrales aplicados (Global o Discriminado)}
        self.df_full = pd.DataFrame()
        self.df_kinetics = pd.DataFrame()
        self.df_preregistration = pd.DataFrame()
        self.df_v0_summary = pd.DataFrame()
        self.df_v0_600_summary = pd.DataFrame()
        self.df_v0_active_summary = pd.DataFrame()
        self.df_categorical_summary = pd.DataFrame()
        self.df_cat_time_windows = pd.DataFrame()
        self.df_cat_hypothesis_tests = pd.DataFrame()
        self.df_morphology_summary = pd.DataFrame()
        self.df_table4 = pd.DataFrame()
        self.df_table5 = pd.DataFrame()
        self.df_saturation = pd.DataFrame()
        self.df_time_bin_overlap = pd.DataFrame()

        self.curve_thresholds = {}
        self._build_main_layout()

        self.root.deiconify()
        self.root.lift()
        self.root.focus_force()

    def _build_main_layout(self):
        # Barra de Menú Principal Nativa (Sesión/Configuración, Modelos Matemáticos y Guía de Uso)
        self.menubar = tk.Menu(self.root)

        menu_session = tk.Menu(self.menubar, tearoff=0)
        menu_session.add_command(
            label="💾 Guardar / Exportar Sesión de Trabajo (.json)...    [Ctrl+S]",
            command=self.save_session_to_file
        )
        menu_session.add_command(
            label="📂 Cargar / Importar Sesión de Trabajo (.json)...     [Ctrl+O]",
            command=self.load_session_from_file
        )
        menu_session.add_separator()
        menu_session.add_command(
            label="⚡ Guardado Rápido en Workspace (last_session_config.json)",
            command=self.quick_save_session
        )
        menu_session.add_command(
            label="🔄 Restaurar Última Sesión Guardada (last_session_config.json)",
            command=self.quick_load_session
        )
        menu_session.add_separator()
        menu_session.add_command(
            label="📋 Exportar Tabla de Cut-Offs por Muestra (.csv)...",
            command=self.export_sample_cutoffs_csv
        )
        menu_session.add_separator()
        menu_session.add_command(label="🚪 Salir", command=self.root.destroy)
        self.menubar.add_cascade(label="💾 Sesión y Configuración", menu=menu_session)

        menu_models = tk.Menu(self.menubar, tearoff=0)
        menu_models.add_command(
            label="📐 1. Modelos Farmacodinámicos y Cinéticos (Hill, Michaelis-Menten, Gompertz, Weibull, Bi-Exp)",
            command=lambda: self.open_mathematical_models_help_dialog(initial_tab=0)
        )
        menu_models.add_command(
            label="🔬 2. Gating (V1L/V1R), Cut-Off de Histograma (V2L/V2R) y Extracción Cinética Continua",
            command=lambda: self.open_mathematical_models_help_dialog(initial_tab=1)
        )
        menu_models.add_command(
            label="⏱️ 3. Hitos de Tiempo (T50 Efecto vs LT50, ΔT10-90), Entropía de Shannon, Hazard y Bimodalidad",
            command=lambda: self.open_mathematical_models_help_dialog(initial_tab=2)
        )
        menu_models.add_command(
            label="📊 4. Pruebas Estadísticas (Fisher, χ², h de Cohen, Wilson IC95%, Cochran-Armitage, Welch + FDR)",
            command=lambda: self.open_mathematical_models_help_dialog(initial_tab=3)
        )
        self.menubar.add_cascade(label="❓ Ayuda: Modelos Matemáticos y Cálculos", menu=menu_models)

        menu_guide = tk.Menu(self.menubar, tearoff=0)
        menu_guide.add_command(
            label="📖 Guía Paso a Paso de Uso del Software    [F1]",
            command=self.open_user_guide_dialog
        )
        menu_guide.add_separator()
        menu_guide.add_command(
            label="📄 Exportar Guía de Uso y Manual Matemático Completo (.md)...",
            command=self.export_markdown_guide_file
        )
        self.menubar.add_cascade(label="📖 Guía de Uso del Software", menu=menu_guide)

        self.root.config(menu=self.menubar)
        self.root.bind("<Control-s>", lambda e: self.save_session_to_file())
        self.root.bind("<Control-o>", lambda e: self.load_session_from_file())
        self.root.bind("<F1>", lambda e: self.open_user_guide_dialog())

        top_bar = tk.Frame(self.root, bg="#002B49", pady=5, padx=10)
        top_bar.pack(fill="x")
        tk.Label(
            top_bar,
            text="🧬 WORKSTATION DE CINÉTICA DE PERMEABILIZACIÓN (.fcs) — Preanálisis, Cut-Offs por Muestra, Modelos y Auditoría",
            bg="#002B49", fg="white", font=("DejaVu Sans", 9, "bold")
        ).pack(side="left")

        btn_bar = tk.Frame(top_bar, bg="#002B49")
        btn_bar.pack(side="right", padx=4)
        tk.Button(
            btn_bar, text="💾 Guardar Sesión", bg="#117864", fg="white",
            font=("DejaVu Sans", 8, "bold"), padx=6, pady=1,
            command=self.save_session_to_file
        ).pack(side="left", padx=2)
        tk.Button(
            btn_bar, text="📂 Cargar Sesión", bg="#1F618D", fg="white",
            font=("DejaVu Sans", 8, "bold"), padx=6, pady=1,
            command=self.load_session_from_file
        ).pack(side="left", padx=2)
        tk.Button(
            btn_bar, text="📐 Modelos y Cálculos", bg="#6C3483", fg="white",
            font=("DejaVu Sans", 8, "bold"), padx=6, pady=1,
            command=lambda: self.open_mathematical_models_help_dialog(initial_tab=0)
        ).pack(side="left", padx=2)
        tk.Button(
            btn_bar, text="📖 Guía de Uso", bg="#B9770E", fg="white",
            font=("DejaVu Sans", 8, "bold"), padx=6, pady=1,
            command=self.open_user_guide_dialog
        ).pack(side="left", padx=2)

        self.lbl_global_status = tk.Label(
            top_bar, text="Estado: Listo (Motor Nativo Activo)",
            bg="#002B49", fg="#A3E4D7", font=("DejaVu Sans", 8, "bold")
        )
        self.lbl_global_status.pack(side="right", padx=8)

        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True, padx=8, pady=6)

        self.tab_ingest = ttk.Frame(self.notebook)
        self.tab_preanalysis = ttk.Frame(self.notebook)
        self.tab_catpanel = ttk.Frame(self.notebook)
        self.tab_prereport = ttk.Frame(self.notebook)
        self.tab_finalreport = ttk.Frame(self.notebook)

        self.notebook.add(self.tab_ingest, text=" 📁 1. Carga de Datos y Guía ")
        self.notebook.add(self.tab_preanalysis, text=" 🔍 2. Pre-Análisis Dinámico (FITC / FSC / SSC) ")
        self.notebook.add(self.tab_catpanel, text=" 🧬 3. Panel Categórico (Afectados vs No Afectados y Ventanas) ")
        self.notebook.add(self.tab_prereport, text=" 📊 4. Pre-Reporte, Modelos Cinéticos y Estadística ")
        self.notebook.add(self.tab_finalreport, text=" 📄 5. Reporte Final y Exportación (600 DPI / PDF) ")

        self._build_tab1_ingest()
        self._build_tab2_preanalysis()
        self._build_tab_categorical_panel()
        self._build_tab3_prereport()
        self._build_tab4_finalreport()

    # --------------------------------------------------------------------------
    # PESTAÑA 1: CARGA DE DATOS, GUÍA DE NOMENCLATURA Y EDICIÓN DE METADATOS
    # --------------------------------------------------------------------------
    def _build_tab1_ingest(self):
        guide_frame = ttk.LabelFrame(
            self.tab_ingest,
            text=" 💡 Guía de Nomenclatura de Archivos .fcs para Autodetección Inmediata ",
            padding=8
        )
        guide_frame.pack(fill="x", padx=10, pady=5)

        guide_txt = (
            "Formato recomendado:   [Compuesto]_[LineaCelular]_[Concentracion_uM]_[CC_opcional]_[Replica].fcs\n"
            "  • Control Basal Previo (-60 a 0 s): Incluir 'CC' o 'BASE'  ➔  Ej:  820_THP-1_CC_4.35_1.fcs   o   Sat1_PBMC_CC_5.70_R1.fcs\n"
            "  • Muestra Experimental (0 a 600 s): Sin 'CC' en el nombre  ➔  Ej:  820_THP-1_4.35_1.fcs      o   Sat1_THP-1_2.80_R2.fcs\n"
            "Nota: Si tus archivos tienen otros nombres, puedes editar cualquier celda con Doble Clic en la tabla inferior."
        )
        ttk.Label(guide_frame, text=guide_txt, foreground="#1B4F72", font=("DejaVu Sans", 9)).pack(anchor="w")

        cfg_frame = ttk.LabelFrame(
            self.tab_ingest,
            text=" Parámetros Celulares Globales y Selección del Canal de Fluorescencia ",
            padding=8
        )
        cfg_frame.pack(fill="x", padx=10, pady=4)

        ttk.Label(cfg_frame, text="Célula Diana (Target):").grid(row=0, column=0, sticky="w", padx=5, pady=3)
        self.ent_target = ttk.Entry(cfg_frame, width=16)
        self.ent_target.insert(0, "THP-1")
        self.ent_target.grid(row=0, column=1, sticky="w", padx=5, pady=3)

        ttk.Label(cfg_frame, text="Célula Control Sano:").grid(row=0, column=2, sticky="w", padx=5, pady=3)
        self.ent_control = ttk.Entry(cfg_frame, width=16)
        self.ent_control.insert(0, "PBMCs")
        self.ent_control.grid(row=0, column=3, sticky="w", padx=5, pady=3)

        ttk.Label(cfg_frame, text="Densidad Celular (cél/mL):").grid(row=0, column=4, sticky="w", padx=5, pady=3)
        self.ent_density = ttk.Entry(cfg_frame, width=14)
        self.ent_density.insert(0, "1.0e6")
        self.ent_density.grid(row=0, column=5, sticky="w", padx=5, pady=3)

        ttk.Label(cfg_frame, text="Canal Permeabilidad:", font=("DejaVu Sans", 9, "bold")).grid(row=0, column=6, sticky="w", padx=5, pady=3)
        self.var_perm_channel = tk.StringVar(value="FITC (Por defecto)")
        self.cmb_tab1_channel = ttk.Combobox(
            cfg_frame, textvariable=self.var_perm_channel, state="readonly", width=22,
            values=["FITC (Por defecto)", "FITC-A", "FL1-A", "FITC-H", "FL1-H", "BL1-A"]
        )
        self.cmb_tab1_channel.grid(row=0, column=7, sticky="w", padx=5, pady=3)

        btn_frame = ttk.LabelFrame(self.tab_ingest, text=" Gestión de Archivos Binarios (.fcs) ", padding=8)
        btn_frame.pack(fill="x", padx=10, pady=4)

        ttk.Button(btn_frame, text="➕ Subir Archivos .fcs", command=self.add_fcs_files).pack(side="left", padx=4)
        ttk.Button(btn_frame, text="📁 Subir Carpeta Completa", command=self.add_fcs_folder).pack(side="left", padx=4)
        ttk.Button(btn_frame, text="🧪 Generar y Cargar Set Demo (.fcs)", command=self.load_demo_dataset).pack(side="left", padx=4)
        ttk.Button(btn_frame, text="✏️ Editar en Lote Seleccionados", command=self.batch_edit_selected).pack(side="left", padx=4)
        ttk.Button(btn_frame, text="🗑️ Quitar Seleccionados", command=self.remove_selected).pack(side="left", padx=4)
        ttk.Button(btn_frame, text="🧹 Limpiar Tabla", command=self.clear_all).pack(side="left", padx=4)

        tk.Button(
            btn_frame, text="🔍 EJECUTAR PRE-ANÁLISIS Y VER GRÁFICAS ➔",
            bg="#004488", fg="white", font=("DejaVu Sans", 9, "bold"),
            padx=12, pady=4, command=self.run_preanalysis_ingestion
        ).pack(side="right", padx=6)

        tbl_frame = ttk.LabelFrame(
            self.tab_ingest,
            text=" Tabla de Muestras (Doble clic en cualquier celda para editar Compuesto, Célula, Rol CC, Concentración µM, Réplica o Exclusión t_start) ",
            padding=8
        )
        tbl_frame.pack(fill="both", expand=True, padx=10, pady=5)

        columns = ("Filename", "Compound", "Cell_Line", "Sample_Role", "Concentration_uM", "Replicate", "Exclusion_s")
        self.tree = ttk.Treeview(tbl_frame, columns=columns, show="headings", selectmode="extended")
        headers = {
            "Filename": ("Archivo .fcs", 230),
            "Compound": ("Compuesto / Péptido", 145),
            "Cell_Line": ("Línea Celular", 115),
            "Sample_Role": ("Rol (Control CC vs Experimental)", 210),
            "Concentration_uM": ("Concentración (µM)", 130),
            "Replicate": ("Réplica (R1..Rn)", 105),
            "Exclusion_s": ("Exclusión Inicial t_start (s)", 150)
        }
        for col, (text, width) in headers.items():
            self.tree.heading(col, text=text)
            self.tree.column(col, width=width, anchor="center" if col != "Filename" else "w")

        vsb = ttk.Scrollbar(tbl_frame, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(tbl_frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        tbl_frame.grid_rowconfigure(0, weight=1)
        tbl_frame.grid_columnconfigure(0, weight=1)

        self.tree.bind("<Double-1>", self.on_double_click_cell)

    def _refresh_available_channels_from_files(self):
        if not self.file_records:
            return
        chans = get_fcs_channels(self.file_records[0]['Filepath'])
        if chans:
            opts = ["FITC (Por defecto)"] + [c for c in chans if c.upper() != "TIME"]
            self.cmb_tab1_channel['values'] = opts
            if hasattr(self, 'cmb_perm_channel'):
                self.cmb_perm_channel['values'] = opts

    def add_fcs_files(self):
        files = filedialog.askopenfilenames(
            title="Seleccionar archivos .fcs",
            filetypes=[("Flow Cytometry Standard (.fcs)", "*.fcs *.FCS"), ("Todos los archivos", "*.*")],
            parent=self.root
        )
        if files:
            self._prompt_batch_assignment_and_insert(files)

    def add_fcs_folder(self):
        folder = filedialog.askdirectory(title="Seleccionar carpeta con archivos .fcs", parent=self.root)
        if not folder:
            return
        fcs_files = []
        for root_dir, _, files in os.walk(folder):
            for f in files:
                if f.lower().endswith('.fcs'):
                    fcs_files.append(os.path.join(root_dir, f))
        if not fcs_files:
            messagebox.showwarning("Sin archivos", "No se encontraron archivos .fcs en la carpeta seleccionada.", parent=self.root)
            return
        self._prompt_batch_assignment_and_insert(sorted(fcs_files))

    def load_demo_dataset(self):
        out_folder = Path.cwd() / "Resultados_Cinetica_Calceina"
        created = generate_demo_fcs_dataset(out_folder)
        existing_paths = {r['Filepath'] for r in self.file_records}
        t_cell = self.ent_target.get().strip() or "THP-1"
        c_cell = self.ent_control.get().strip() or "PBMCs"
        for fp in created:
            if str(fp) in existing_paths:
                continue
            meta = guess_metadata_from_filename(fp, t_cell, c_cell)
            self.file_records.append(meta)
            self.tree.insert("", "end", iid=str(fp), values=(
                meta['Filename'], meta['Compound'], meta['Cell_Line'],
                meta['Sample_Role'], f"{meta['Concentration_uM']:.2f}",
                meta['Replicate'], f"{meta['Exclusion_s']:.1f}"
            ))
        self._refresh_available_channels_from_files()
        self.lbl_global_status.config(text=f"Estado: {len(self.file_records)} archivos Demo cargados")
        self.run_preanalysis_ingestion()

    def _prompt_batch_assignment_and_insert(self, filepaths):
        dlg = tk.Toplevel(self.root)
        dlg.title(f"Configurar Lote ({len(filepaths)} archivos .fcs)")
        dlg.geometry("460x270")
        dlg.transient(self.root)
        dlg.grab_set()

        ttk.Label(dlg, text=f"Se detectaron {len(filepaths)} archivos .fcs", font=("DejaVu Sans", 9, "bold")).pack(pady=6)
        frm = ttk.Frame(dlg, padding=10)
        frm.pack(fill="both", expand=True)

        raw_folder = os.path.basename(os.path.dirname(filepaths[0]))
        default_comp = PEPTIDE_ALIASES.get(raw_folder.lower(), raw_folder)
        ttk.Label(frm, text="Compuesto (dejar '(Auto)' para deducir):").grid(row=0, column=0, sticky="w", pady=4)
        ent_c = ttk.Entry(frm, width=22)
        ent_c.insert(0, default_comp if default_comp else "(Auto)")
        ent_c.grid(row=0, column=1, pady=4)

        ttk.Label(frm, text="Línea Celular:").grid(row=1, column=0, sticky="w", pady=4)
        cmb_cell = ttk.Combobox(frm, values=["(Auto-detectar)", self.ent_target.get(), self.ent_control.get()], width=20)
        cmb_cell.current(0)
        cmb_cell.grid(row=1, column=1, pady=4)

        ttk.Label(frm, text="Exclusión Inicial t_start (s) ('Auto' o seg):").grid(row=2, column=0, sticky="w", pady=4)
        ent_excl = ttk.Entry(frm, width=22)
        ent_excl.insert(0, "Auto")
        ent_excl.grid(row=2, column=1, pady=4)

        def confirm():
            comp = ent_c.get().strip()
            cell_choice = cmb_cell.get().strip()
            excl_raw = ent_excl.get().strip()

            existing = {r['Filepath'] for r in self.file_records}
            for fp in filepaths:
                if str(fp) in existing:
                    continue
                meta = guess_metadata_from_filename(fp, self.ent_target.get().strip(), self.ent_control.get().strip())
                if comp and comp != "(Auto)":
                    meta['Compound'] = PEPTIDE_ALIASES.get(comp.lower(), comp)
                if cell_choice != "(Auto-detectar)":
                    meta['Cell_Line'] = cell_choice
                if excl_raw.lower() != "auto":
                    try:
                        meta['Exclusion_s'] = float(excl_raw)
                    except ValueError:
                        pass

                self.file_records.append(meta)
                self.tree.insert("", "end", iid=str(fp), values=(
                    meta['Filename'], meta['Compound'], meta['Cell_Line'],
                    meta['Sample_Role'], f"{meta['Concentration_uM']:.2f}",
                    meta['Replicate'], f"{meta['Exclusion_s']:.1f}"
                ))
            self._refresh_available_channels_from_files()
            self.lbl_global_status.config(text=f"Estado: {len(self.file_records)} archivos en lista")
            dlg.destroy()

        ttk.Button(dlg, text="Aceptar y Cargar a la Tabla", command=confirm).pack(pady=8)

    def on_double_click_cell(self, event):
        region = self.tree.identify("region", event.x, event.y)
        if region != "cell":
            return
        col = self.tree.identify_column(event.x)
        row_id = self.tree.identify_row(event.y)
        if not row_id or col == "#1":
            return

        col_idx = int(col.replace("#", "")) - 1
        col_name = self.tree["columns"][col_idx]
        x, y, width, height = self.tree.bbox(row_id, col)
        current_val = self.tree.item(row_id, "values")[col_idx]

        if col_name == "Sample_Role":
            editor = ttk.Combobox(self.tree, values=["CC (Baseline -60 to 0s)", "Experimental (0 to 600s)"], state="readonly")
            editor.set(current_val)
        elif col_name == "Cell_Line":
            editor = ttk.Combobox(self.tree, values=[self.ent_target.get().strip(), self.ent_control.get().strip()])
            editor.set(current_val)
        else:
            editor = ttk.Entry(self.tree)
            editor.insert(0, current_val)

        editor.place(x=x, y=y, width=width, height=height)
        editor.focus()

        def save_edit(ev=None):
            new_val = editor.get().strip()
            vals = list(self.tree.item(row_id, "values"))
            if col_name in ("Concentration_uM", "Exclusion_s"):
                try:
                    num = float(new_val)
                    new_val = f"{num:.2f}" if col_name == "Concentration_uM" else f"{num:.1f}"
                except ValueError:
                    editor.destroy()
                    return
            vals[col_idx] = new_val
            self.tree.item(row_id, values=vals)
            for rec in self.file_records:
                if rec['Filepath'] == row_id:
                    if col_name in ("Concentration_uM", "Exclusion_s"):
                        rec[col_name] = float(new_val)
                    else:
                        rec[col_name] = new_val
                    break
            editor.destroy()

        editor.bind("<Return>", save_edit)
        editor.bind("<FocusOut>", save_edit)

    def batch_edit_selected(self):
        selected = self.tree.selection()
        if not selected:
            messagebox.showinfo("Selección", "Selecciona primero una o más filas en la tabla.", parent=self.root)
            return

        dlg = tk.Toplevel(self.root)
        dlg.title(f"Edición en Lote ({len(selected)} archivos)")
        dlg.geometry("410x290")
        dlg.transient(self.root)
        dlg.grab_set()

        frm = ttk.Frame(dlg, padding=10)
        frm.pack(fill="both", expand=True)

        fields = [
            ("Compuesto:", "Compound"),
            ("Línea Celular:", "Cell_Line"),
            ("Rol:", "Sample_Role"),
            ("Dosis (µM):", "Concentration_uM"),
            ("Réplica:", "Replicate"),
            ("Exclusión t_start (s):", "Exclusion_s")
        ]
        entries = {}
        for idx, (lbl, key) in enumerate(fields):
            ttk.Label(frm, text=lbl).grid(row=idx, column=0, sticky="w", pady=3)
            if key == "Sample_Role":
                w = ttk.Combobox(frm, values=["", "CC (Baseline -60 to 0s)", "Experimental (0 to 600s)"], width=22)
            elif key == "Cell_Line":
                w = ttk.Combobox(frm, values=["", self.ent_target.get().strip(), self.ent_control.get().strip()], width=22)
            else:
                w = ttk.Entry(frm, width=24)
            w.grid(row=idx, column=1, pady=3)
            entries[key] = w

        def apply_batch():
            for row_id in selected:
                vals = list(self.tree.item(row_id, "values"))
                rec = next((r for r in self.file_records if r['Filepath'] == row_id), None)
                if not rec:
                    continue
                for col_idx, key in enumerate(["Filename", "Compound", "Cell_Line", "Sample_Role", "Concentration_uM", "Replicate", "Exclusion_s"]):
                    if key in entries and entries[key].get().strip() != "":
                        raw_v = entries[key].get().strip()
                        if key in ("Concentration_uM", "Exclusion_s"):
                            rec[key] = float(raw_v)
                            vals[col_idx] = f"{rec[key]:.2f}" if key == "Concentration_uM" else f"{rec[key]:.1f}"
                        else:
                            rec[key] = raw_v
                            vals[col_idx] = raw_v
                self.tree.item(row_id, values=vals)
            dlg.destroy()

        ttk.Button(dlg, text="Aplicar Cambios", command=apply_batch).pack(pady=8)

    def remove_selected(self):
        for row_id in self.tree.selection():
            self.tree.delete(row_id)
            self.file_records = [r for r in self.file_records if r['Filepath'] != row_id]

    def clear_all(self):
        for row_id in self.tree.get_children():
            self.tree.delete(row_id)
        self.file_records.clear()

    # --------------------------------------------------------------------------
    # PESTAÑA 2: PANEL DE PRE-ANÁLISIS DINÁMICO (HISTOGRAMA CUT-OFF + FITC/FSC/SSC)
    # --------------------------------------------------------------------------
    def _build_tab2_preanalysis(self):
        paned = ttk.PanedWindow(self.tab_preanalysis, orient="horizontal")
        paned.pack(fill="both", expand=True, padx=6, pady=6)

        left_panel = ttk.Frame(paned, width=460)
        paned.add(left_panel, weight=0)

        f_gate = ttk.LabelFrame(left_panel, text=" 1. Canal, Prefiltros (FSC/SSC/V1L) y Cut-Off V2L/V2R ", padding=5)
        f_gate.pack(fill="x", padx=4, pady=2)

        ttk.Label(f_gate, text="Canal Permeabilidad:", font=("DejaVu Sans", 8, "bold")).grid(row=0, column=0, columnspan=2, sticky="w", pady=1)
        self.cmb_perm_channel = ttk.Combobox(
            f_gate, textvariable=self.var_perm_channel, state="readonly", width=22,
            values=["FITC (Por defecto)", "FITC-A", "FL1-A", "FITC-H", "FL1-H", "BL1-A"]
        )
        self.cmb_perm_channel.grid(row=0, column=2, columnspan=2, sticky="ew", padx=3, pady=1)
        self.cmb_perm_channel.bind("<<ComboboxSelected>>", lambda e: self.run_preanalysis_ingestion())

        ttk.Label(f_gate, text="FSC-A Mín (>):").grid(row=1, column=0, sticky="w", pady=1)
        self.ent_fsc_min = ttk.Entry(f_gate, width=9)
        self.ent_fsc_min.insert(0, "50000")
        self.ent_fsc_min.grid(row=1, column=1, padx=2, pady=1)

        ttk.Label(f_gate, text="FSC-A Máx:").grid(row=1, column=2, sticky="w", pady=1)
        self.ent_fsc_max = ttk.Entry(f_gate, width=9)
        self.ent_fsc_max.insert(0, "Sin límite")
        self.ent_fsc_max.grid(row=1, column=3, padx=2, pady=1)

        ttk.Label(f_gate, text="SSC-A Mín (>):").grid(row=2, column=0, sticky="w", pady=1)
        self.ent_ssc_min = ttk.Entry(f_gate, width=9)
        self.ent_ssc_min.insert(0, "10000")
        self.ent_ssc_min.grid(row=2, column=1, padx=2, pady=1)

        ttk.Label(f_gate, text="SSC-A Máx:").grid(row=2, column=2, sticky="w", pady=1)
        self.ent_ssc_max = ttk.Entry(f_gate, width=9)
        self.ent_ssc_max.insert(0, "Sin límite")
        self.ent_ssc_max.grid(row=2, column=3, padx=2, pady=1)

        # Cut-Off 1: Exclusión de fluorescencia V1L (omite células sin cargar / detritos)
        ttk.Label(f_gate, text="Exclusión V1L (<log10):", foreground="#4A235A", font=("DejaVu Sans", 8, "bold")).grid(row=3, column=0, sticky="w", pady=1)
        self.var_log_excl_v1l = tk.DoubleVar(value=3.00)
        self.ent_log_excl_v1l = ttk.Entry(f_gate, textvariable=self.var_log_excl_v1l, width=9)
        self.ent_log_excl_v1l.grid(row=3, column=1, padx=2, pady=1)
        self.ent_log_excl_v1l.bind("<Return>", lambda e: self.apply_prefilter_by_current_scope())

        # Cut-Off 2: Discriminación V2L (Afectados) vs V2R (No Afectados)
        ttk.Label(f_gate, text="Cut-Off V2L/V2R (log10):", foreground="#C0392B", font=("DejaVu Sans", 8, "bold")).grid(row=3, column=2, sticky="w", pady=1)
        self.var_log_cutoff = tk.DoubleVar(value=4.85)
        self.ent_log_cutoff = ttk.Entry(f_gate, textvariable=self.var_log_cutoff, width=9)
        self.ent_log_cutoff.grid(row=3, column=3, padx=2, pady=1)
        self.ent_log_cutoff.bind("<Return>", lambda e: self.apply_prefilter_by_current_scope())

        ttk.Label(f_gate, text="Bin Temporal (s):").grid(row=4, column=0, sticky="w", pady=1)
        self.ent_bin_size = ttk.Entry(f_gate, width=9)
        self.ent_bin_size.insert(0, "15")
        self.ent_bin_size.grid(row=4, column=1, padx=2, pady=1)

        ttk.Label(f_gate, text="Eje Y Hist.:").grid(row=4, column=2, sticky="w", pady=1)
        self.var_hist_ymode = tk.StringVar(value="Porcentaje (%)")
        cmb_hym = ttk.Combobox(f_gate, textvariable=self.var_hist_ymode, state="readonly", width=11, values=["Porcentaje (%)", "N° Eventos"])
        cmb_hym.grid(row=4, column=3, sticky="ew", padx=2, pady=1)
        cmb_hym.bind("<<ComboboxSelected>>", lambda e: self.update_preanalysis_plots())

        ttk.Label(f_gate, text="Alcance Filtro/Cut-Off:", font=("DejaVu Sans", 8, "bold")).grid(row=5, column=0, sticky="w", pady=1)
        self.var_cutoff_scope = tk.StringVar(value="Global (Todas las muestras)")
        cmb_scope = ttk.Combobox(
            f_gate, textvariable=self.var_cutoff_scope, state="readonly", width=26,
            values=["Global (Todas las muestras)", "Discriminado (Muestras seleccionadas)"]
        )
        cmb_scope.grid(row=5, column=1, columnspan=3, sticky="ew", padx=2, pady=1)

        btn_apply_frm = ttk.Frame(f_gate)
        btn_apply_frm.grid(row=6, column=0, columnspan=4, sticky="ew", pady=2)
        tk.Button(
            btn_apply_frm, text="🌐 Aplicar a TODAS",
            bg="#1B4F72", fg="white", font=("DejaVu Sans", 8, "bold"),
            command=lambda: self.apply_prefilter_and_cutoffs(scope="global")
        ).pack(side="left", expand=True, fill="x", padx=1)
        tk.Button(
            btn_apply_frm, text="🎯 Aplicar a SELECCIONADAS",
            bg="#633974", fg="white", font=("DejaVu Sans", 8, "bold"),
            command=lambda: self.apply_prefilter_and_cutoffs(scope="selected")
        ).pack(side="left", expand=True, fill="x", padx=1)
        ttk.Button(
            btn_apply_frm, text="Auto CC", width=8,
            command=self.auto_calibrate_cutoff_from_cc
        ).pack(side="left", padx=1)

        ttk.Label(f_gate, text="Modo Gráfica 1:").grid(row=7, column=0, sticky="w", pady=1)
        self.var_perm_display_mode = tk.StringVar(value="Dual (Continua + Afectados V2L %)")
        cmb_pdisp = ttk.Combobox(
            f_gate, textvariable=self.var_perm_display_mode, state="readonly", width=26,
            values=[
                "Continua (Mediana FITC %)",
                "Categórica (Afectados V2L % por Cut-Off)",
                "Dual (Continua + Afectados V2L %)"
            ]
        )
        cmb_pdisp.grid(row=7, column=1, columnspan=3, sticky="ew", padx=2, pady=1)
        cmb_pdisp.bind("<<ComboboxSelected>>", lambda e: self.update_preanalysis_plots())

        self.var_scatter_norm = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            f_gate, text="Mostrar FSC/SSC Normalizados (% vs Basal)",
            variable=self.var_scatter_norm, command=self.update_preanalysis_plots
        ).grid(row=8, column=0, columnspan=4, sticky="w", pady=1)

        self.lbl_gate_stats = ttk.Label(
            f_gate, text="Canal: FITC | V1R Retenidos: -- | V1L Omitidos: --",
            foreground="#117864", font=("DejaVu Sans", 8, "bold")
        )
        self.lbl_gate_stats.grid(row=9, column=0, columnspan=4, sticky="w")

        f_sel = ttk.LabelFrame(left_panel, text=" 2. Selección Dinámica de Muestras a Graficar ", padding=4)
        f_sel.pack(fill="x", padx=4, pady=1)

        sel_grid = ttk.Frame(f_sel)
        sel_grid.pack(fill="x")
        sel_grid.columnconfigure(0, weight=1)
        sel_grid.columnconfigure(1, weight=1)
        sel_grid.columnconfigure(2, weight=1)

        ttk.Label(sel_grid, text="Compuestos:", font=("DejaVu Sans", 8, "bold")).grid(row=0, column=0, sticky="w", padx=2)
        ttk.Label(sel_grid, text="Líneas:", font=("DejaVu Sans", 8, "bold")).grid(row=0, column=1, sticky="w", padx=2)
        ttk.Label(sel_grid, text="Dosis (µM):", font=("DejaVu Sans", 8, "bold")).grid(row=0, column=2, sticky="w", padx=2)

        self.lst_compounds = tk.Listbox(sel_grid, selectmode="extended", height=3, width=13, exportselection=False)
        self.lst_compounds.grid(row=1, column=0, sticky="ew", padx=2, pady=1)
        self.lst_compounds.bind("<<ListboxSelect>>", lambda e: self.update_preanalysis_plots())

        self.lst_cells = tk.Listbox(sel_grid, selectmode="extended", height=3, width=11, exportselection=False)
        self.lst_cells.grid(row=1, column=1, sticky="ew", padx=2, pady=1)
        self.lst_cells.bind("<<ListboxSelect>>", lambda e: self.update_preanalysis_plots())

        self.lst_doses = tk.Listbox(sel_grid, selectmode="extended", height=3, width=10, exportselection=False)
        self.lst_doses.grid(row=1, column=2, sticky="ew", padx=2, pady=1)
        self.lst_doses.bind("<<ListboxSelect>>", lambda e: self.update_preanalysis_plots())

        btn_sel_frame = ttk.Frame(f_sel)
        btn_sel_frame.pack(fill="x", pady=1)
        ttk.Button(btn_sel_frame, text="Seleccionar Todas", command=self.select_all_filters).pack(side="left", expand=True, fill="x", padx=2)
        ttk.Button(btn_sel_frame, text="Refrescar Gráficas", command=self.update_preanalysis_plots).pack(side="left", expand=True, fill="x", padx=2)

        f_thresh = ttk.LabelFrame(left_panel, text=" 3. Ajuste de T_lag y T_max por Curva (Arrastrable) ", padding=4)
        f_thresh.pack(fill="x", padx=4, pady=1)

        curv_row = ttk.Frame(f_thresh)
        curv_row.pack(fill="x", pady=1)
        ttk.Label(curv_row, text="Curva Activa:", font=("DejaVu Sans", 8, "bold")).pack(side="left")
        self.cmb_active_curve = ttk.Combobox(curv_row, state="readonly", width=29)
        self.cmb_active_curve.pack(side="left", fill="x", expand=True, padx=3)
        self.cmb_active_curve.bind("<<ComboboxSelected>>", self.on_active_curve_selected)

        self.var_thresh_mode = tk.StringVar(value="auto")
        mode_frm = ttk.Frame(f_thresh)
        mode_frm.pack(fill="x", pady=1)
        ttk.Radiobutton(
            mode_frm, text="Automático", value="auto",
            variable=self.var_thresh_mode, command=self.on_threshold_mode_toggle
        ).pack(side="left", padx=2)
        ttk.Radiobutton(
            mode_frm, text="Manual / Arrastrar", value="manual",
            variable=self.var_thresh_mode, command=self.on_threshold_mode_toggle
        ).pack(side="left", padx=2)
        tk.Button(
            mode_frm, text="⚡ Auto (Todas)", bg="#117864", fg="white",
            font=("DejaVu Sans", 7, "bold"), padx=4, pady=0,
            command=self.reset_all_curves_to_auto_tlag_tmax
        ).pack(side="right", padx=1)
        tk.Button(
            mode_frm, text="⚡ Auto (Activa)", bg="#1B4F72", fg="white",
            font=("DejaVu Sans", 7, "bold"), padx=4, pady=0,
            command=self.reset_active_curve_to_auto_tlag_tmax
        ).pack(side="right", padx=1)

        sliders_grid = ttk.Frame(f_thresh)
        sliders_grid.pack(fill="x")
        sliders_grid.columnconfigure(0, weight=1)
        sliders_grid.columnconfigure(1, weight=1)

        ttk.Label(sliders_grid, text="T_lag (s) [Azul :]:", font=("DejaVu Sans", 8)).grid(row=0, column=0, sticky="w")
        ttk.Label(sliders_grid, text="T_max (s) [Rojo -.]:", font=("DejaVu Sans", 8)).grid(row=0, column=1, sticky="w")

        self.var_tlag_slider = tk.DoubleVar(value=60.0)
        self.scl_tlag = tk.Scale(
            sliders_grid, from_=0, to=600, resolution=15, orient="horizontal",
            variable=self.var_tlag_slider, length=160
        )
        self.scl_tlag.grid(row=1, column=0, sticky="ew", padx=2)
        self.scl_tlag.bind("<B1-Motion>", self.on_slider_thresh_change)
        self.scl_tlag.bind("<ButtonRelease-1>", lambda e: self.on_slider_thresh_change(e, final=True))

        self.var_tmax_slider = tk.DoubleVar(value=300.0)
        self.scl_tmax = tk.Scale(
            sliders_grid, from_=15, to=600, resolution=15, orient="horizontal",
            variable=self.var_tmax_slider, length=160
        )
        self.scl_tmax.grid(row=1, column=1, sticky="ew", padx=2)
        self.scl_tmax.bind("<B1-Motion>", self.on_slider_thresh_change)
        self.scl_tmax.bind("<ButtonRelease-1>", lambda e: self.on_slider_thresh_change(e, final=True))

        self.lbl_active_metrics = tk.Label(
            f_thresh,
            text="T_lag: -- s | T_max: -- s\nV_active: -- %/s | V2L Máx: -- % | LT50: -- s",
            bg="#EBF5FB", fg="#002B49", font=("DejaVu Sans", 8, "bold"),
            relief="groove", pady=3
        )
        self.lbl_active_metrics.pack(fill="x", pady=2)

        # Caja 4: Tabla Interactiva de Cut-Offs y Ventanas por Muestra (Clic = Seleccionar, Doble Clic = Editar)
        f_cut_tbl = ttk.LabelFrame(
            left_panel,
            text=" 4. Tabla Interactiva de Cut-Offs por Muestra (Doble Clic en Celda = Editar) ",
            padding=4
        )
        f_cut_tbl.pack(fill="both", expand=True, padx=4, pady=2)

        tbl_container = ttk.Frame(f_cut_tbl)
        tbl_container.pack(fill="both", expand=True)

        cut_cols = (
            "Muestra", "Scope", "V1L_Log", "V2L_Log",
            "FSC_Min", "SSC_Min", "T_lag_s", "T_max_s", "Modo_T", "V2L_Max_Pct"
        )
        self.tree_sample_cutoffs = ttk.Treeview(
            tbl_container, columns=cut_cols, show="headings", height=6, selectmode="browse"
        )
        col_specs = [
            ("Muestra", "Muestra (Comp | Cél | µM)", 150, "w"),
            ("Scope", "Alcance", 74, "center"),
            ("V1L_Log", "V1L <log", 58, "center"),
            ("V2L_Log", "CutOff V2L", 68, "center"),
            ("FSC_Min", "FSC Mín", 62, "center"),
            ("SSC_Min", "SSC Mín", 62, "center"),
            ("T_lag_s", "T_lag(s)", 56, "center"),
            ("T_max_s", "T_max(s)", 56, "center"),
            ("Modo_T", "Modo T", 56, "center"),
            ("V2L_Max_Pct", "%V2L Máx", 64, "center")
        ]
        for cid, heading_txt, w_px, anc in col_specs:
            self.tree_sample_cutoffs.heading(cid, text=heading_txt)
            self.tree_sample_cutoffs.column(cid, width=w_px, anchor=anc, stretch=False)

        vsb_cut = ttk.Scrollbar(tbl_container, orient="vertical", command=self.tree_sample_cutoffs.yview)
        hsb_cut = ttk.Scrollbar(tbl_container, orient="horizontal", command=self.tree_sample_cutoffs.xview)
        self.tree_sample_cutoffs.configure(yscrollcommand=vsb_cut.set, xscrollcommand=hsb_cut.set)

        self.tree_sample_cutoffs.grid(row=0, column=0, sticky="nsew")
        vsb_cut.grid(row=0, column=1, sticky="ns")
        hsb_cut.grid(row=1, column=0, sticky="ew")
        tbl_container.rowconfigure(0, weight=1)
        tbl_container.columnconfigure(0, weight=1)

        self.tree_sample_cutoffs.bind("<<TreeviewSelect>>", self.on_select_sample_cutoff_table)
        self.tree_sample_cutoffs.bind("<Double-1>", self.on_double_click_sample_cutoff_table)

        cut_btn_row = ttk.Frame(f_cut_tbl)
        cut_btn_row.pack(fill="x", pady=2)
        tk.Button(
            cut_btn_row, text="✏️ Editar Muestra Sel.",
            bg="#1F618D", fg="white", font=("DejaVu Sans", 8, "bold"),
            command=self.open_selected_sample_cutoff_editor_modal
        ).pack(side="left", expand=True, fill="x", padx=1)
        tk.Button(
            cut_btn_row, text="📋 Exportar Cut-Offs (.csv)",
            bg="#117864", fg="white", font=("DejaVu Sans", 8, "bold"),
            command=self.export_sample_cutoffs_csv
        ).pack(side="left", expand=True, fill="x", padx=1)
        tk.Button(
            cut_btn_row, text="🧬 Panel Categórico ➔",
            bg="#004488", fg="white", font=("DejaVu Sans", 8, "bold"),
            command=lambda: self.notebook.select(self.tab_catpanel)
        ).pack(side="left", expand=True, fill="x", padx=1)

        right_panel = ttk.PanedWindow(paned, orient="vertical")
        paned.add(right_panel, weight=1)

        # Panel Superior Derecho: Histograma de Fluorescencia con Cut-Off Arrastrable (V2L vs V2R)
        self.hist_cutoff_plot = NativeHistogramCutoffCanvas(
            right_panel, on_cutoff_change_callback=self.apply_new_log_cutoff
        )
        right_panel.add(self.hist_cutoff_plot, weight=1)

        # Panel Inferior Derecho: 3 Gráficas Temporales Sincronizadas
        self.multi_plot = NativeInteractiveMultiPlot(
            right_panel, on_threshold_drag_callback=self.on_canvas_drag_threshold
        )
        right_panel.add(self.multi_plot, weight=3)

    # --------------------------------------------------------------------------
    # PESTAÑA 3: PANEL DE ANÁLISIS DE PROPORCIONES CATEGÓRICAS Y VENTANAS DE TIEMPO
    # --------------------------------------------------------------------------
    def _build_tab_categorical_panel(self):
        ctrl_cat = ttk.LabelFrame(
            self.tab_catpanel,
            text=" Configuración del Panel Categórico: Proporciones Afectados (V2L) vs. No Afectados (V2R), Ventanas de Tiempo, T50, T_efecto_máx y Pruebas de Afectación ",
            padding=7
        )
        ctrl_cat.pack(fill="x", padx=10, pady=4)

        ttk.Label(ctrl_cat, text="Compuesto:").grid(row=0, column=0, sticky="w", padx=4, pady=2)
        self.var_cat_pep = tk.StringVar(value="Todos")
        self.cmb_cat_pep = ttk.Combobox(ctrl_cat, textvariable=self.var_cat_pep, state="readonly", width=18, values=["Todos"])
        self.cmb_cat_pep.grid(row=0, column=1, sticky="w", padx=4, pady=2)
        self.cmb_cat_pep.bind("<<ComboboxSelected>>", lambda e: self.update_categorical_panel_view())

        ttk.Label(ctrl_cat, text="Línea Celular:").grid(row=0, column=2, sticky="w", padx=4, pady=2)
        self.var_cat_cell = tk.StringVar(value="Todas")
        self.cmb_cat_cell = ttk.Combobox(ctrl_cat, textvariable=self.var_cat_cell, state="readonly", width=14, values=["Todas", "THP-1", "PBMCs"])
        self.cmb_cat_cell.grid(row=0, column=3, sticky="w", padx=4, pady=2)
        self.cmb_cat_cell.bind("<<ComboboxSelected>>", lambda e: self.update_categorical_panel_view())

        ttk.Label(ctrl_cat, text="Curva Destacada:").grid(row=0, column=4, sticky="w", padx=4, pady=2)
        self.var_cat_curve = tk.StringVar(value="")
        self.cmb_cat_curve = ttk.Combobox(ctrl_cat, textvariable=self.var_cat_curve, state="readonly", width=28)
        self.cmb_cat_curve.grid(row=0, column=5, sticky="w", padx=4, pady=2)
        self.cmb_cat_curve.bind("<<ComboboxSelected>>", lambda e: self.update_categorical_panel_view())

        ttk.Label(ctrl_cat, text="Visualización Gráfica:", font=("DejaVu Sans", 8, "bold")).grid(row=1, column=0, sticky="w", padx=4, pady=2)
        self.var_cat_vis_mode = tk.StringVar(value="Balance Temporal: % Afectados (V2L) vs % No Afectados (V2R) y Hitos")
        self.cmb_cat_vis = ttk.Combobox(
            ctrl_cat, textvariable=self.var_cat_vis_mode, state="readonly", width=56,
            values=[
                "Balance Temporal: % Afectados (V2L) vs % No Afectados (V2R) y Hitos",
                "Comparativa de Hitos: T_onset, T50 del Efecto y T_efecto_máx + % Máx/Mín por Dosis",
                "Mapa de Calor Espaciotemporal: % Eventos Afectados (V2L) por Compuesto, Dosis y Tiempo",
                "Velocidad Instantánea de Conversión d(%V2L)/dt y Entropía de Shannon H(t)"
            ]
        )
        self.cmb_cat_vis.grid(row=1, column=1, columnspan=3, sticky="ew", padx=4, pady=2)
        self.cmb_cat_vis.bind("<<ComboboxSelected>>", lambda e: self.update_categorical_panel_view())

        btn_cat_box = ttk.Frame(ctrl_cat)
        btn_cat_box.grid(row=1, column=4, columnspan=2, sticky="e", padx=4, pady=2)
        tk.Button(
            btn_cat_box, text="🔄 Recalcular Panel Categórico",
            bg="#004488", fg="white", font=("DejaVu Sans", 8, "bold"),
            padx=8, command=self.recompute_all_models_and_tables
        ).pack(side="left", padx=3)
        tk.Button(
            btn_cat_box, text="📄 Exportar en Reporte Final ➔",
            bg="#117733", fg="white", font=("DejaVu Sans", 8, "bold"),
            padx=8, command=lambda: self.notebook.select(self.tab_finalreport)
        ).pack(side="left", padx=3)

        cat_paned = ttk.PanedWindow(self.tab_catpanel, orient="vertical")
        cat_paned.pack(fill="both", expand=True, padx=10, pady=4)

        self.cat_prop_canvas = NativeCategoricalProportionsCanvas(cat_paned)
        cat_paned.add(self.cat_prop_canvas, weight=2)

        cat_nb = ttk.Notebook(cat_paned)
        cat_paned.add(cat_nb, weight=3)

        self.tree_cat_milestones = self._create_dynamic_treeview(
            cat_nb, " 1. Hitos del Efecto (% Afectados Máx/Mín, T_efecto_máx, T50 Efecto, LT50 y Definición de Afectación) "
        )
        self.tree_cat_windows = self._create_dynamic_treeview(
            cat_nb, " 2. Proporciones por Ventana de Tiempo (Afectados V2L vs No Afectados V2R, IC95%, Odds y Fisher/Chi²) "
        )
        self.tree_cat_tests = self._create_dynamic_treeview(
            cat_nb, " 3. Pruebas Estadísticas (Tendencia Cochran-Armitage, Log-Rank, Isocronas EC50(t) y Contrastes) "
        )

        frm_insights = ttk.Frame(cat_nb, padding=6)
        cat_nb.add(frm_insights, text=" 4. Interpretación Biofísica Automática y Análisis Sugeridos del Fenómeno ")
        self.txt_cat_insights = tk.Text(frm_insights, wrap="word", font=("Courier", 9), bg="#F8F9FA")
        vsb_ins = ttk.Scrollbar(frm_insights, orient="vertical", command=self.txt_cat_insights.yview)
        self.txt_cat_insights.configure(yscrollcommand=vsb_ins.set)
        self.txt_cat_insights.pack(side="left", fill="both", expand=True)
        vsb_ins.pack(side="right", fill="y")

    def update_categorical_panel_view(self):
        if self.df_categorical_summary.empty and self.df_cat_time_windows.empty:
            return

        pep_sel = self.var_cat_pep.get()
        cell_sel = self.var_cat_cell.get()

        sum_sub = self.df_categorical_summary.copy()
        win_sub = self.df_cat_time_windows.copy()
        tst_sub = self.df_cat_hypothesis_tests.copy()

        if pep_sel and pep_sel != "Todos":
            sum_sub = sum_sub[sum_sub['Peptide'] == pep_sel]
            win_sub = win_sub[win_sub['Peptide'] == pep_sel]
            if not tst_sub.empty and 'Peptide' in tst_sub.columns:
                tst_sub = tst_sub[tst_sub['Peptide'].astype(str).str.contains(pep_sel, regex=False)]

        if cell_sel and cell_sel != "Todas":
            sum_sub = sum_sub[sum_sub['Cell_Line'] == cell_sel]
            win_sub = win_sub[win_sub['Cell_Line'] == cell_sel]

        if not sum_sub.empty:
            cids = sorted(sum_sub.apply(lambda r: f"{r['Peptide']} | {r['Cell_Line']} | {r['Concentration_uM']:.2f} µM", axis=1).unique())
            self.cmb_cat_curve['values'] = cids
            if self.var_cat_curve.get() not in cids:
                target_c = self.ent_target.get().strip()
                pref = [c for c in cids if target_c in c]
                self.var_cat_curve.set(pref[-1] if pref else cids[-1])

        self._populate_treeview_from_df(self.tree_cat_milestones, sum_sub)
        self._populate_treeview_from_df(self.tree_cat_windows, win_sub)
        self._populate_treeview_from_df(self.tree_cat_tests, tst_sub)

        self.cat_prop_canvas.set_panel_data(
            self.var_cat_vis_mode.get(),
            win_sub, sum_sub,
            active_label=self.var_cat_curve.get()
        )

    # --------------------------------------------------------------------------
    # PESTAÑA 4: PRE-REPORTE, SELECCIÓN DE MODELOS CINÉTICOS Y ESTADÍSTICA
    # --------------------------------------------------------------------------
    def _build_tab3_prereport(self):
        ctrl_bar = ttk.LabelFrame(
            self.tab_prereport,
            text=" Configuración de Analizadores Cinéticos / Farmacodinámicos y Diseño de Comparación Estadística ",
            padding=8
        )
        ctrl_bar.pack(fill="x", padx=10, pady=5)

        ttk.Label(ctrl_bar, text="Analizador Cinético / Modelo:").grid(row=0, column=0, sticky="w", padx=4, pady=2)
        self.var_kinetic_analyzer = tk.StringVar(value="Allosteric Hill (Regularizado nH <= 6.0)")
        self.cmb_analyzer = ttk.Combobox(
            ctrl_bar, textvariable=self.var_kinetic_analyzer, state="readonly", width=36,
            values=[
                "Allosteric Hill (Regularizado nH <= 6.0)",
                "Michaelis-Menten Hiperbólico (nH = 1.0)",
                "Gompertz Sigmoidal Asimétrico",
                "Weibull / Exponencial Estirado",
                "Bi-Exponencial Secuencial en Tiempo (k_ads + k_lysis)"
            ]
        )
        self.cmb_analyzer.grid(row=0, column=1, sticky="w", padx=4, pady=2)
        self.cmb_analyzer.bind("<<ComboboxSelected>>", lambda e: self.recompute_all_models_and_tables())

        ttk.Label(ctrl_bar, text="Modo de Comparación Estadística:").grid(row=0, column=2, sticky="w", padx=4, pady=2)
        self.var_stat_mode = tk.StringVar(value="Vs. Control Basal (0 µM) y Línea Control Sano")
        self.cmb_stat_mode = ttk.Combobox(
            ctrl_bar, textvariable=self.var_stat_mode, state="readonly", width=36,
            values=[
                "Vs. Control Basal (0 µM) y Línea Control Sano",
                "Entre Compuestos (a Dosis Equivalentes / Orden de Potencia)",
                "Entre Dosis Sucesivas (2x Step-Up dentro de cada Compuesto)",
                "Entre Líneas Celulares (Selectividad por Concentración)"
            ]
        )
        self.cmb_stat_mode.grid(row=0, column=3, sticky="w", padx=4, pady=2)
        self.cmb_stat_mode.bind("<<ComboboxSelected>>", lambda e: self.recompute_all_models_and_tables())

        ttk.Label(ctrl_bar, text="Corrección Múltiple:").grid(row=1, column=0, sticky="w", padx=4, pady=2)
        self.var_fdr_method = tk.StringVar(value="Benjamini-Hochberg (FDR)")
        self.cmb_fdr = ttk.Combobox(
            ctrl_bar, textvariable=self.var_fdr_method, state="readonly", width=28,
            values=["Benjamini-Hochberg (FDR)", "Holm-Bonferroni", "Sin ajuste (Welch p directo)"]
        )
        self.cmb_fdr.grid(row=1, column=1, sticky="w", padx=4, pady=2)
        self.cmb_fdr.bind("<<ComboboxSelected>>", lambda e: self.recompute_all_models_and_tables())

        ttk.Label(ctrl_bar, text="Réplicas Bootstrap (IC 95%):").grid(row=1, column=2, sticky="w", padx=4, pady=2)
        self.ent_bootstrap_n = ttk.Entry(ctrl_bar, width=10)
        self.ent_bootstrap_n.insert(0, "200")
        self.ent_bootstrap_n.grid(row=1, column=3, sticky="w", padx=4, pady=2)

        tk.Button(
            ctrl_bar, text="🔄 Recalcular Modelos y Tablas",
            bg="#004488", fg="white", font=("DejaVu Sans", 9, "bold"),
            padx=10, command=self.recompute_all_models_and_tables
        ).grid(row=0, column=4, padx=6, pady=2)

        tk.Button(
            ctrl_bar, text="❓ Explicación del Modelo y Ecuaciones",
            bg="#6C3483", fg="white", font=("DejaVu Sans", 8, "bold"),
            padx=8, command=lambda: self.open_mathematical_models_help_dialog(initial_tab=0)
        ).grid(row=1, column=4, padx=6, pady=2)

        mid_paned = ttk.PanedWindow(self.tab_prereport, orient="horizontal")
        mid_paned.pack(fill="both", expand=True, padx=10, pady=4)

        tbl_notebook = ttk.Notebook(mid_paned)
        mid_paned.add(tbl_notebook, weight=3)

        self.tree_t0 = self._create_dynamic_treeview(tbl_notebook, "Tabla 0: Prerregistro Pre-Análisis")
        self.tree_t1 = self._create_dynamic_treeview(tbl_notebook, "Tabla 1: Tasa 0-120s")
        self.tree_t2 = self._create_dynamic_treeview(tbl_notebook, "Tabla 2: Tasa 0-600s")
        self.tree_t3 = self._create_dynamic_treeview(tbl_notebook, "Tabla 3: V_active [T_lag-T_max]")
        self.tree_tcat = self._create_dynamic_treeview(tbl_notebook, "Tabla Cat: Cut-Off V2L%, LT50, Entropía y Modo")
        self.tree_t4 = self._create_dynamic_treeview(tbl_notebook, "Tabla 4: Selectividad y Odds Ratio")
        self.tree_t5 = self._create_dynamic_treeview(tbl_notebook, "Tabla 5: Contrastes Dinámicos / Dosis")
        self.tree_t6 = self._create_dynamic_treeview(tbl_notebook, "Tabla 6: Parámetros del Modelo (Hill / Otros)")
        self.tree_ttimebin = self._create_dynamic_treeview(tbl_notebook, "Tabla Bin: Datos por Bin de Tiempo (Diana vs Control)")

        sat_frame = ttk.LabelFrame(
            mid_paned,
            text=" 📐 Panel de Previsualización: Gráfica de la Función, Ecuación de los Modelos y R² del Grupo de Muestras ",
            padding=4
        )
        mid_paned.add(sat_frame, weight=2)

        prev_ctrl = ttk.Frame(sat_frame)
        prev_ctrl.pack(fill="x", padx=2, pady=2)

        ttk.Label(prev_ctrl, text="Grupo de Muestras:", font=("DejaVu Sans", 8, "bold")).pack(side="left", padx=2)
        self.var_preview_group = tk.StringVar(value="Todos los Grupos")
        self.cmb_preview_group = ttk.Combobox(
            prev_ctrl, textvariable=self.var_preview_group, state="readonly", width=22,
            values=["Todos los Grupos"]
        )
        self.cmb_preview_group.pack(side="left", padx=2)
        self.cmb_preview_group.bind("<<ComboboxSelected>>", lambda e: self.recompute_all_models_and_tables())

        ttk.Label(prev_ctrl, text="Vista:", font=("DejaVu Sans", 8, "bold")).pack(side="left", padx=2)
        self.var_preview_mode = tk.StringVar(value="Modelo Seleccionado (Por Grupo)")
        self.cmb_preview_mode = ttk.Combobox(
            prev_ctrl, textvariable=self.var_preview_mode, state="readonly", width=28,
            values=[
                "Modelo Seleccionado (Por Grupo)",
                "Comparar Todos los Modelos en el Grupo"
            ]
        )
        self.cmb_preview_mode.pack(side="left", fill="x", expand=True, padx=2)
        self.cmb_preview_mode.bind("<<ComboboxSelected>>", lambda e: self.recompute_all_models_and_tables())

        self.sat_preview = NativeSaturationPreview(sat_frame)
        self.sat_preview.pack(fill="both", expand=True, pady=2)

        eq_tbl_frame = ttk.LabelFrame(
            sat_frame,
            text=" Ecuación Parametrizada de la Función y Bondad de Ajuste (R² del Grupo de Muestras) ",
            padding=2
        )
        eq_tbl_frame.pack(fill="x", pady=2)

        eq_cols = ("Grupo", "Modelo", "Ecuacion_Evaluada", "R2_Grupo", "R2_Replicas", "R2_Adj", "RMSE")
        self.tree_model_eq_r2 = ttk.Treeview(eq_tbl_frame, columns=eq_cols, show="headings", height=4, selectmode="browse")
        eq_col_specs = [
            ("Grupo", "Grupo de Muestras", 125, "w"),
            ("Modelo", "Modelo", 115, "w"),
            ("Ecuacion_Evaluada", "Ecuación Evaluada de la Función f(x)", 245, "w"),
            ("R2_Grupo", "R² Grupo (Medias)", 95, "center"),
            ("R2_Replicas", "R² Réplicas (Todas)", 98, "center"),
            ("R2_Adj", "R² Ajustado", 75, "center"),
            ("RMSE", "RMSE", 62, "center")
        ]
        for cid, h_txt, w_px, anc in eq_col_specs:
            self.tree_model_eq_r2.heading(cid, text=h_txt)
            self.tree_model_eq_r2.column(cid, width=w_px, anchor=anc, stretch=(cid == "Ecuacion_Evaluada"))

        vsb_eq = ttk.Scrollbar(eq_tbl_frame, orient="vertical", command=self.tree_model_eq_r2.yview)
        hsb_eq = ttk.Scrollbar(eq_tbl_frame, orient="horizontal", command=self.tree_model_eq_r2.xview)
        self.tree_model_eq_r2.configure(yscrollcommand=vsb_eq.set, xscrollcommand=hsb_eq.set)
        self.tree_model_eq_r2.grid(row=0, column=0, sticky="nsew")
        vsb_eq.grid(row=0, column=1, sticky="ns")
        hsb_eq.grid(row=1, column=0, sticky="ew")
        eq_tbl_frame.rowconfigure(0, weight=1)
        eq_tbl_frame.columnconfigure(0, weight=1)

        self.lbl_model_preview_summary = tk.Label(
            sat_frame,
            text="Modelo Activo: Allosteric Hill | Ecuación: V_active(C) = V_max·C^nH / (Kd^nH + C^nH)",
            bg="#EBF5FB", fg="#002B49", font=("DejaVu Sans", 8, "bold"),
            relief="groove", pady=4, justify="left"
        )
        self.lbl_model_preview_summary.pack(fill="x", pady=2)

    def _create_dynamic_treeview(self, parent_notebook, tab_title):
        frame = ttk.Frame(parent_notebook, padding=4)
        parent_notebook.add(frame, text=tab_title)
        tree = ttk.Treeview(frame, show="headings")
        vsb = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        return tree

    def _populate_treeview_from_df(self, tree, df):
        tree.delete(*tree.get_children())
        if df is None or df.empty:
            return
        cols = list(df.columns)
        tree["columns"] = cols
        for c in cols:
            tree.heading(c, text=str(c))
            tree.column(c, width=max(75, min(140, len(str(c)) * 9)), anchor="center")
        for _, row in df.iterrows():
            vals = [
                f"{v:.4f}" if isinstance(v, (float, np.floating)) and abs(v) < 10 and v != 0
                else (f"{v:.2f}" if isinstance(v, (float, np.floating)) else str(v))
                for v in row.values
            ]
            tree.insert("", "end", values=vals)

    # --------------------------------------------------------------------------
    # PESTAÑA 4: REPORTE FINAL Y EXPORTACIÓN (600 DPI + PDF POR CATEGORÍAS)
    # --------------------------------------------------------------------------
    def _build_tab4_finalreport(self):
        top_f = ttk.LabelFrame(
            self.tab_finalreport,
            text=" Configuración de Exportación del Reporte Auditable por Categorías (600 DPI PNG + CSV/TSV + PDF) ",
            padding=10
        )
        top_f.pack(fill="x", padx=12, pady=8)

        ttk.Label(top_f, text="Carpeta de Destino:").grid(row=0, column=0, sticky="w", padx=5, pady=4)
        self.ent_outdir = ttk.Entry(top_f, width=58)
        self.ent_outdir.insert(0, str(Path.cwd() / "Resultados_Cinetica_Calceina"))
        self.ent_outdir.grid(row=0, column=1, sticky="w", padx=5, pady=4)
        ttk.Button(top_f, text="Examinar...", command=self.select_outdir).grid(row=0, column=2, padx=5, pady=4)

        ttk.Label(top_f, text="Título del Reporte PDF:").grid(row=1, column=0, sticky="w", padx=5, pady=4)
        self.ent_rep_title = ttk.Entry(top_f, width=58)
        self.ent_rep_title.insert(0, "The Ultimate Biostatistical Audit Report: Membrane Disruption Analytics")
        self.ent_rep_title.grid(row=1, column=1, sticky="w", padx=5, pady=4)

        ttk.Label(top_f, text="Autor / Laboratorio:").grid(row=2, column=0, sticky="w", padx=5, pady=4)
        self.ent_rep_author = ttk.Entry(top_f, width=58)
        self.ent_rep_author.insert(0, "Flow Cytometry Biophysics & Pharmacodynamics Unit")
        self.ent_rep_author.grid(row=2, column=1, sticky="w", padx=5, pady=4)

        tk.Button(
            top_f,
            text="🚀 GENERAR Y EXPORTAR REPORTE FINAL AUDITABLE (PDF CON IMÁGENES + 600 DPI + CSV/TSV)",
            bg="#117733", fg="white", font=("DejaVu Sans", 10, "bold"),
            padx=16, pady=8, command=self.export_final_report_bundle
        ).grid(row=3, column=0, columnspan=3, pady=10)

        preview_f = ttk.LabelFrame(
            self.tab_finalreport,
            text=" Bitácora Ejecutiva del Reporte Final (El programa permanece abierto tras exportar) ",
            padding=8
        )
        preview_f.pack(fill="both", expand=True, padx=12, pady=6)

        self.txt_report_summary = tk.Text(preview_f, wrap="word", font=("Courier", 9), bg="#F8F9FA")
        vsb = ttk.Scrollbar(preview_f, orient="vertical", command=self.txt_report_summary.yview)
        self.txt_report_summary.configure(yscrollcommand=vsb.set)
        self.txt_report_summary.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")

    def select_outdir(self):
        d = filedialog.askdirectory(title="Seleccionar Carpeta de Salida", parent=self.root)
        if d:
            self.ent_outdir.delete(0, tk.END)
            self.ent_outdir.insert(0, d)

    # --------------------------------------------------------------------------
    # MOTOR DE INGESTIÓN, CUT-OFF DE HISTOGRAMA Y CÁLCULO DE PERMEABILIZACIÓN
    # --------------------------------------------------------------------------
    def _parse_bound_value(self, raw_str):
        txt = str(raw_str).strip().upper()
        if not txt or any(k in txt for k in ("SIN", "INF", "NONE", "MAX")):
            return np.inf
        return float(txt)

    def auto_calibrate_cutoff_from_cc(self):
        """
        Establece automáticamente:
          - Cut-Off 1 (Exclusión V1L) en el percentil 0.5% de las muestras control (CC).
          - Cut-Off 2 (Discriminación V2L/V2R) en el percentil 5.0% de las células cargadas CC.
        """
        if not self.tube_event_cache:
            return
        cc_logs = []
        for meta in self.file_records:
            if meta['Sample_Role'].startswith('CC') and meta['Filepath'] in self.tube_event_cache:
                l_arr = self.tube_event_cache[meta['Filepath']]['log10_fitc']
                if len(l_arr) > 0:
                    cc_logs.append(l_arr)
        if not cc_logs:
            for fp, ev_c in self.tube_event_cache.items():
                mask0 = ev_c['time_s'] <= 0
                if np.any(mask0):
                    cc_logs.append(ev_c['log10_fitc'][mask0])
        if cc_logs:
            all_cc = np.concatenate(cc_logs)
            excl_val = float(self.var_log_excl_v1l.get())
            loaded_cc = all_cc[all_cc >= excl_val] if np.any(all_cc >= excl_val) else all_cc
            new_cut = round(float(np.percentile(loaded_cc, 5.0)), 2)
            self.var_log_cutoff.set(new_cut)
            self.apply_prefilter_by_current_scope(final=True)

    def _get_selected_sample_filters(self):
        """Devuelve las listas (sel_peps, sel_cells, sel_doses) activas en la Caja 2."""
        sel_peps = [self.lst_compounds.get(i) for i in self.lst_compounds.curselection()]
        sel_cells = [self.lst_cells.get(i) for i in self.lst_cells.curselection()]
        sel_doses = [float(self.lst_doses.get(i)) for i in self.lst_doses.curselection()]

        if not sel_peps and not self.df_full.empty:
            sel_peps = list(self.df_full['Peptide'].unique())
        if not sel_cells and not self.df_full.empty:
            sel_cells = list(self.df_full['Cell_Line'].unique())
        if not sel_doses and not self.df_full.empty:
            sel_doses = list(self.df_full['Concentration_uM'].unique())
        return sel_peps, sel_cells, sel_doses

    def apply_prefilter_by_current_scope(self, final=True):
        scope = "selected" if "Discriminado" in self.var_cutoff_scope.get() else "global"
        self.apply_prefilter_and_cutoffs(scope=scope, final=final, preserve_individual=(scope != "global"))

    def apply_new_log_cutoff(self, new_cutoff_val, new_excl_val=None, final=True):
        """Callback invocado al arrastrar cualquiera de las 2 líneas de Cut-Off en el histograma."""
        try:
            self.var_log_cutoff.set(round(float(new_cutoff_val), 2))
            if new_excl_val is not None:
                self.var_log_excl_v1l.set(round(float(new_excl_val), 2))
        except ValueError:
            return
        scope = "selected" if "Discriminado" in self.var_cutoff_scope.get() else "global"
        self.apply_prefilter_and_cutoffs(scope=scope, final=final, preserve_individual=(scope != "global"))

    def _sync_all_samples_with_active_user_cutoffs(self):
        """
        Sincroniza y recalcula todos los tubos desde la caché de eventos aplicando exactamente
        los Cut-Offs configurados por el usuario (ya sean Individuales en la Tabla del Panel 2,
        Discriminados por selección en Caja 2, o Globales en los controles del Panel 2), así como
        los umbrales T_lag y T_max (Manuales o Automáticos) configurados en self.curve_thresholds.
        """
        if not self.tube_event_cache:
            return

        try:
            fsc_min = float(self.ent_fsc_min.get().strip())
            fsc_max = self._parse_bound_value(self.ent_fsc_max.get())
            ssc_min = float(self.ent_ssc_min.get().strip())
            ssc_max = self._parse_bound_value(self.ent_ssc_max.get())
            bin_size = int(float(self.ent_bin_size.get().strip()))
            cells_per_ml = float(self.ent_density.get().strip())
            log_excl_v1l = float(self.var_log_excl_v1l.get())
            log_cutoff = float(self.var_log_cutoff.get())
        except ValueError:
            return

        if log_excl_v1l >= log_cutoff:
            log_excl_v1l = round(log_cutoff - 0.20, 2)
            self.var_log_excl_v1l.set(log_excl_v1l)

        records = []
        tot_raw, tot_v1l_omitted, tot_ret = 0, 0, 0
        used_channel_name = self.var_perm_channel.get()

        for meta in self.file_records:
            fp = meta['Filepath']
            if fp not in self.tube_event_cache:
                continue

            existing_ov = self.sample_filter_overrides.get(fp, {})
            existing_scope = existing_ov.get('scope_label', 'Global')

            # Respetar íntegramente el override de cada muestra en self.sample_filter_overrides
            # (ya sea Individual, Discriminado o Global previamente aplicado por el usuario).
            if existing_ov:
                ov = {
                    'bin_size': int(existing_ov.get('bin_size', bin_size)),
                    'fsc_min': float(existing_ov.get('fsc_min', fsc_min)),
                    'fsc_max': float(existing_ov['fsc_max']) if (existing_ov.get('fsc_max') is not None and np.isfinite(existing_ov.get('fsc_max'))) else fsc_max,
                    'ssc_min': float(existing_ov.get('ssc_min', ssc_min)),
                    'ssc_max': float(existing_ov['ssc_max']) if (existing_ov.get('ssc_max') is not None and np.isfinite(existing_ov.get('ssc_max'))) else ssc_max,
                    'log_excl_v1l': float(existing_ov.get('log_excl_v1l', log_excl_v1l)),
                    'log_cutoff': float(existing_ov.get('log_cutoff', log_cutoff)),
                    'scope_label': existing_scope
                }
            else:
                ov = {
                    'bin_size': bin_size,
                    'fsc_min': fsc_min, 'fsc_max': fsc_max,
                    'ssc_min': ssc_min, 'ssc_max': ssc_max,
                    'log_excl_v1l': log_excl_v1l,
                    'log_cutoff': log_cutoff,
                    'scope_label': 'Global'
                }
            self.sample_filter_overrides[fp] = ov

            try:
                df_t = apply_gate_and_cutoffs_to_cache(
                    fp, self.tube_event_cache[fp],
                    bin_size=ov['bin_size'],
                    fsc_min=ov['fsc_min'], fsc_max=ov['fsc_max'],
                    ssc_min=ov['ssc_min'], ssc_max=ov['ssc_max'],
                    log_excl_v1l=ov['log_excl_v1l'], log_cutoff=ov['log_cutoff'],
                    scope_label=ov['scope_label']
                )
                df_t['Applied_Bin_s'] = int(ov['bin_size'])
                df_t['Is_CC'] = meta['Sample_Role'].startswith('CC')
                used_channel_name = str(df_t['Perm_Channel_Used'].iloc[0])
                tot_raw += int(df_t['Total_Raw_Events'].iloc[0])
                tot_v1l_omitted += int(df_t['Excluded_V1L_Events'].iloc[0])
                tot_ret += int(df_t['Retained_Events'].iloc[0])
                records.append(df_t)
            except Exception as ex:
                print(f"[!] Aviso al sincronizar cutoffs de {meta['Filename']}: {ex}")

        if not records:
            return

        pct_ret = (tot_ret / max(1, tot_raw)) * 100.0
        self.lbl_gate_stats.config(
            text=f"Canal: {used_channel_name} | V1R Analizados: {tot_ret:,} ({pct_ret:.1f}%) | V1L Omitidos: {tot_v1l_omitted:,}"
        )

        self._rebuild_normalized_full_dataframe(records, cells_per_ml)
        self.recompute_kinetics_dataframe()
        self.update_preanalysis_plots()

    def apply_prefilter_and_cutoffs(self, scope="global", final=True, preserve_individual=False):
        """
        Aplica los prefiltros (FSC, SSC, Cut-Off 1 de Exclusión V1L) y el Cut-Off 2 de
        Discriminación (V2L/V2R) de manera Global (a todas las muestras) o Discriminada
        (únicamente a las muestras seleccionadas en la Caja 2), recalculando normalización y cinética.
        """
        if not self.tube_event_cache:
            return

        try:
            fsc_min = float(self.ent_fsc_min.get().strip())
            fsc_max = self._parse_bound_value(self.ent_fsc_max.get())
            ssc_min = float(self.ent_ssc_min.get().strip())
            ssc_max = self._parse_bound_value(self.ent_ssc_max.get())
            bin_size = int(float(self.ent_bin_size.get().strip()))
            cells_per_ml = float(self.ent_density.get().strip())
            log_excl_v1l = float(self.var_log_excl_v1l.get())
            log_cutoff = float(self.var_log_cutoff.get())
        except ValueError:
            return

        if log_excl_v1l >= log_cutoff:
            log_excl_v1l = round(log_cutoff - 0.20, 2)
            self.var_log_excl_v1l.set(log_excl_v1l)

        sel_peps, sel_cells, sel_doses = self._get_selected_sample_filters()
        if scope == "selected":
            self.var_cutoff_scope.set("Discriminado (Muestras seleccionadas)")
        else:
            self.var_cutoff_scope.set("Global (Todas las muestras)")

        records = []
        tot_raw, tot_v1l_omitted, tot_ret = 0, 0, 0
        used_channel_name = self.var_perm_channel.get()

        for meta in self.file_records:
            fp = meta['Filepath']
            if fp not in self.tube_event_cache:
                continue

            is_selected_sample = (
                (meta['Compound'] in sel_peps) and
                (meta['Cell_Line'] in sel_cells) and
                any(np.isclose(meta['Concentration_uM'], d, atol=1e-2) for d in sel_doses)
            )
            existing_scope = self.sample_filter_overrides.get(fp, {}).get('scope_label', 'Global')

            if scope == "global" and not (preserve_individual and existing_scope == "Individual"):
                self.sample_filter_overrides[fp] = {
                    'bin_size': bin_size,
                    'fsc_min': fsc_min, 'fsc_max': fsc_max,
                    'ssc_min': ssc_min, 'ssc_max': ssc_max,
                    'log_excl_v1l': log_excl_v1l,
                    'log_cutoff': log_cutoff,
                    'scope_label': 'Global'
                }
            elif fp not in self.sample_filter_overrides or is_selected_sample:
                lbl_sc = "Discriminado" if (scope == "selected" and is_selected_sample) else existing_scope
                self.sample_filter_overrides[fp] = {
                    'bin_size': bin_size,
                    'fsc_min': fsc_min, 'fsc_max': fsc_max,
                    'ssc_min': ssc_min, 'ssc_max': ssc_max,
                    'log_excl_v1l': log_excl_v1l,
                    'log_cutoff': log_cutoff,
                    'scope_label': lbl_sc
                }

            ov = self.sample_filter_overrides[fp]
            try:
                df_t = apply_gate_and_cutoffs_to_cache(
                    fp, self.tube_event_cache[fp],
                    bin_size=ov['bin_size'],
                    fsc_min=ov['fsc_min'], fsc_max=ov['fsc_max'],
                    ssc_min=ov['ssc_min'], ssc_max=ov['ssc_max'],
                    log_excl_v1l=ov['log_excl_v1l'], log_cutoff=ov['log_cutoff'],
                    scope_label=ov['scope_label']
                )
                df_t['Applied_Bin_s'] = int(ov['bin_size'])
                df_t['Is_CC'] = meta['Sample_Role'].startswith('CC')
                used_channel_name = str(df_t['Perm_Channel_Used'].iloc[0])
                tot_raw += int(df_t['Total_Raw_Events'].iloc[0])
                tot_v1l_omitted += int(df_t['Excluded_V1L_Events'].iloc[0])
                tot_ret += int(df_t['Retained_Events'].iloc[0])
                records.append(df_t)
            except Exception as ex:
                print(f"[!] Aviso al filtrar {meta['Filename']}: {ex}")

        if not records:
            return

        pct_ret = (tot_ret / max(1, tot_raw)) * 100.0
        self.lbl_gate_stats.config(
            text=f"Canal: {used_channel_name} | V1R Analizados: {tot_ret:,} ({pct_ret:.1f}%) | V1L Omitidos: {tot_v1l_omitted:,}"
        )

        self._rebuild_normalized_full_dataframe(records, cells_per_ml)
        self.recompute_kinetics_dataframe()
        self.update_preanalysis_plots()
        if final:
            self.recompute_all_models_and_tables()

    def _rebuild_normalized_full_dataframe(self, records, cells_per_ml):
        df_raw = pd.concat(records, ignore_index=True)
        df_meta = pd.DataFrame(self.file_records).drop_duplicates(subset=['Filepath'])
        df_merged = pd.merge(df_raw, df_meta, on='Filepath', how='left')
        df_merged.rename(columns={'Compound': 'Peptide'}, inplace=True)

        df_cc = df_merged[df_merged['Is_CC']].copy()
        df_exp = df_merged[~df_merged['Is_CC']].copy()

        if not df_cc.empty:
            basals_exact = df_cc.groupby(['Cell_Line', 'Peptide', 'Concentration_uM', 'Replicate']).agg(
                FITC_Basal=('Median_FITC', 'median'),
                FSC_Basal=('Mean_FSC', 'mean'),
                SSC_Basal=('Mean_SSC', 'mean')
            ).reset_index()
            df_exp = pd.merge(df_exp, basals_exact, on=['Cell_Line', 'Peptide', 'Concentration_uM', 'Replicate'], how='left')

            basals_rep = df_cc.groupby(['Cell_Line', 'Peptide', 'Replicate']).agg(
                FITC_Basal_Rep=('Median_FITC', 'median'),
                FSC_Basal_Rep=('Mean_FSC', 'mean'),
                SSC_Basal_Rep=('Mean_SSC', 'mean')
            ).reset_index()
            df_exp = pd.merge(df_exp, basals_rep, on=['Cell_Line', 'Peptide', 'Replicate'], how='left')
            df_exp['FITC_Basal'] = df_exp['FITC_Basal'].fillna(df_exp['FITC_Basal_Rep'])
            df_exp['FSC_Basal'] = df_exp['FSC_Basal'].fillna(df_exp['FSC_Basal_Rep'])
            df_exp['SSC_Basal'] = df_exp['SSC_Basal'].fillna(df_exp['SSC_Basal_Rep'])
            df_exp.drop(columns=['FITC_Basal_Rep', 'FSC_Basal_Rep', 'SSC_Basal_Rep'], inplace=True)
        else:
            df_exp['FITC_Basal'] = np.nan
            df_exp['FSC_Basal'] = np.nan
            df_exp['SSC_Basal'] = np.nan

        df_exp['FITC_Basal'] = df_exp['FITC_Basal'].fillna(df_exp.groupby('Filepath')['Median_FITC'].transform('first'))
        df_exp['FSC_Basal'] = df_exp['FSC_Basal'].fillna(df_exp.groupby('Filepath')['Mean_FSC'].transform('first'))
        df_exp['SSC_Basal'] = df_exp['SSC_Basal'].fillna(df_exp.groupby('Filepath')['Mean_SSC'].transform('first'))

        df_exp['Retention_Pct'] = (df_exp['Median_FITC'] / df_exp['FITC_Basal']) * 100.0
        df_exp['Permeabilization_Pct'] = np.clip(100.0 - df_exp['Retention_Pct'], 0.0, 100.0)
        df_exp['Norm_FSC_Pct'] = (df_exp['Mean_FSC'] / df_exp['FSC_Basal']) * 100.0
        df_exp['Norm_SSC_Pct'] = (df_exp['Mean_SSC'] / df_exp['SSC_Basal']) * 100.0

        if not df_cc.empty:
            df_cc_plot = df_cc.copy()
            df_cc_plot['Retention_Pct'] = 100.0
            df_cc_plot['Permeabilization_Pct'] = 0.0
            df_cc_plot['Norm_FSC_Pct'] = 100.0
            df_cc_plot['Norm_SSC_Pct'] = 100.0
            self.df_full = pd.concat([df_cc_plot, df_exp], ignore_index=True)
        else:
            self.df_full = df_exp.copy()

        self.df_full['Molecules_Per_Cell'] = (self.df_full['Concentration_uM'] * 1e-9 * AVOGADRO) / cells_per_ml
        self.df_full['Dose_Label'] = self.df_full['Concentration_uM'].apply(lambda c: f"{c:.2f} µM")

    def run_preanalysis_ingestion(self):
        if not self.file_records:
            messagebox.showwarning(
                "Sin Archivos",
                "Primero sube tus archivos .fcs o haz clic en '🧪 Generar y Cargar Set Demo (.fcs)'.",
                parent=self.root
            )
            return

        try:
            fsc_min = float(self.ent_fsc_min.get().strip())
            fsc_max = self._parse_bound_value(self.ent_fsc_max.get())
            ssc_min = float(self.ent_ssc_min.get().strip())
            ssc_max = self._parse_bound_value(self.ent_ssc_max.get())
            bin_size = int(float(self.ent_bin_size.get().strip()))
            cells_per_ml = float(self.ent_density.get().strip())
            perm_chan = self.var_perm_channel.get().strip() or "FITC (Por defecto)"
            log_excl_v1l = float(self.var_log_excl_v1l.get())
            log_cutoff = float(self.var_log_cutoff.get())
        except ValueError as e:
            messagebox.showerror("Error de Parámetros", f"Verifica que los umbrales FSC/SSC, Bin y Cut-Offs sean numéricos:\n{e}", parent=self.root)
            return

        records = []
        self.tube_event_cache.clear()
        self.sample_filter_overrides.clear()
        tot_raw, tot_v1l, tot_ret = 0, 0, 0
        used_channel_name = perm_chan
        for meta in self.file_records:
            is_cc = meta['Sample_Role'].startswith('CC')
            try:
                df_t, ev_cache = extract_tube_kinetics(
                    meta['Filepath'], is_cc=is_cc, bin_size=bin_size,
                    fsc_min=fsc_min, fsc_max=fsc_max, ssc_min=ssc_min, ssc_max=ssc_max,
                    perm_channel=perm_chan, log_excl_v1l=log_excl_v1l, log_cutoff=log_cutoff
                )
                self.tube_event_cache[meta['Filepath']] = ev_cache
                self.sample_filter_overrides[meta['Filepath']] = {
                    'bin_size': bin_size,
                    'fsc_min': fsc_min, 'fsc_max': fsc_max,
                    'ssc_min': ssc_min, 'ssc_max': ssc_max,
                    'log_excl_v1l': log_excl_v1l,
                    'log_cutoff': log_cutoff,
                    'scope_label': 'Global'
                }
                df_t['Is_CC'] = is_cc
                used_channel_name = str(df_t['Perm_Channel_Used'].iloc[0])
                tot_raw += int(df_t['Total_Raw_Events'].iloc[0])
                tot_v1l += int(df_t['Excluded_V1L_Events'].iloc[0])
                tot_ret += int(df_t['Retained_Events'].iloc[0])
                records.append(df_t)
            except Exception as ex:
                print(f"[!] Aviso en {meta['Filename']}: {ex}")

        if not records:
            messagebox.showerror("Sin eventos viables", "Ningún archivo superó el cut-off de FSC/SSC/V1L actual.", parent=self.root)
            return

        pct_ret = (tot_ret / max(1, tot_raw)) * 100.0
        self.lbl_gate_stats.config(
            text=f"Canal: {used_channel_name} | V1R Analizados: {tot_ret:,} ({pct_ret:.1f}%) | V1L Omitidos: {tot_v1l:,}"
        )

        self._rebuild_normalized_full_dataframe(records, cells_per_ml)
        self._populate_preanalysis_selectors()
        self.auto_calibrate_cutoff_from_cc()
        self.recompute_kinetics_dataframe()
        self.on_active_curve_selected()
        self.recompute_all_models_and_tables()
        self.notebook.select(self.tab_preanalysis)
        self.lbl_global_status.config(text=f"Estado: Pre-Análisis Activo [{used_channel_name}] ({tot_ret:,} eventos V1R analizados)")

    def _populate_preanalysis_selectors(self):
        self.lst_compounds.delete(0, tk.END)
        for p in sorted(self.df_full['Peptide'].unique()):
            self.lst_compounds.insert(tk.END, p)
        self.lst_compounds.selection_set(0, tk.END)

        self.lst_cells.delete(0, tk.END)
        for c in sorted(self.df_full['Cell_Line'].unique()):
            self.lst_cells.insert(tk.END, c)
        self.lst_cells.selection_set(0, tk.END)

        self.lst_doses.delete(0, tk.END)
        for d in sorted(self.df_full['Concentration_uM'].unique()):
            self.lst_doses.insert(tk.END, f"{d:.2f}")
        self.lst_doses.selection_set(0, tk.END)

        curves = []
        for (pep, cell, conc), grp in self.df_full.groupby(['Peptide', 'Cell_Line', 'Concentration_uM']):
            if grp['Time_Window'].max() > 0:
                curves.append(f"{pep} | {cell} | {conc:.2f} µM")
        self.cmb_active_curve['values'] = curves
        if curves:
            target_c = self.ent_target.get().strip()
            pref = [c for c in curves if target_c in c]
            self.cmb_active_curve.set(pref[-1] if pref else curves[0])

        if hasattr(self, 'cmb_cat_pep'):
            peps_all = ["Todos"] + sorted(list(self.df_full['Peptide'].unique()))
            self.cmb_cat_pep['values'] = peps_all
            if self.var_cat_pep.get() not in peps_all:
                self.var_cat_pep.set("Todos")
        if hasattr(self, 'cmb_cat_cell'):
            cells_all = ["Todas"] + sorted(list(self.df_full['Cell_Line'].unique()))
            self.cmb_cat_cell['values'] = cells_all
            if self.var_cat_cell.get() not in cells_all:
                self.var_cat_cell.set("Todas")

    def select_all_filters(self):
        self.lst_compounds.selection_set(0, tk.END)
        self.lst_cells.selection_set(0, tk.END)
        self.lst_doses.selection_set(0, tk.END)
        self.update_preanalysis_plots()

    def _parse_active_curve_key(self):
        val = self.cmb_active_curve.get().strip()
        if not val or "|" not in val:
            return None
        parts = [p.strip() for p in val.split("|")]
        if len(parts) != 3:
            return None
        conc = float(parts[2].replace("µM", "").strip())
        return parts[0], parts[1], conc

    def recompute_kinetics_dataframe(self):
        if self.df_full.empty:
            return
        control_cell = self.ent_control.get().strip()
        tube_params = []

        for keys, grp in self.df_full.groupby(['Peptide', 'Cell_Line', 'Concentration_uM', 'Replicate']):
            pep, cell, conc, rep = keys
            if grp['Time_Window'].max() <= 0:
                continue

            grp_sorted = grp.sort_values('Time_Window')
            t_arr = grp_sorted['Time_Window'].values
            p_arr = grp_sorted['Permeabilization_Pct'].values
            v2l_arr = grp_sorted['Affected_V2L_Pct'].values if 'Affected_V2L_Pct' in grp_sorted else np.zeros_like(p_arr)
            ev_arr = grp_sorted['Events'].values if 'Events' in grp_sorted else np.full_like(p_arr, 100, dtype=int)
            v2l_ev_arr = grp_sorted['Affected_V2L_Events'].values if 'Affected_V2L_Events' in grp_sorted else np.round(ev_arr * v2l_arr / 100.0).astype(int)
            sarle_bc = float(grp_sorted['Sarle_Bimodality_BC'].iloc[0]) if 'Sarle_Bimodality_BC' in grp_sorted else 0.33
            excl_s = float(grp_sorted['Exclusion_s'].iloc[0])

            sub_120 = grp_sorted[(grp_sorted['Time_Window'] >= 0) & (grp_sorted['Time_Window'] <= 120)]
            v0_120 = max(0.0, float(stats.linregress(sub_120['Time_Window'], sub_120['Permeabilization_Pct']).slope)) if len(sub_120) >= 3 else 0.0

            sub_600 = grp_sorted[(grp_sorted['Time_Window'] >= 0) & (grp_sorted['Time_Window'] <= 600)]
            v0_600 = max(0.0, float(stats.linregress(sub_600['Time_Window'], sub_600['Permeabilization_Pct']).slope)) if len(sub_600) >= 3 else 0.0

            is_ctrl = (cell == control_cell)
            curve_key = next((k for k in self.curve_thresholds if k[0] == pep and k[1] == cell and np.isclose(k[2], conc, atol=1e-2)), None)
            manual_tlag, manual_tmax = None, None
            if curve_key and self.curve_thresholds[curve_key]['mode'] == 'manual':
                manual_tlag = self.curve_thresholds[curve_key]['t_lag']
                manual_tmax = self.curve_thresholds[curve_key]['t_max']

            t_lag, t_max, v_active, max_p = extract_slope_based_parameters(
                t_arr, p_arr, is_control_cell=is_ctrl, start_time=excl_s,
                manual_tlag=manual_tlag, manual_tmax=manual_tmax,
                v2l_arr_full=v2l_arr
            )

            exp_mask = t_arr >= 0
            auc_val = trapezoid(y=p_arr[exp_mask], x=t_arr[exp_mask]) if np.sum(exp_mask) >= 2 else 0.0
            norm_auc = (auc_val / (600.0 * 100.0)) * 100.0
            mol_cell = float(grp['Molecules_Per_Cell'].iloc[0])

            # Métricas categóricas integrales basadas en Cut-Off del Histograma (V2L vs V2R) y ventanas de tiempo
            cat_metrics = compute_categorical_kinetics_metrics(
                t_arr, v2l_arr, p_arr, start_time=excl_s,
                events_arr=ev_arr, v2l_events_arr=v2l_ev_arr, sarle_bc=sarle_bc
            )

            # Cambios morfológicos en dispersión frontal (FSC) y lateral (SSC)
            fsc_end = float(grp_sorted[exp_mask]['Norm_FSC_Pct'].iloc[-3:].mean()) if np.sum(exp_mask) >= 3 else 100.0
            ssc_end = float(grp_sorted[exp_mask]['Norm_SSC_Pct'].iloc[-3:].mean()) if np.sum(exp_mask) >= 3 else 100.0

            # Metadatos de pre-registro y Cut-Offs específicos aplicados a este tubo experimental por el usuario
            grp_exp = grp_sorted[exp_mask] if np.any(exp_mask) else grp_sorted
            applied_v1l = float(grp_exp['Applied_V1L_Excl_Log10'].iloc[0]) if 'Applied_V1L_Excl_Log10' in grp_exp else float(self.var_log_excl_v1l.get())
            applied_v2l = float(grp_exp['Applied_V2L_Cutoff_Log10'].iloc[0]) if 'Applied_V2L_Cutoff_Log10' in grp_exp else float(self.var_log_cutoff.get())
            applied_scope = str(grp_exp['Applied_Scope'].iloc[0]) if 'Applied_Scope' in grp_exp else "Global"
            applied_bin_s = int(grp_exp['Applied_Bin_s'].iloc[0]) if 'Applied_Bin_s' in grp_exp else int(float(self.ent_bin_size.get().strip() or 15))
            applied_fsc_min = float(grp_exp['Applied_FSC_Min'].iloc[0]) if 'Applied_FSC_Min' in grp_exp else float(self.ent_fsc_min.get().strip() or 50000.0)
            applied_fsc_max = float(grp_exp['Applied_FSC_Max'].iloc[0]) if 'Applied_FSC_Max' in grp_exp else self._parse_bound_value(self.ent_fsc_max.get())
            applied_ssc_min = float(grp_exp['Applied_SSC_Min'].iloc[0]) if 'Applied_SSC_Min' in grp_exp else float(self.ent_ssc_min.get().strip() or 10000.0)
            applied_ssc_max = float(grp_exp['Applied_SSC_Max'].iloc[0]) if 'Applied_SSC_Max' in grp_exp else self._parse_bound_value(self.ent_ssc_max.get())
            raw_ev = int(grp_exp['Total_Raw_Events'].iloc[0]) if 'Total_Raw_Events' in grp_exp else 0
            omitted_v1l = int(grp_exp['Excluded_V1L_Events'].iloc[0]) if 'Excluded_V1L_Events' in grp_exp else 0
            retained_v1r = int(grp_exp['Retained_Events'].iloc[0]) if 'Retained_Events' in grp_exp else 0
            tlag_mode_lbl = "Manual" if (manual_tlag is not None) else "Auto"

            tube_params.append({
                'Peptide': pep, 'Cell_Line': cell, 'Concentration_uM': conc,
                'Dose_Label': f"{conc:.2f} µM", 'Molecules_Per_Cell': mol_cell,
                'Replicate': rep, 'V0_Early': v0_120, 'V0_600': v0_600,
                'V0_Active': v_active, 'T_lag': t_lag, 'T_max': t_max,
                'Normalized_AUC_Pct': norm_auc, 'Max_Perm': max_p,
                'Basal_V2L_Pct': cat_metrics['Basal_V2L_Pct'],
                'Min_V2L_Pct': cat_metrics['Min_V2L_Pct'],
                'Max_V2L_Pct': cat_metrics['Max_V2L_Pct'],
                'End_V2L_Pct': cat_metrics['End_V2L_Pct'],
                'Min_Unaffected_V2R_Pct': cat_metrics['Min_Unaffected_V2R_Pct'],
                'Max_Unaffected_V2R_Pct': cat_metrics['Max_Unaffected_V2R_Pct'],
                'T_Effect_Min_s': cat_metrics['T_Effect_Min_s'],
                'T_Effect_Max_s': cat_metrics['T_Effect_Max_s'],
                'T10_Effect_s': cat_metrics['T10_Effect_s'],
                'T50_Effect_s': cat_metrics['T50_Effect_s'],
                'T90_Effect_s': cat_metrics['T90_Effect_s'],
                'Span_T10_T90_s': cat_metrics['Span_T10_T90_s'],
                'LT50_s': cat_metrics['LT50_s'],
                'T_Peak_Conversion_s': cat_metrics['T_Peak_Conversion_s'],
                'V_Cat_Max_Inst': cat_metrics['V_Cat_Max_Inst'],
                'T_Stat_Onset_s': cat_metrics['T_Stat_Onset_s'],
                'AUC_V2L_Pct': cat_metrics['AUC_V2L_Pct'],
                'V_Cat_Active': cat_metrics['V_Cat_Active'],
                'Max_Entropy_Bits': cat_metrics['Max_Entropy_Bits'],
                'Hazard_Rate_s': cat_metrics['Hazard_Rate_s'],
                'Sync_Index': cat_metrics['Sync_Index'],
                'Bimodality_BC': cat_metrics['Bimodality_BC'],
                'Cohen_h_Max': cat_metrics['Cohen_h_Max'],
                'OR_vs_Basal': cat_metrics['OR_vs_Basal'],
                'RR_vs_Basal': cat_metrics['RR_vs_Basal'],
                'Chi2_Max': cat_metrics['Chi2_Max'],
                'p_Fisher_Max': cat_metrics['p_Fisher_Max'],
                'Affectation_State': cat_metrics['Affectation_State'],
                'Mode_Delta': cat_metrics['Mode_Delta'],
                'Mechanism_Class': cat_metrics['Mechanism_Class'],
                'End_Norm_FSC_Pct': fsc_end,
                'End_Norm_SSC_Pct': ssc_end,
                'Filter_Scope': applied_scope,
                'V1L_Excl_Log10': applied_v1l,
                'V2L_Cutoff_Log10': applied_v2l,
                'Applied_Bin_s': applied_bin_s,
                'Applied_FSC_Min': applied_fsc_min,
                'Applied_FSC_Max': applied_fsc_max,
                'Applied_SSC_Min': applied_ssc_min,
                'Applied_SSC_Max': applied_ssc_max,
                'Tlag_Mode': tlag_mode_lbl,
                'Exclusion_s': excl_s,
                'Raw_Events': raw_ev,
                'V1L_Omitted': omitted_v1l,
                'V1R_Analyzed': retained_v1r
            })

        self.df_kinetics = pd.DataFrame(tube_params)

        if not self.df_kinetics.empty:
            for (pep, cell, conc), grp in self.df_kinetics.groupby(['Peptide', 'Cell_Line', 'Concentration_uM']):
                existing_k = next((k for k in self.curve_thresholds if k[0] == pep and k[1] == cell and np.isclose(k[2], conc, atol=1e-2)), None)
                if not existing_k:
                    self.curve_thresholds[(pep, cell, float(conc))] = {
                        'mode': 'auto',
                        't_lag': float(grp['T_lag'].mean()),
                        't_max': float(grp['T_max'].mean())
                    }
                elif self.curve_thresholds[existing_k]['mode'] == 'auto':
                    self.curve_thresholds[existing_k]['t_lag'] = float(grp['T_lag'].mean())
                    self.curve_thresholds[existing_k]['t_max'] = float(grp['T_max'].mean())

    def on_active_curve_selected(self, event=None):
        key = self._parse_active_curve_key()
        if not key:
            return
        match_k = next((k for k in self.curve_thresholds if k[0] == key[0] and k[1] == key[1] and np.isclose(k[2], key[2], atol=1e-2)), None)
        if match_k:
            info = self.curve_thresholds[match_k]
            self.var_thresh_mode.set(info['mode'])
            self.var_tlag_slider.set(round(info['t_lag']))
            self.var_tmax_slider.set(round(info['t_max']))
        self.update_preanalysis_plots()

    def on_threshold_mode_toggle(self):
        key = self._parse_active_curve_key()
        if not key:
            return
        match_k = next((k for k in self.curve_thresholds if k[0] == key[0] and k[1] == key[1] and np.isclose(k[2], key[2], atol=1e-2)), key)
        mode = self.var_thresh_mode.get()
        if match_k not in self.curve_thresholds:
            self.curve_thresholds[match_k] = {'mode': mode, 't_lag': 60.0, 't_max': 300.0}
        else:
            self.curve_thresholds[match_k]['mode'] = mode
        self.recompute_kinetics_dataframe()
        if mode == 'auto' and match_k in self.curve_thresholds:
            self.var_tlag_slider.set(round(self.curve_thresholds[match_k]['t_lag']))
            self.var_tmax_slider.set(round(self.curve_thresholds[match_k]['t_max']))
        self.update_preanalysis_plots()
        self.recompute_all_models_and_tables()

    def reset_active_curve_to_auto_tlag_tmax(self):
        """Fuerza el recálculo automático de T_lag y T_max para la curva activa."""
        self.var_thresh_mode.set("auto")
        self.on_threshold_mode_toggle()

    def reset_all_curves_to_auto_tlag_tmax(self):
        """Recalcula automáticamente T_lag y T_max para todas las curvas cargadas."""
        for k in list(self.curve_thresholds.keys()):
            self.curve_thresholds[k]['mode'] = 'auto'
        self.var_thresh_mode.set("auto")
        self.recompute_kinetics_dataframe()
        key = self._parse_active_curve_key()
        if key:
            match_k = next((k for k in self.curve_thresholds if k[0] == key[0] and k[1] == key[1] and np.isclose(k[2], key[2], atol=1e-2)), None)
            if match_k:
                self.var_tlag_slider.set(round(self.curve_thresholds[match_k]['t_lag']))
                self.var_tmax_slider.set(round(self.curve_thresholds[match_k]['t_max']))
        self.update_preanalysis_plots()
        self.recompute_all_models_and_tables()

    def on_slider_thresh_change(self, event=None, final=False):
        key = self._parse_active_curve_key()
        if not key:
            return
        tlag = float(self.var_tlag_slider.get())
        tmax = float(self.var_tmax_slider.get())
        if tmax <= tlag:
            tmax = min(600.0, tlag + 15.0)
            self.var_tmax_slider.set(tmax)
        match_k = next((k for k in self.curve_thresholds if k[0] == key[0] and k[1] == key[1] and np.isclose(k[2], key[2], atol=1e-2)), key)
        self.var_thresh_mode.set("manual")
        self.curve_thresholds[match_k] = {'mode': 'manual', 't_lag': tlag, 't_max': tmax}
        self.recompute_kinetics_dataframe()
        self.update_preanalysis_plots()
        if final:
            self.recompute_all_models_and_tables()

    def on_canvas_drag_threshold(self, new_tlag, new_tmax, final=False):
        key = self._parse_active_curve_key()
        if not key:
            return
        match_k = next((k for k in self.curve_thresholds if k[0] == key[0] and k[1] == key[1] and np.isclose(k[2], key[2], atol=1e-2)), key)
        self.var_thresh_mode.set("manual")
        self.curve_thresholds[match_k] = {'mode': 'manual', 't_lag': float(new_tlag), 't_max': float(new_tmax)}
        self.var_tlag_slider.set(round(new_tlag))
        self.var_tmax_slider.set(round(new_tmax))
        self.recompute_kinetics_dataframe()
        self.update_preanalysis_plots()
        if final:
            self.recompute_all_models_and_tables()

    def update_preanalysis_plots(self):
        if self.df_full.empty:
            return

        sel_peps = [self.lst_compounds.get(i) for i in self.lst_compounds.curselection()]
        sel_cells = [self.lst_cells.get(i) for i in self.lst_cells.curselection()]
        sel_doses = [float(self.lst_doses.get(i)) for i in self.lst_doses.curselection()]

        if not sel_peps:
            sel_peps = list(self.df_full['Peptide'].unique())
        if not sel_cells:
            sel_cells = list(self.df_full['Cell_Line'].unique())
        if not sel_doses:
            sel_doses = list(self.df_full['Concentration_uM'].unique())

        sub = self.df_full[
            (self.df_full['Peptide'].isin(sel_peps)) &
            (self.df_full['Cell_Line'].isin(sel_cells)) &
            (np.isclose(self.df_full['Concentration_uM'].values[:, None], sel_doses, atol=1e-2).any(axis=1))
        ].copy()

        if sub.empty:
            self.multi_plot.set_data([], None, None, "Sin datos seleccionados", self.var_scatter_norm.get())
            return

        sub['Curve_ID'] = sub.apply(lambda r: f"{r['Peptide']} | {r['Cell_Line']} | {r['Concentration_uM']:.2f} µM", axis=1)
        unique_curves = sorted(sub['Curve_ID'].unique())
        control_cell = self.ent_control.get().strip()

        use_norm_scatter = self.var_scatter_norm.get()
        fsc_col_plot = 'Norm_FSC_Pct' if use_norm_scatter else 'Mean_FSC'
        ssc_col_plot = 'Norm_SSC_Pct' if use_norm_scatter else 'Mean_SSC'

        series_list = []
        for idx, cid in enumerate(unique_curves):
            c_df = sub[sub['Curve_ID'] == cid]
            agg = c_df.groupby('Time_Window').agg(
                perm_m=('Permeabilization_Pct', 'mean'), perm_se=('Permeabilization_Pct', 'sem'),
                v2l_m=('Affected_V2L_Pct', 'mean'), v2l_se=('Affected_V2L_Pct', 'sem'),
                fsc_m=(fsc_col_plot, 'mean'), fsc_se=(fsc_col_plot, 'sem'),
                ssc_m=(ssc_col_plot, 'mean'), ssc_se=(ssc_col_plot, 'sem')
            ).reset_index()

            series_list.append({
                'label': cid,
                'color': CURVE_PALETTE[idx % len(CURVE_PALETTE)],
                'is_ctrl': f"| {control_cell} |" in cid,
                't': agg['Time_Window'].values,
                'perm_m': agg['perm_m'].values,
                'perm_se': agg['perm_se'].fillna(0).values,
                'v2l_m': agg['v2l_m'].values,
                'v2l_se': agg['v2l_se'].fillna(0).values,
                'fsc_m': agg['fsc_m'].values,
                'fsc_se': agg['fsc_se'].fillna(0).values,
                'ssc_m': agg['ssc_m'].values,
                'ssc_se': agg['ssc_se'].fillna(0).values
            })

        active_key = self._parse_active_curve_key()
        active_title_str = self.cmb_active_curve.get()
        tlag, tmax = None, None
        v2l_summary_str = ""
        if active_key and not self.df_kinetics.empty:
            match_k = next((k for k in self.curve_thresholds if k[0] == active_key[0] and k[1] == active_key[1] and np.isclose(k[2], active_key[2], atol=1e-2)), None)
            if match_k:
                tlag = self.curve_thresholds[match_k]['t_lag']
                tmax = self.curve_thresholds[match_k]['t_max']
                k_sub = self.df_kinetics[
                    (self.df_kinetics['Peptide'] == match_k[0]) &
                    (self.df_kinetics['Cell_Line'] == match_k[1]) &
                    (np.isclose(self.df_kinetics['Concentration_uM'], match_k[2], atol=1e-2))
                ]
                if not k_sub.empty:
                    v_act_m = k_sub['V0_Active'].mean()
                    v2l_max_m = k_sub['Max_V2L_Pct'].mean()
                    t50_eff_m = k_sub['T50_Effect_s'].mean()
                    tmax_eff_m = k_sub['T_Effect_Max_s'].mean()
                    mech_str = k_sub['Mechanism_Class'].iloc[0]
                    mode_lbl = "AUTO" if self.curve_thresholds[match_k]['mode'] == 'auto' else "MANUAL"
                    self.lbl_active_metrics.config(
                        text=f"[{mode_lbl}] T_lag: {tlag:.0f}s | T_max: {tmax:.0f}s | V_act: {v_act_m:.3f}%/s\n"
                             f"V2L Máx: {v2l_max_m:.1f}% (en {tmax_eff_m:.0f}s) | T50 Efecto: {t50_eff_m:.0f}s | {mech_str}"
                    )
                    v2l_summary_str = f"V2L Máx: {v2l_max_m:.1f}% | T50: {t50_eff_m:.0f}s"

        self.multi_plot.set_data(
            series_list, tlag, tmax, active_title_str,
            use_norm_scatter, self.var_perm_display_mode.get()
        )
        self._update_histogram_canvas(sel_peps, sel_cells, sel_doses, unique_curves, v2l_summary_str)
        self._refresh_sample_cutoffs_table()

    # --------------------------------------------------------------------------
    # TABLA INTERACTIVA DE CUT-OFFS Y VENTANAS POR MUESTRA EN EL PANEL 2
    # --------------------------------------------------------------------------
    def _refresh_sample_cutoffs_table(self):
        if not hasattr(self, 'tree_sample_cutoffs'):
            return
        self._updating_cutoff_tree = True
        try:
            for item in self.tree_sample_cutoffs.get_children():
                self.tree_sample_cutoffs.delete(item)
            if self.df_full.empty:
                return

            active_curve_str = self.cmb_active_curve.get().strip()
            selected_iid = None

            for (pep, cell, conc), grp in self.df_full.groupby(['Peptide', 'Cell_Line', 'Concentration_uM']):
                if grp['Time_Window'].max() <= 0:
                    continue
                cid_label = f"{pep} | {cell} | {conc:.2f} µM"
                iid = f"{pep}||{cell}||{conc:.4f}"

                # Extraer override de cut-offs de los archivos de esta muestra
                sample_fps = [
                    m['Filepath'] for m in self.file_records
                    if (not m['Sample_Role'].startswith('CC') and
                        m['Compound'] == pep and m['Cell_Line'] == cell and
                        np.isclose(m['Concentration_uM'], conc, atol=1e-2))
                ]
                ov = {}
                if sample_fps and sample_fps[0] in self.sample_filter_overrides:
                    ov = self.sample_filter_overrides[sample_fps[0]]

                scope_lbl = ov.get('scope_label', str(grp['Applied_Scope'].iloc[0]) if 'Applied_Scope' in grp else 'Global')
                v1l_val = float(ov.get('log_excl_v1l', grp['Applied_V1L_Excl_Log10'].iloc[0] if 'Applied_V1L_Excl_Log10' in grp else self.var_log_excl_v1l.get()))
                v2l_val = float(ov.get('log_cutoff', grp['Applied_V2L_Cutoff_Log10'].iloc[0] if 'Applied_V2L_Cutoff_Log10' in grp else self.var_log_cutoff.get()))
                fsc_min_val = float(ov.get('fsc_min', float(self.ent_fsc_min.get() or 50000)))
                ssc_min_val = float(ov.get('ssc_min', float(self.ent_ssc_min.get() or 10000)))

                match_k = next((k for k in self.curve_thresholds if k[0] == pep and k[1] == cell and np.isclose(k[2], conc, atol=1e-2)), None)
                if match_k:
                    tlag_val = float(self.curve_thresholds[match_k]['t_lag'])
                    tmax_val = float(self.curve_thresholds[match_k]['t_max'])
                    tmode_lbl = "Manual" if self.curve_thresholds[match_k]['mode'] == 'manual' else "Auto"
                else:
                    tlag_val, tmax_val, tmode_lbl = 60.0, 300.0, "Auto"

                v2l_max_pct = 0.0
                if not self.df_kinetics.empty:
                    k_sub = self.df_kinetics[
                        (self.df_kinetics['Peptide'] == pep) &
                        (self.df_kinetics['Cell_Line'] == cell) &
                        (np.isclose(self.df_kinetics['Concentration_uM'], conc, atol=1e-2))
                    ]
                    if not k_sub.empty:
                        v2l_max_pct = float(k_sub['Max_V2L_Pct'].mean())

                self.tree_sample_cutoffs.insert(
                    "", "end", iid=iid,
                    values=(
                        cid_label,
                        scope_lbl,
                        f"{v1l_val:.2f}",
                        f"{v2l_val:.2f}",
                        f"{fsc_min_val:.0f}",
                        f"{ssc_min_val:.0f}",
                        f"{tlag_val:.0f}",
                        f"{tmax_val:.0f}",
                        tmode_lbl,
                        f"{v2l_max_pct:.1f}%"
                    )
                )
                if cid_label == active_curve_str:
                    selected_iid = iid

            if selected_iid and self.tree_sample_cutoffs.exists(selected_iid):
                self.tree_sample_cutoffs.selection_set(selected_iid)
                self.tree_sample_cutoffs.see(selected_iid)
        finally:
            self._updating_cutoff_tree = False

    def on_select_sample_cutoff_table(self, event=None):
        if getattr(self, '_updating_cutoff_tree', False):
            return
        sel = self.tree_sample_cutoffs.selection()
        if not sel:
            return
        iid = sel[0]
        vals = self.tree_sample_cutoffs.item(iid, "values")
        if not vals:
            return
        cid_label = vals[0]
        if cid_label in self.cmb_active_curve['values'] and self.cmb_active_curve.get() != cid_label:
            self.cmb_active_curve.set(cid_label)
            try:
                self.var_log_excl_v1l.set(round(float(vals[2]), 2))
                self.var_log_cutoff.set(round(float(vals[3]), 2))
            except Exception:
                pass
            self.on_active_curve_selected()

    def on_double_click_sample_cutoff_table(self, event):
        """Permite editar en línea cualquier celda de Cut-Off (V1L, V2L, FSC, SSC, T_lag, T_max) con doble clic."""
        region = self.tree_sample_cutoffs.identify("region", event.x, event.y)
        if region != "cell":
            return
        row_id = self.tree_sample_cutoffs.identify_row(event.y)
        col_id = self.tree_sample_cutoffs.identify_column(event.x)
        if not row_id or not col_id:
            return

        col_idx = int(col_id.replace("#", "")) - 1
        # Columnas editables: 2=V1L_Log, 3=V2L_Log, 4=FSC_Min, 5=SSC_Min, 6=T_lag_s, 7=T_max_s
        if col_idx not in (2, 3, 4, 5, 6, 7):
            self.open_selected_sample_cutoff_editor_modal()
            return

        bbox = self.tree_sample_cutoffs.bbox(row_id, col_id)
        if not bbox:
            return
        bx, by, bw, bh = bbox
        cur_vals = list(self.tree_sample_cutoffs.item(row_id, "values"))
        cur_text = str(cur_vals[col_idx])

        entry_ed = tk.Entry(
            self.tree_sample_cutoffs, font=("DejaVu Sans", 8, "bold"),
            bg="#FEF9E7", fg="#002B49", justify="center", relief="solid", bd=1
        )
        entry_ed.place(x=bx, y=by, width=max(52, bw), height=max(20, bh))
        entry_ed.insert(0, cur_text)
        entry_ed.select_range(0, tk.END)
        entry_ed.focus_set()

        committed = {'done': False}

        def commit_edit(ev=None):
            if committed['done']:
                return
            committed['done'] = True
            new_str = entry_ed.get().strip()
            entry_ed.destroy()
            try:
                new_val = float(new_str)
            except ValueError:
                return
            cur_vals[col_idx] = f"{new_val:.2f}" if col_idx in (2, 3) else f"{new_val:.0f}"
            self._apply_single_sample_cutoff_override(
                row_id,
                v1l_log=float(cur_vals[2]),
                v2l_log=float(cur_vals[3]),
                fsc_min=float(cur_vals[4]),
                ssc_min=float(cur_vals[5]),
                tlag_s=float(cur_vals[6]),
                tmax_s=float(cur_vals[7]),
                edited_col_idx=col_idx
            )

        entry_ed.bind("<Return>", commit_edit)
        entry_ed.bind("<FocusOut>", commit_edit)
        entry_ed.bind("<Escape>", lambda e: entry_ed.destroy())

    def _apply_single_sample_cutoff_override(self, row_iid, v1l_log, v2l_log, fsc_min, ssc_min, tlag_s, tmax_s, edited_col_idx=None):
        parts = row_iid.split("||")
        if len(parts) != 3:
            return
        pep, cell, conc = parts[0], parts[1], float(parts[2])

        v1l_log = float(np.clip(v1l_log, 1.50, 5.20))
        v2l_log = float(np.clip(v2l_log, v1l_log + 0.15, 6.20))
        fsc_min = max(0.0, float(fsc_min))
        ssc_min = max(0.0, float(ssc_min))
        tlag_s = float(np.clip(tlag_s, 0.0, 585.0))
        tmax_s = float(np.clip(max(tlag_s + 15.0, tmax_s), 15.0, 600.0))

        # Si se editó T_lag o T_max, actualizar umbral de esa curva a modo manual
        match_k = next((k for k in self.curve_thresholds if k[0] == pep and k[1] == cell and np.isclose(k[2], conc, atol=1e-2)), (pep, cell, conc))
        if edited_col_idx in (6, 7, None):
            self.curve_thresholds[match_k] = {'mode': 'manual', 't_lag': tlag_s, 't_max': tmax_s}

        # Si se editó V1L, V2L, FSC o SSC, actualizar override individual de los tubos de esa muestra
        if edited_col_idx in (2, 3, 4, 5, None):
            bin_size = int(float(self.ent_bin_size.get().strip() or 15))
            fsc_max = self._parse_bound_value(self.ent_fsc_max.get())
            ssc_max = self._parse_bound_value(self.ent_ssc_max.get())
            for meta in self.file_records:
                if (meta['Compound'] == pep and meta['Cell_Line'] == cell and
                        np.isclose(meta['Concentration_uM'], conc, atol=1e-2)):
                    fp = meta['Filepath']
                    self.sample_filter_overrides[fp] = {
                        'bin_size': bin_size,
                        'fsc_min': fsc_min, 'fsc_max': fsc_max,
                        'ssc_min': ssc_min, 'ssc_max': ssc_max,
                        'log_excl_v1l': v1l_log,
                        'log_cutoff': v2l_log,
                        'scope_label': 'Individual'
                    }
            self.var_log_excl_v1l.set(round(v1l_log, 2))
            self.var_log_cutoff.set(round(v2l_log, 2))

        # Reconstruir registros desde la caché en memoria con los overrides de cada muestra
        records = []
        tot_raw, tot_v1l_omitted, tot_ret = 0, 0, 0
        used_channel_name = self.var_perm_channel.get()
        for meta in self.file_records:
            fp = meta['Filepath']
            if fp not in self.tube_event_cache:
                continue
            ov = self.sample_filter_overrides.get(fp, {
                'bin_size': int(float(self.ent_bin_size.get().strip() or 15)),
                'fsc_min': fsc_min, 'fsc_max': np.inf,
                'ssc_min': ssc_min, 'ssc_max': np.inf,
                'log_excl_v1l': v1l_log, 'log_cutoff': v2l_log,
                'scope_label': 'Global'
            })
            try:
                f_max = float(ov['fsc_max']) if ov.get('fsc_max') is not None else np.inf
                s_max = float(ov['ssc_max']) if ov.get('ssc_max') is not None else np.inf
                df_t = apply_gate_and_cutoffs_to_cache(
                    fp, self.tube_event_cache[fp],
                    bin_size=ov['bin_size'],
                    fsc_min=ov['fsc_min'], fsc_max=f_max,
                    ssc_min=ov['ssc_min'], ssc_max=s_max,
                    log_excl_v1l=ov['log_excl_v1l'], log_cutoff=ov['log_cutoff'],
                    scope_label=ov['scope_label']
                )
                df_t['Applied_Bin_s'] = int(ov['bin_size'])
                df_t['Is_CC'] = meta['Sample_Role'].startswith('CC')
                used_channel_name = str(df_t['Perm_Channel_Used'].iloc[0])
                tot_raw += int(df_t['Total_Raw_Events'].iloc[0])
                tot_v1l_omitted += int(df_t['Excluded_V1L_Events'].iloc[0])
                tot_ret += int(df_t['Retained_Events'].iloc[0])
                records.append(df_t)
            except Exception:
                pass

        if records:
            cells_per_ml = float(self.ent_density.get().strip() or 1e5)
            pct_ret = (tot_ret / max(1, tot_raw)) * 100.0
            self.lbl_gate_stats.config(
                text=f"Canal: {used_channel_name} | V1R Analizados: {tot_ret:,} ({pct_ret:.1f}%) | V1L Omitidos: {tot_v1l_omitted:,}"
            )
            self._rebuild_normalized_full_dataframe(records, cells_per_ml)

        self.recompute_kinetics_dataframe()
        self.update_preanalysis_plots()
        self.recompute_all_models_and_tables()

    def open_selected_sample_cutoff_editor_modal(self):
        """Abre un diálogo modal para editar cómodamente todos los Cut-Offs y T_lag/T_max de la muestra seleccionada."""
        sel = self.tree_sample_cutoffs.selection()
        if not sel:
            messagebox.showinfo("Seleccionar Muestra", "Selecciona primero una fila en la Tabla de Cut-Offs por Muestra.", parent=self.root)
            return
        row_id = sel[0]
        vals = self.tree_sample_cutoffs.item(row_id, "values")
        if not vals:
            return

        dlg = tk.Toplevel(self.root)
        dlg.title(f"Editar Cut-Offs de Muestra: {vals[0]}")
        dlg.geometry("420x310")
        dlg.transient(self.root)
        dlg.grab_set()

        frm = ttk.LabelFrame(dlg, text=f" Parámetros Individuales: {vals[0]} ", padding=10)
        frm.pack(fill="both", expand=True, padx=10, pady=8)

        fields = [
            ("Exclusión Detritos V1L (< log10):", vals[2]),
            ("Cut-Off Afectados V2L/V2R (log10):", vals[3]),
            ("Umbral Morfológico FSC-A Mín (>):", vals[4]),
            ("Umbral Morfológico SSC-A Mín (>):", vals[5]),
            ("Inicio Fase Activa T_lag (s):", vals[6]),
            ("Fin Fase Activa / Meseta T_max (s):", vals[7]),
        ]
        entries = []
        for idx_f, (lbl_txt, init_v) in enumerate(fields):
            ttk.Label(frm, text=lbl_txt, font=("DejaVu Sans", 8, "bold")).grid(row=idx_f, column=0, sticky="w", pady=4, padx=4)
            ent = ttk.Entry(frm, width=14)
            ent.insert(0, str(init_v))
            ent.grid(row=idx_f, column=1, sticky="ew", pady=4, padx=4)
            entries.append(ent)

        def save_modal():
            try:
                v1l_v = float(entries[0].get().strip())
                v2l_v = float(entries[1].get().strip())
                fsc_v = float(entries[2].get().strip())
                ssc_v = float(entries[3].get().strip())
                tlag_v = float(entries[4].get().strip())
                tmax_v = float(entries[5].get().strip())
            except ValueError as err:
                messagebox.showerror("Valor Inválido", f"Todos los campos deben ser numéricos:\n{err}", parent=dlg)
                return
            dlg.destroy()
            self._apply_single_sample_cutoff_override(row_id, v1l_v, v2l_v, fsc_v, ssc_v, tlag_v, tmax_v, edited_col_idx=None)

        btn_f = ttk.Frame(frm)
        btn_f.grid(row=len(fields), column=0, columnspan=2, pady=10)
        tk.Button(btn_f, text="✅ Aplicar a Esta Muestra", bg="#117864", fg="white", font=("DejaVu Sans", 8, "bold"), padx=10, command=save_modal).pack(side="left", padx=4)
        ttk.Button(btn_f, text="Cancelar", command=dlg.destroy).pack(side="left", padx=4)

    def export_sample_cutoffs_csv(self):
        if not hasattr(self, 'tree_sample_cutoffs') or not self.tree_sample_cutoffs.get_children():
            messagebox.showwarning("Sin Cut-Offs", "Ejecuta primero el Pre-Análisis para poblar la tabla de Cut-Offs por muestra.", parent=self.root)
            return
        out_dir = Path(self.ent_outdir.get().strip() or DEFAULT_OUTPUT_DIR)
        out_dir.mkdir(parents=True, exist_ok=True)
        default_path = out_dir / "Sample_CutOffs_Table_Panel2.csv"
        fpath = filedialog.asksaveasfilename(
            parent=self.root,
            title="Exportar Tabla de Cut-Offs por Muestra (.csv)",
            initialdir=str(out_dir),
            initialfile=default_path.name,
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not fpath:
            return
        rows = []
        cols = ["Muestra", "Scope", "V1L_Excl_Log10", "V2L_Cutoff_Log10", "FSC_Min", "SSC_Min", "T_lag_s", "T_max_s", "Modo_T", "V2L_Max_Pct"]
        for iid in self.tree_sample_cutoffs.get_children():
            rows.append(dict(zip(cols, self.tree_sample_cutoffs.item(iid, "values"))))
        pd.DataFrame(rows).to_csv(fpath, index=False)
        messagebox.showinfo("Cut-Offs Exportados", f"Tabla de Cut-Offs por muestra exportada exitosamente en:\n{fpath}", parent=self.root)

    def _build_histogram_series_data(self, sel_peps=None, sel_cells=None, sel_doses=None, unique_curves=None, y_mode="pct"):
        """
        Construye dinámicamente las curvas del histograma log10(FITC) sincronizadas con
        las muestras seleccionadas en la Caja 2 ('2. Selección Dinámica de Muestras a Graficar').
        """
        bins_edges = np.linspace(1.5, 6.4, 85)
        hist_series = []

        if sel_peps is None or sel_cells is None or sel_doses is None:
            sel_peps, sel_cells, sel_doses = self._get_selected_sample_filters()

        # 1. Pool de eventos del Control Basal CC (-60 a 0 s) de las muestras seleccionadas
        cc_logs = []
        for meta in self.file_records:
            if (meta['Sample_Role'].startswith('CC') and
                    meta['Compound'] in sel_peps and
                    meta['Cell_Line'] in sel_cells and
                    meta['Filepath'] in self.tube_event_cache):
                ev_c = self.tube_event_cache[meta['Filepath']]
                ov = self.sample_filter_overrides.get(meta['Filepath'], {})
                f_min = ov.get('fsc_min', 50000.0)
                s_min = ov.get('ssc_min', 10000.0)
                m_mask = (ev_c['fsc'] > f_min) & (ev_c['ssc'] > s_min)
                cc_logs.append(ev_c['log10_fitc'][m_mask])

        if cc_logs:
            all_cc = np.concatenate(cc_logs)
            if len(all_cc) > 0:
                counts_cc, _ = np.histogram(all_cc, bins=bins_edges)
                y_cc = (counts_cc / max(1, len(all_cc)) * 100.0) if y_mode == "pct" else counts_cc.astype(float)
                hist_series.append({'label': 'Control Basal CC (-60 a 0s)', 'color': '#117733', 'dashed': True, 'y': y_cc})

        # 2. Distribuciones de cada muestra seleccionada en la Caja 2 (con los mismos colores de las gráficas 1-3)
        if not unique_curves and not self.df_full.empty:
            sub = self.df_full[
                (self.df_full['Peptide'].isin(sel_peps)) &
                (self.df_full['Cell_Line'].isin(sel_cells)) &
                (np.isclose(self.df_full['Concentration_uM'].values[:, None], sel_doses, atol=1e-2).any(axis=1))
            ]
            unique_curves = sorted(
                sub.apply(lambda r: f"{r['Peptide']} | {r['Cell_Line']} | {r['Concentration_uM']:.2f} µM", axis=1).unique()
            )

        control_cell = self.ent_control.get().strip()
        for idx, cid in enumerate(unique_curves or []):
            parts = [p.strip() for p in cid.split("|")]
            if len(parts) != 3:
                continue
            pep, cell, conc_str = parts[0], parts[1], parts[2]
            conc = float(conc_str.replace("µM", "").strip())
            col_hex = CURVE_PALETTE[idx % len(CURVE_PALETTE)]

            curve_logs = []
            sample_v1l = float(self.var_log_excl_v1l.get())
            sample_v2l = float(self.var_log_cutoff.get())
            sample_scope = "Global"
            for meta in self.file_records:
                if (not meta['Sample_Role'].startswith('CC') and
                        meta['Compound'] == pep and meta['Cell_Line'] == cell and
                        np.isclose(meta['Concentration_uM'], conc, atol=1e-2) and
                        meta['Filepath'] in self.tube_event_cache):
                    ev_c = self.tube_event_cache[meta['Filepath']]
                    ov = self.sample_filter_overrides.get(meta['Filepath'], {})
                    f_min = ov.get('fsc_min', 50000.0)
                    s_min = ov.get('ssc_min', 10000.0)
                    sample_v1l = float(ov.get('log_excl_v1l', sample_v1l))
                    sample_v2l = float(ov.get('log_cutoff', sample_v2l))
                    sample_scope = str(ov.get('scope_label', sample_scope))
                    m_mask = (ev_c['fsc'] > f_min) & (ev_c['ssc'] > s_min)
                    curve_logs.append(ev_c['log10_fitc'][m_mask])

            if curve_logs:
                arr_c = np.concatenate(curve_logs)
                if len(arr_c) > 0:
                    c_hist, _ = np.histogram(arr_c, bins=bins_edges)
                    y_c = (c_hist / len(arr_c) * 100.0) if y_mode == "pct" else c_hist.astype(float)
                    hist_series.append({
                        'label': cid,
                        'color': col_hex,
                        'dashed': (cell == control_cell),
                        'y': y_c,
                        'v1l_log': sample_v1l,
                        'v2l_log': sample_v2l,
                        'scope': sample_scope
                    })

        return bins_edges, hist_series

    def _update_histogram_canvas(self, sel_peps=None, sel_cells=None, sel_doses=None, unique_curves=None, v2l_summary_str=""):
        y_mode = "pct" if "Porcentaje" in self.var_hist_ymode.get() else "count"
        bins_edges, hist_series = self._build_histogram_series_data(
            sel_peps, sel_cells, sel_doses, unique_curves, y_mode=y_mode
        )
        scope_short = "Discriminado" if "Discriminado" in self.var_cutoff_scope.get() else "Global"
        self.hist_cutoff_plot.set_histogram_data(
            bins_edges, hist_series,
            cutoff_log=self.var_log_cutoff.get(),
            excl_log=self.var_log_excl_v1l.get(),
            y_mode=y_mode,
            v2l_summary_text=v2l_summary_str,
            scope_text=scope_short
        )

    # --------------------------------------------------------------------------
    # MOTOR DE ANÁLISIS CATEGÓRICO POR VENTANA DE TIEMPO, PRUEBAS E INSIGHTS
    # --------------------------------------------------------------------------
    def _compute_categorical_time_windows_and_tests(self, fdr_method, target_cell, control_cell):
        """
        Construye las 3 tablas maestras del Panel de Análisis Categórico:
          1) self.df_categorical_summary: Hitos de efecto por compuesto y concentración
             (% Afectados Mín/Máx, % No Afectados Mín, T_efecto_mín, T_efecto_máx, T50 del efecto,
              LT50 poblacional, Span T10-T90, T_onset estadístico, V_cat, Entropía, BC Sarle,
              h de Cohen, OR, p-Fisher, q-FDR y Decisión de Afectación).
          2) self.df_cat_time_windows: Proporciones entre eventos afectados (V2L) y no afectados (V2R)
             para cada compuesto, concentración y ventana de tiempo t (-60 a 600 s).
          3) self.df_cat_hypothesis_tests: Pruebas estadísticas de tendencia dosis-respuesta de
             Cochran-Armitage, Log-Rank de supervivencia celular, Isocronas EC50_V2L(t) y contrastes.
        """
        if self.df_full.empty or self.df_kinetics.empty:
            return

        # 1. Tabla Resumen de Hitos Categóricos por Compuesto y Concentración (integrando Cut-Offs configurados)
        cat_sum = self.df_kinetics.groupby(['Peptide', 'Cell_Line', 'Concentration_uM']).agg(
            Filter_Scope=('Filter_Scope', 'first'),
            V1L_Excl_Log10=('V1L_Excl_Log10', 'first'),
            V2L_Cutoff_Log10=('V2L_Cutoff_Log10', 'first'),
            T_lag_s=('T_lag', 'mean'),
            T_max_s=('T_max', 'mean'),
            Tlag_Mode=('Tlag_Mode', 'first'),
            Basal_V2L_Pct=('Basal_V2L_Pct', 'mean'),
            Min_V2L_Pct=('Min_V2L_Pct', 'mean'),
            Max_V2L_Pct=('Max_V2L_Pct', 'mean'),
            Min_Unaffected_V2R_Pct=('Min_Unaffected_V2R_Pct', 'mean'),
            End_V2L_Pct=('End_V2L_Pct', 'mean'),
            T_Effect_Min_s=('T_Effect_Min_s', 'mean'),
            T_Effect_Max_s=('T_Effect_Max_s', 'mean'),
            T50_Effect_s=('T50_Effect_s', 'mean'),
            LT50_s=('LT50_s', 'mean'),
            Span_T10_T90_s=('Span_T10_T90_s', 'mean'),
            T_Stat_Onset_s=('T_Stat_Onset_s', 'mean'),
            T_Peak_Conversion_s=('T_Peak_Conversion_s', 'mean'),
            V_Cat_Max_Inst=('V_Cat_Max_Inst', 'mean'),
            V_Cat_Active=('V_Cat_Active', 'mean'),
            AUC_V2L_Pct=('AUC_V2L_Pct', 'mean'),
            Max_Entropy_Bits=('Max_Entropy_Bits', 'mean'),
            Hazard_Rate_s=('Hazard_Rate_s', 'mean'),
            Sync_Index=('Sync_Index', 'mean'),
            Bimodality_BC=('Bimodality_BC', 'mean'),
            Cohen_h_Max=('Cohen_h_Max', 'mean'),
            OR_vs_Basal=('OR_vs_Basal', 'mean'),
            p_Fisher_Max=('p_Fisher_Max', 'min'),
            Mode_Delta=('Mode_Delta', 'mean'),
            Mechanism_Class=('Mechanism_Class', 'first')
        ).reset_index()

        if not cat_sum.empty:
            q_cat = adjust_pvalues(cat_sum['p_Fisher_Max'].values, fdr_method)
            cat_sum['q_FDR_Basal'] = q_cat
            cat_sum['Affectation_State'] = [
                classify_affectation_status(qv, hv, mv)
                for qv, hv, mv in zip(cat_sum['q_FDR_Basal'], cat_sum['Cohen_h_Max'], cat_sum['Max_V2L_Pct'])
            ]
        self.df_categorical_summary = cat_sum

        # 2. Tabla de Proporciones Afectados vs No Afectados por Ventana de Tiempo (con Cut-Offs aplicados)
        df_f = self.df_full.copy()
        if 'Affected_V2L_Events' not in df_f.columns:
            df_f['Affected_V2L_Events'] = np.round(df_f['Events'] * df_f['Affected_V2L_Pct'] / 100.0).astype(int)
        if 'Unaffected_V2R_Events' not in df_f.columns:
            df_f['Unaffected_V2R_Events'] = (df_f['Events'] - df_f['Affected_V2L_Events']).astype(int)
        if 'Applied_Scope' not in df_f.columns:
            df_f['Applied_Scope'] = "Global"
        if 'Applied_V1L_Excl_Log10' not in df_f.columns:
            df_f['Applied_V1L_Excl_Log10'] = float(self.var_log_excl_v1l.get())
        if 'Applied_V2L_Cutoff_Log10' not in df_f.columns:
            df_f['Applied_V2L_Cutoff_Log10'] = float(self.var_log_cutoff.get())

        win_grp = df_f.groupby(['Peptide', 'Cell_Line', 'Concentration_uM', 'Time_Window']).agg(
            Filter_Scope=('Applied_Scope', 'first'),
            V1L_Excl_Log10=('Applied_V1L_Excl_Log10', 'first'),
            V2L_Cutoff_Log10=('Applied_V2L_Cutoff_Log10', 'first'),
            N_Rep=('Replicate', 'nunique'),
            Total_V1R_Events=('Events', 'sum'),
            Affected_V2L_Events=('Affected_V2L_Events', 'sum'),
            Unaffected_V2R_Events=('Unaffected_V2R_Events', 'sum'),
            Mean_Affected_V2L_Pct=('Affected_V2L_Pct', 'mean'),
            SEM_V2L_Pct=('Affected_V2L_Pct', 'sem'),
            Mean_Unaffected_V2R_Pct=('Unaffected_V2R_Pct', 'mean')
        ).reset_index()
        win_grp['SEM_V2L_Pct'] = win_grp['SEM_V2L_Pct'].fillna(0.0)

        win_rows = []
        for (pep, cell, conc), c_grp in win_grp.groupby(['Peptide', 'Cell_Line', 'Concentration_uM']):
            c_sorted = c_grp.sort_values('Time_Window').copy()
            base_sub = c_sorted[c_sorted['Time_Window'] < 0]
            if not base_sub.empty:
                b_aff = int(base_sub['Affected_V2L_Events'].sum())
                b_unaff = int(base_sub['Unaffected_V2R_Events'].sum())
                b_pct = float(base_sub['Mean_Affected_V2L_Pct'].mean())
            else:
                first_r = c_sorted.iloc[0]
                b_aff = int(first_r['Affected_V2L_Events'])
                b_unaff = int(first_r['Unaffected_V2R_Events'])
                b_pct = float(first_r['Mean_Affected_V2L_Pct'])

            t_vec = c_sorted['Time_Window'].values.astype(float)
            v2l_vec = c_sorted['Mean_Affected_V2L_Pct'].values.astype(float)
            dt_vec = np.gradient(t_vec) if len(t_vec) >= 2 else np.array([15.0])
            dt_vec[dt_vec == 0] = 15.0
            d_v2l_dt = np.gradient(v2l_vec) / dt_vec if len(v2l_vec) >= 2 else np.zeros_like(v2l_vec)

            raw_ps = []
            temp_list = []
            for idx_w, (_, r_w) in enumerate(c_sorted.iterrows()):
                tw = float(r_w['Time_Window'])
                tot_ev = int(r_w['Total_V1R_Events'])
                aff_ev = int(r_w['Affected_V2L_Events'])
                unaff_ev = int(r_w['Unaffected_V2R_Events'])
                m_aff = float(r_w['Mean_Affected_V2L_Pct'])
                se_aff = float(r_w['SEM_V2L_Pct'])
                m_unaff = float(r_w['Mean_Unaffected_V2R_Pct'])

                ci_lo, ci_hi = wilson_binomial_ci(aff_ev, tot_ev)
                odds_t = float((aff_ev + 0.5) / (unaff_ev + 0.5))
                log2_odds = float(math.log2(max(1e-6, odds_t)))
                p_prob = float(np.clip(m_aff / 100.0, 1e-6, 1.0 - 1e-6))
                shannon_w = float(-(p_prob * math.log2(p_prob) + (1.0 - p_prob) * math.log2(1.0 - p_prob))) if 1.0 <= m_aff <= 99.0 else 0.0

                chi2_w, pv_w, or_b, _, _, _ = contingency_2x2_p_and_or(aff_ev, unaff_ev, b_aff, b_unaff)
                h_w = cohen_h_proportions(m_aff, b_pct)
                raw_ps.append(pv_w)
                temp_list.append({
                    'Peptide': pep, 'Cell_Line': cell, 'Concentration_uM': conc,
                    'Filter_Scope': str(r_w['Filter_Scope']),
                    'V1L_Excl_Log10': round(float(r_w['V1L_Excl_Log10']), 2),
                    'V2L_Cutoff_Log10': round(float(r_w['V2L_Cutoff_Log10']), 2),
                    'Time_Window': tw, 'N_Rep': int(r_w['N_Rep']),
                    'Total_V1R_Events': tot_ev,
                    'Affected_V2L_Events': aff_ev,
                    'Unaffected_V2R_Events': unaff_ev,
                    'Mean_Affected_V2L_Pct': m_aff,
                    'SEM_V2L_Pct': se_aff,
                    'Mean_Unaffected_V2R_Pct': m_unaff,
                    'Wilson_95CI_V2L': f"[{ci_lo:.1f}%, {ci_hi:.1f}%]",
                    'Ratio_V2L_to_V2R': odds_t,
                    'Log2_Odds': log2_odds,
                    'Inst_Conversion_Rate_pct_s': float(d_v2l_dt[idx_w]),
                    'Shannon_H_Bits': shannon_w,
                    'Cohen_h_vs_Basal': h_w,
                    'OR_vs_Basal': or_b,
                    'Chi2_Stat': chi2_w,
                    'p_Fisher': pv_w
                })

            q_ws = adjust_pvalues(raw_ps, fdr_method)
            for item, qw in zip(temp_list, q_ws):
                item['q_FDR'] = float(qw)
                if item['Time_Window'] < 0:
                    item['Window_Status'] = "Control Basal CC"
                else:
                    item['Window_Status'] = classify_affectation_status(qw, item['Cohen_h_vs_Basal'], item['Mean_Affected_V2L_Pct'])
                win_rows.append(item)

        self.df_cat_time_windows = pd.DataFrame(win_rows)

        # 3. Tabla de Pruebas Estadísticas de Afectación (Tendencia Cochran-Armitage, Log-Rank, Isocronas EC50(t) y Contrastes)
        test_rows = []
        if not self.df_cat_time_windows.empty:
            # 3A. Prueba de Tendencia Dosis-Respuesta de Cochran-Armitage e Isocronas EC50_V2L(t) en tiempos clave
            key_times = [60.0, 120.0, 180.0, 300.0, 450.0, 600.0]
            for (pep, cell), grp_pc in self.df_cat_time_windows.groupby(['Peptide', 'Cell_Line']):
                avail_times = sorted(grp_pc[grp_pc['Time_Window'] >= 0]['Time_Window'].unique())
                if not avail_times:
                    continue
                for kt in key_times:
                    nearest_t = float(avail_times[int(np.argmin(np.abs(np.asarray(avail_times) - kt)))])
                    sub_t = grp_pc[np.isclose(grp_pc['Time_Window'], nearest_t, atol=1.0)].sort_values('Concentration_uM')
                    if len(sub_t) < 2:
                        continue
                    d_arr = sub_t['Concentration_uM'].values
                    aff_arr = sub_t['Affected_V2L_Events'].values
                    tot_arr = sub_t['Total_V1R_Events'].values
                    pct_arr = sub_t['Mean_Affected_V2L_Pct'].values

                    z_ca, p_ca, slope_ca = cochran_armitage_trend_test(d_arr, aff_arr, tot_arr)
                    # Estimación de Isocrona EC50_V2L(t) (concentración para lograr 50% de eventos afectados en el tiempo t)
                    if np.max(pct_arr) >= 50.0:
                        ec50_t = interpolate_crossing_time(d_arr, pct_arr, 50.0, start_time=0.0)
                        ec50_str = f"{ec50_t:.2f} µM"
                    else:
                        ec50_str = f"> {np.max(d_arr):.2f} µM"

                    test_rows.append({
                        'Analysis_Type': 'Tendencia Cochran-Armitage + Isocrona EC50(t)',
                        'Peptide': pep,
                        'Comparison_Scope': f"{cell} @ t={nearest_t:.0f}s",
                        'Metric_Focus': f"Max V2L={np.max(pct_arr):.1f}% | EC50(t)={ec50_str}",
                        'Stat_Name': f"Z_trend={z_ca:.2f} (Pend={slope_ca:+.1f}%/µM)",
                        'Effect_Size': f"EC50_V2L={ec50_str}",
                        'raw_p': p_ca,
                        'Decision': "Tendencia Dosis-Respuesta Significativa" if p_ca < 0.05 and z_ca > 0 else "Sin Tendencia Monotónica (NS)"
                    })

            # 3B. Contraste de Proporciones en el Efecto Máximo: Línea Diana vs Línea Control Sano y vs Basal
            for _, r_sum in self.df_categorical_summary.iterrows():
                pep = r_sum['Peptide']
                cell = r_sum['Cell_Line']
                conc = float(r_sum['Concentration_uM'])
                t_max_e = float(r_sum['T_Effect_Max_s'])
                mv2l = float(r_sum['Max_V2L_Pct'])
                t50_e = float(r_sum['T50_Effect_s'])
                h_val = float(r_sum['Cohen_h_Max'])
                or_val = float(r_sum['OR_vs_Basal'])
                pv_b = float(r_sum['p_Fisher_Max'])

                test_rows.append({
                    'Analysis_Type': 'Prueba Exacta/Chi² en Efecto Máximo (vs Basal)',
                    'Peptide': pep,
                    'Comparison_Scope': f"{cell} | {conc:.2f} µM (T_máx={t_max_e:.0f}s, T50={t50_e:.0f}s)",
                    'Metric_Focus': f"Afectados Máx={mv2l:.1f}% vs Basal={r_sum['Basal_V2L_Pct']:.1f}%",
                    'Stat_Name': f"OR={or_val:.1f}x",
                    'Effect_Size': f"h Cohen={h_val:.2f}",
                    'raw_p': pv_b,
                    'Decision': str(r_sum['Affectation_State'])
                })

        df_tests = pd.DataFrame(test_rows)
        if not df_tests.empty:
            q_t = adjust_pvalues(df_tests['raw_p'].values, fdr_method)
            df_tests['p_value'] = [f"{p:.4f}" if p >= 0.0001 else "< 0.0001" for p in df_tests['raw_p']]
            df_tests['q_FDR'] = [f"{q:.4f}" if q >= 0.0001 else "< 0.0001" for q in q_t]
            df_tests.drop(columns=['raw_p'], inplace=True)
        self.df_cat_hypothesis_tests = df_tests

        self._generate_categorical_biophysical_insights_text(target_cell, control_cell)
        if hasattr(self, 'cat_prop_canvas'):
            self.update_categorical_panel_view()

    def _generate_categorical_biophysical_insights_text(self, target_cell, control_cell):
        """
        Genera un informe interpretativo automático y sugerencias metodológicas en la Sub-pestaña 4
        del Panel Categórico, integrando también sus conclusiones en el Reporte Final.
        """
        if not hasattr(self, 'txt_cat_insights'):
            return
        self.txt_cat_insights.delete("1.0", tk.END)
        if self.df_categorical_summary.empty:
            return

        lines = [
            "==========================================================================================",
            " SÍNTESIS BIOFÍSICA DEL PANEL DE PROPORCIONES CATEGÓRICAS (AFECTADOS V2L vs NO AFECTADOS V2R)",
            "==========================================================================================",
            ""
        ]

        t_sub = self.df_categorical_summary[self.df_categorical_summary['Cell_Line'] == target_cell].copy()
        if t_sub.empty:
            t_sub = self.df_categorical_summary.copy()

        for pep, grp_p in t_sub.groupby('Peptide'):
            grp_s = grp_p.sort_values('Concentration_uM')
            best_row = grp_s.loc[grp_s['Max_V2L_Pct'].idxmax()]
            fastest_t50_row = grp_s.loc[grp_s['T50_Effect_s'].idxmin()]
            lines.append(f"▶ COMPUESTO: {pep} (Línea Celular: {best_row['Cell_Line']})")
            lines.append(
                f"  • Rango de % Eventos Afectados (V2L): desde un mínimo post-inyección de {grp_s['Min_V2L_Pct'].min():.1f}% "
                f"hasta un máximo de {best_row['Max_V2L_Pct']:.1f}% (remanente no afectado V2R = {best_row['Min_Unaffected_V2R_Pct']:.1f}%) "
                f"a {best_row['Concentration_uM']:.2f} µM."
            )
            lines.append(
                f"  • Tiempo para lograr el Efecto Máximo (T_efecto_máx): {best_row['T_Effect_Max_s']:.1f} s "
                f"(a {best_row['Concentration_uM']:.2f} µM), con inicio estadístico de afectación en T_onset = {best_row['T_Stat_Onset_s']:.1f} s."
            )
            lines.append(
                f"  • T50 del Efecto (tiempo al 50% del efecto máximo propio de cada dosis) vs LT50 (50% poblacional absoluto): "
                f"El menor T50 del efecto es {fastest_t50_row['T50_Effect_s']:.1f} s ({fastest_t50_row['Concentration_uM']:.2f} µM), "
                f"mientras que el LT50 poblacional es {fastest_t50_row['LT50_s']:.1f} s."
            )
            lines.append(
                f"  • Sincronía y Mecanismo de Acción: Ventana de transición ΔT(10%→90%) = {best_row['Span_T10_T90_s']:.1f} s, "
                f"Índice de Bimodalidad de Sarle BC = {best_row['Bimodality_BC']:.3f}, "
                f"Entropía máxima H = {best_row['Max_Entropy_Bits']:.2f} bits ➔ Clasificación: {best_row['Mechanism_Class']}."
            )
            lines.append(
                f"  • Pruebas Estadísticas de Definición de Afectación: h de Cohen = {best_row['Cohen_h_Max']:.2f}, "
                f"Odds Ratio vs Basal = {best_row['OR_vs_Basal']:.1f}x (q-FDR = {best_row['q_FDR_Basal']:.4f}) ➔ {best_row['Affectation_State']}."
            )
            lines.append("")

        lines.extend([
            "------------------------------------------------------------------------------------------",
            " ANÁLISIS ADICIONALES INTEGRADOS PARA COMPRENDER MEJOR EL FENÓMENO BIOFÍSICO",
            "------------------------------------------------------------------------------------------",
            "1. Distinción entre T50 del Efecto (ET50 relativo) y LT50 Poblacional (Kaplan-Meier):",
            "   - A concentraciones sub-líticas donde el máximo de eventos afectados no supera el 50% (ej. 15-35%),",
            "     el LT50 absoluto permanece en >600 s, pero el T50 del Efecto mide con exactitud qué tan rápido",
            "     esa dosis alcanza la mitad de su capacidad de permeabilización.",
            "2. Amplitud de Transición Temporal ΔT(10%→90%) e Índice de Sincronía (T50 / ΔT10-90):",
            "   - Un ΔT(10%→90%) corto (<90 s) con velocidad instantánea d(%V2L)/dt alta indica un colapso",
            "     cooperativo en avalancha tras superar la tensión crítica de membrana.",
            "3. Coeficiente de Bimodalidad de Sarle (BC > 0.555) y Δ_Modo (AUC_Continua - AUC_V2L):",
            "   - Permite discriminar si las células en el cuadrante No Afectado (V2R) conservan intacta su",
            "     fluorescencia (mecanismo Todo-o-Nada por poro toroidal/barril) o sufren fuga parcial gradual.",
            "4. Isocronas Dosis-Respuesta EC50_V2L(t) y Prueba de Tendencia de Cochran-Armitage:",
            "   - Cuantifica cómo disminuye la concentración necesaria para afectar al 50% de la población",
            "     conforme progresan las ventanas de tiempo (60 s, 120 s, 180 s, 300 s, 450 s, 600 s)."
        ])
        self.txt_cat_insights.insert(tk.END, "\n".join(lines))

    # --------------------------------------------------------------------------
    # MOTOR DE MODELADO CINÉTICO MÚLTIPLE Y ESTADÍSTICA (PESTAÑA 4)
    # --------------------------------------------------------------------------
    def recompute_all_models_and_tables(self):
        if self.df_kinetics.empty:
            return

        target_cell = self.ent_target.get().strip()
        control_cell = self.ent_control.get().strip()
        fdr_method = self.var_fdr_method.get()
        stat_mode = self.var_stat_mode.get()
        analyzer_name = self.var_kinetic_analyzer.get()
        try:
            n_boot = max(20, min(1000, int(self.ent_bootstrap_n.get().strip())))
        except ValueError:
            n_boot = 200

        # Tabla 0: Descripción del Prerregistro y Parámetros Específicos de Pre-Análisis Usados por Muestra
        if 'Applied_Bin_s' not in self.df_kinetics.columns:
            self.df_kinetics['Applied_Bin_s'] = int(float(self.ent_bin_size.get().strip() or 15))
        if 'Applied_FSC_Min' not in self.df_kinetics.columns:
            self.df_kinetics['Applied_FSC_Min'] = float(self.ent_fsc_min.get().strip() or 50000.0)
        if 'Applied_FSC_Max' not in self.df_kinetics.columns:
            self.df_kinetics['Applied_FSC_Max'] = self._parse_bound_value(self.ent_fsc_max.get())
        if 'Applied_SSC_Min' not in self.df_kinetics.columns:
            self.df_kinetics['Applied_SSC_Min'] = float(self.ent_ssc_min.get().strip() or 10000.0)
        if 'Applied_SSC_Max' not in self.df_kinetics.columns:
            self.df_kinetics['Applied_SSC_Max'] = self._parse_bound_value(self.ent_ssc_max.get())

        prereg_grp = self.df_kinetics.groupby(['Peptide', 'Cell_Line', 'Concentration_uM']).agg(
            Filter_Scope=('Filter_Scope', 'first'),
            V1L_Excl_Log10=('V1L_Excl_Log10', 'first'),
            V2L_Cutoff_Log10=('V2L_Cutoff_Log10', 'first'),
            Bin_s=('Applied_Bin_s', 'first'),
            FSC_Min=('Applied_FSC_Min', 'first'),
            FSC_Max=('Applied_FSC_Max', 'first'),
            SSC_Min=('Applied_SSC_Min', 'first'),
            SSC_Max=('Applied_SSC_Max', 'first'),
            Tlag_Mode=('Tlag_Mode', 'first'),
            Exclusion_s=('Exclusion_s', 'first'),
            Mean_Tlag=('T_lag', 'mean'),
            Mean_Tmax=('T_max', 'mean'),
            Max_V2L_Pct=('Max_V2L_Pct', 'mean'),
            Max_Perm_Pct=('Max_Perm', 'mean'),
            Raw_Events=('Raw_Events', 'sum'),
            V1L_Omitted=('V1L_Omitted', 'sum'),
            V1R_Analyzed=('V1R_Analyzed', 'sum')
        ).reset_index()
        prereg_grp['Channel'] = self.var_perm_channel.get()
        prereg_grp['FSC_Gate'] = [
            f">{fmin:.0f}..{fmax:.0f}" if np.isfinite(fmax) else f">{fmin:.0f}..Sin límite"
            for fmin, fmax in zip(prereg_grp['FSC_Min'], prereg_grp['FSC_Max'])
        ]
        prereg_grp['SSC_Gate'] = [
            f">{smin:.0f}..{smax:.0f}" if np.isfinite(smax) else f">{smin:.0f}..Sin límite"
            for smin, smax in zip(prereg_grp['SSC_Min'], prereg_grp['SSC_Max'])
        ]
        prereg_grp['T_lag_s'] = prereg_grp['Mean_Tlag'].round(1)
        prereg_grp['T_max_s'] = prereg_grp['Mean_Tmax'].round(1)
        prereg_grp['Retained_Pct'] = (prereg_grp['V1R_Analyzed'] / np.maximum(1, prereg_grp['Raw_Events'])) * 100.0
        self.df_preregistration = prereg_grp[[
            'Peptide', 'Cell_Line', 'Concentration_uM', 'Channel', 'Bin_s',
            'FSC_Gate', 'SSC_Gate', 'Filter_Scope', 'V1L_Excl_Log10', 'V2L_Cutoff_Log10',
            'T_lag_s', 'T_max_s', 'Tlag_Mode', 'Exclusion_s', 'Raw_Events', 'V1L_Omitted', 'V1R_Analyzed', 'Retained_Pct'
        ]]

        # Tabla 0B: Tabla de Cut-Offs Configurados por el Usuario por Muestra (Panel 2)
        self.df_user_cutoffs_table = prereg_grp.apply(lambda r: pd.Series({
            'Sample_ID': f"{r['Peptide']} | {r['Cell_Line']} | {r['Concentration_uM']:.2f} µM",
            'Peptide': r['Peptide'],
            'Cell_Line': r['Cell_Line'],
            'Concentration_uM': round(float(r['Concentration_uM']), 2),
            'Filter_Scope': r['Filter_Scope'],
            'V1L_Excl_Log10': round(float(r['V1L_Excl_Log10']), 2),
            'V2L_Cutoff_Log10': round(float(r['V2L_Cutoff_Log10']), 2),
            'FSC_Min': round(float(r['FSC_Min']), 0),
            'FSC_Max': round(float(r['FSC_Max']), 0) if np.isfinite(r['FSC_Max']) else "Sin límite",
            'SSC_Min': round(float(r['SSC_Min']), 0),
            'SSC_Max': round(float(r['SSC_Max']), 0) if np.isfinite(r['SSC_Max']) else "Sin límite",
            'Bin_s': int(r['Bin_s']),
            'T_lag_s': round(float(r['Mean_Tlag']), 1),
            'T_max_s': round(float(r['Mean_Tmax']), 1),
            'Tlag_Mode': r['Tlag_Mode'],
            'Max_V2L_Pct': round(float(r['Max_V2L_Pct']), 2),
            'Max_Perm_Pct': round(float(r['Max_Perm_Pct']), 2)
        }), axis=1) if not prereg_grp.empty else pd.DataFrame()

        self.df_v0_summary = self.df_kinetics.groupby(['Peptide', 'Cell_Line', 'Concentration_uM', 'Molecules_Per_Cell']).agg(
            Filter_Scope=('Filter_Scope', 'first'),
            V1L_Excl_Log10=('V1L_Excl_Log10', 'first'),
            V2L_Cutoff_Log10=('V2L_Cutoff_Log10', 'first'),
            Mean_V0_Early=('V0_Early', 'mean'), SEM_V0_Early=('V0_Early', 'sem')
        ).reset_index()

        self.df_v0_600_summary = self.df_kinetics.groupby(['Peptide', 'Cell_Line', 'Concentration_uM', 'Molecules_Per_Cell']).agg(
            Filter_Scope=('Filter_Scope', 'first'),
            V1L_Excl_Log10=('V1L_Excl_Log10', 'first'),
            V2L_Cutoff_Log10=('V2L_Cutoff_Log10', 'first'),
            Mean_V0_600=('V0_600', 'mean'), SEM_V0_600=('V0_600', 'sem')
        ).reset_index()

        self.df_v0_active_summary = self.df_kinetics.groupby(['Peptide', 'Cell_Line', 'Concentration_uM', 'Molecules_Per_Cell']).agg(
            Filter_Scope=('Filter_Scope', 'first'),
            V1L_Excl_Log10=('V1L_Excl_Log10', 'first'),
            V2L_Cutoff_Log10=('V2L_Cutoff_Log10', 'first'),
            T_lag_s=('T_lag', 'mean'),
            T_max_s=('T_max', 'mean'),
            Tlag_Mode=('Tlag_Mode', 'first'),
            Mean_V0_Active=('V0_Active', 'mean'), SEM_V0_Active=('V0_Active', 'sem')
        ).reset_index()

        # Calcular Tablas del Panel Categórico (Hitos T50/T_max, Proporciones por Ventana y Pruebas Estadísticas)
        self._compute_categorical_time_windows_and_tests(fdr_method, target_cell, control_cell)

        # Resumen Morfológico (FSC-A y SSC-A)
        self.df_morphology_summary = self.df_kinetics.groupby(['Peptide', 'Cell_Line', 'Concentration_uM']).agg(
            End_Norm_FSC_Pct=('End_Norm_FSC_Pct', 'mean'),
            SEM_FSC=('End_Norm_FSC_Pct', 'sem'),
            End_Norm_SSC_Pct=('End_Norm_SSC_Pct', 'mean'),
            SEM_SSC=('End_Norm_SSC_Pct', 'sem')
        ).reset_index()

        pbmc_contrasts = []
        ctrl_pool_df = self.df_kinetics[self.df_kinetics['Cell_Line'] == control_cell]
        target_df = self.df_kinetics[self.df_kinetics['Cell_Line'] == target_cell]

        raw_pvals_t4 = []
        for pep, pep_grp in target_df.groupby('Peptide'):
            ctrl_sub = ctrl_pool_df[ctrl_pool_df['Peptide'] == pep]
            for conc, grp in pep_grp.groupby('Concentration_uM'):
                thp_auc = grp['Normalized_AUC_Pct'].values
                p_thp_v2l = np.clip(grp['End_V2L_Pct'].mean() / 100.0, 0.005, 0.995)
                if not ctrl_sub.empty:
                    pbmc_auc = ctrl_sub['Normalized_AUC_Pct'].values
                    p_ctrl_v2l = np.clip(ctrl_sub['End_V2L_Pct'].mean() / 100.0, 0.005, 0.995)
                    diff_mean = float(np.mean(thp_auc) - np.mean(pbmc_auc))
                    ttest_res = stats.ttest_ind(thp_auc, pbmc_auc, equal_var=False)
                    pval = float(ttest_res.pvalue) if not np.isnan(ttest_res.pvalue) else 1.0
                    sem_diff = float(np.sqrt(np.var(thp_auc, ddof=1)/max(1, len(thp_auc)) + np.var(pbmc_auc, ddof=1)/max(1, len(pbmc_auc))))
                    df_welch = max(1, len(thp_auc) + len(pbmc_auc) - 2)
                    ci_margin = stats.t.ppf(0.975, df=df_welch) * sem_diff
                    pbmc_mean = float(np.mean(pbmc_auc))
                    odds_ratio = (p_thp_v2l / (1.0 - p_thp_v2l)) / (p_ctrl_v2l / (1.0 - p_ctrl_v2l))
                else:
                    diff_mean = float(np.mean(thp_auc))
                    pval = 1.0
                    ci_margin = 0.0
                    pbmc_mean = 0.0
                    odds_ratio = (p_thp_v2l / (1.0 - p_thp_v2l)) / (0.01 / 0.99)

                raw_pvals_t4.append(pval)
                pbmc_contrasts.append({
                    'Peptide': pep, 'Dose_uM': conc,
                    'THP1_Tlag': grp['T_lag'].mean(),
                    'THP1_Tmax': grp['T_max'].mean(),
                    'THP1_Vact': grp['V0_Active'].mean(),
                    'THP1_AUC': float(np.mean(thp_auc)),
                    'PBMC_AUC': pbmc_mean,
                    'Delta_AUC': diff_mean,
                    'Odds_Ratio_V2L': odds_ratio,
                    'CI_95': f"[{diff_mean - ci_margin:.2f}, {diff_mean + ci_margin:.2f}]",
                    'raw_p': pval
                })

        self.df_table4 = pd.DataFrame(pbmc_contrasts)
        if not self.df_table4.empty:
            q_vals = adjust_pvalues(self.df_table4['raw_p'].values, fdr_method)
            self.df_table4['p_value'] = [f"{p:.4f}" if p >= 0.0001 else "< 0.0001" for p in self.df_table4['raw_p']]
            self.df_table4['q_Adj'] = [f"{q:.4f}" if q >= 0.0001 else "< 0.0001" for q in q_vals]

        t5_records = []
        raw_p_t5 = []
        if "Entre Compuestos" in stat_mode:
            peps = sorted(target_df['Peptide'].unique())
            for i in range(len(peps)):
                for j in range(i + 1, len(peps)):
                    p1, p2 = peps[i], peps[j]
                    g1 = target_df[target_df['Peptide'] == p1].sort_values('Concentration_uM')
                    g2 = target_df[target_df['Peptide'] == p2].sort_values('Concentration_uM')
                    v1 = g1['Normalized_AUC_Pct'].values
                    v2 = g2['Normalized_AUC_Pct'].values
                    if len(v1) > 1 and len(v2) > 1:
                        tt = stats.ttest_ind(v1, v2, equal_var=False)
                        pv = float(tt.pvalue) if not np.isnan(tt.pvalue) else 1.0
                    else:
                        pv = 1.0
                    diff = float(np.mean(v1) - np.mean(v2))
                    raw_p_t5.append(pv)
                    t5_records.append({
                        'Peptide': f"{p1} vs {p2}", 'Dose_uM': 0.0,
                        'T_lag': g1['T_lag'].mean(), 'T_max': g1['T_max'].mean(),
                        'V_active': g1['V0_Active'].mean() - g2['V0_Active'].mean(),
                        'Mean_AUC': float(np.mean(v1)), 'Delta_vs_Ref': diff,
                        'CI_95': f"[{diff - 5.0:.2f}, {diff + 5.0:.2f}]", 'raw_p': pv
                    })
        elif "Entre Dosis Sucesivas" in stat_mode:
            for pep, grp_pep in target_df.groupby('Peptide'):
                doses = sorted(grp_pep['Concentration_uM'].unique())
                for k in range(1, len(doses)):
                    d_prev, d_curr = doses[k - 1], doses[k]
                    g_prev = grp_pep[grp_pep['Concentration_uM'] == d_prev]
                    g_curr = grp_pep[grp_pep['Concentration_uM'] == d_curr]
                    v_prev = g_prev['Normalized_AUC_Pct'].values
                    v_curr = g_curr['Normalized_AUC_Pct'].values
                    tt = stats.ttest_ind(v_curr, v_prev, equal_var=False) if len(v_curr) > 1 and len(v_prev) > 1 else None
                    pv = float(tt.pvalue) if tt and not np.isnan(tt.pvalue) else 1.0
                    diff = float(np.mean(v_curr) - np.mean(v_prev))
                    raw_p_t5.append(pv)
                    t5_records.append({
                        'Peptide': pep, 'Dose_uM': d_curr,
                        'T_lag': g_curr['T_lag'].mean(), 'T_max': g_curr['T_max'].mean(),
                        'V_active': g_curr['V0_Active'].mean(),
                        'Mean_AUC': float(np.mean(v_curr)), 'Delta_vs_Ref': diff,
                        'CI_95': f"(vs {d_prev:.2f} µM)", 'raw_p': pv
                    })
        else:
            for pep, grp_pep in target_df.groupby('Peptide'):
                for conc, grp_dose in grp_pep.groupby('Concentration_uM'):
                    thp_auc = grp_dose['Normalized_AUC_Pct'].values
                    diff_mean = float(np.mean(thp_auc))
                    if len(thp_auc) > 1 and np.std(thp_auc, ddof=1) > 0:
                        ttest_res = stats.ttest_1samp(thp_auc, 0.0)
                        pval = float(ttest_res.pvalue)
                        sem_val = float(stats.sem(thp_auc))
                        ci_margin = stats.t.ppf(0.975, df=len(thp_auc) - 1) * sem_val
                    else:
                        pval = 1.0
                        ci_margin = 0.0
                    raw_p_t5.append(pval)
                    t5_records.append({
                        'Peptide': pep, 'Dose_uM': conc,
                        'T_lag': grp_dose['T_lag'].mean(), 'T_max': grp_dose['T_max'].mean(),
                        'V_active': grp_dose['V0_Active'].mean(),
                        'Mean_AUC': diff_mean, 'Delta_vs_Ref': diff_mean,
                        'CI_95': f"[{max(0.0, diff_mean - ci_margin):.2f}, {diff_mean + ci_margin:.2f}]",
                        'raw_p': pval
                    })

        self.df_table5 = pd.DataFrame(t5_records)
        if not self.df_table5.empty:
            q5 = adjust_pvalues(self.df_table5['raw_p'].values, fdr_method)
            self.df_table5['p_value'] = [f"{p:.4f}" if p >= 0.0001 else "< 0.0001" for p in self.df_table5['raw_p']]
            self.df_table5['q_Adj'] = [f"{q:.4f}" if q >= 0.0001 else "< 0.0001" for q in q5]

        self._fit_selected_kinetic_analyzer(analyzer_name, control_cell, n_boot)
        self._build_time_bin_overlap_dataframe(target_cell, control_cell)

        self._populate_treeview_from_df(self.tree_t0, self.df_preregistration)
        self._populate_treeview_from_df(self.tree_t1, self.df_v0_summary)
        self._populate_treeview_from_df(self.tree_t2, self.df_v0_600_summary)
        self._populate_treeview_from_df(self.tree_t3, self.df_v0_active_summary)
        self._populate_treeview_from_df(self.tree_tcat, self.df_categorical_summary)
        if not self.df_table4.empty:
            self._populate_treeview_from_df(self.tree_t4, self.df_table4.drop(columns=['raw_p'], errors='ignore'))
        if not self.df_table5.empty:
            self._populate_treeview_from_df(self.tree_t5, self.df_table5.drop(columns=['raw_p'], errors='ignore'))
        self._populate_treeview_from_df(self.tree_t6, self.df_saturation)
        if hasattr(self, 'tree_ttimebin'):
            self._populate_treeview_from_df(self.tree_ttimebin, self.df_time_bin_overlap)

    def _build_time_bin_overlap_dataframe(self, target_cell, control_cell):
        """
        Construye la Tabla de Datos por Bin de Tiempo (Time_Window) superponiendo/emparejando
        la Célula Diana (target_cell) y la Célula Control (control_cell) para cada Compuesto,
        Concentración (µM) y Bin de Tiempo (s), incluyendo permeabilización continua (% ± SEM),
        proporción de eventos afectados V2L (% ± SEM), no afectados V2R (%), morfología FSC (%)
        y el diferencial de selectividad (Delta Diana - Control).
        """
        if self.df_full is None or self.df_full.empty:
            self.df_time_bin_overlap = pd.DataFrame()
            return

        try:
            bin_sz = int(float(self.ent_bin_size.get()))
        except Exception:
            bin_sz = 15

        df_f = self.df_full.copy()
        if 'Applied_Scope' not in df_f.columns:
            df_f['Applied_Scope'] = 'Global'
        if 'Applied_V1L_Excl_Log10' not in df_f.columns:
            df_f['Applied_V1L_Excl_Log10'] = float(self.var_log_excl_v1l.get())
        if 'Applied_V2L_Cutoff_Log10' not in df_f.columns:
            df_f['Applied_V2L_Cutoff_Log10'] = float(self.var_log_cutoff.get())

        agg_cell = df_f.groupby(['Peptide', 'Concentration_uM', 'Time_Window', 'Cell_Line']).agg(
            Filter_Scope=('Applied_Scope', 'first'),
            V1L_Excl_Log10=('Applied_V1L_Excl_Log10', 'first'),
            V2L_Cutoff_Log10=('Applied_V2L_Cutoff_Log10', 'first'),
            N_Events=('Events', 'sum'),
            Perm_Mean=('Permeabilization_Pct', 'mean'),
            Perm_SEM=('Permeabilization_Pct', 'sem'),
            V2L_Mean=('Affected_V2L_Pct', 'mean'),
            V2L_SEM=('Affected_V2L_Pct', 'sem'),
            V2R_Mean=('Unaffected_V2R_Pct', 'mean'),
            FSC_Mean=('Norm_FSC_Pct', 'mean'),
            SSC_Mean=('Norm_SSC_Pct', 'mean')
        ).reset_index()
        agg_cell['Perm_SEM'] = agg_cell['Perm_SEM'].fillna(0.0)
        agg_cell['V2L_SEM'] = agg_cell['V2L_SEM'].fillna(0.0)

        # Mapa rápido de T_lag y T_max por (Peptide, Cell_Line, Concentration_uM)
        kin_map = {}
        if not self.df_kinetics.empty:
            for (kp, kc, kd), kgrp in self.df_kinetics.groupby(['Peptide', 'Cell_Line', 'Concentration_uM']):
                kin_map[(str(kp), str(kc), round(float(kd), 2))] = (
                    round(float(kgrp['T_lag'].mean()), 1),
                    round(float(kgrp['T_max'].mean()), 1),
                    str(kgrp['Tlag_Mode'].iloc[0])
                )

        sub_t = agg_cell[agg_cell['Cell_Line'] == target_cell].copy()
        sub_c = agg_cell[agg_cell['Cell_Line'] == control_cell].copy()

        keys = ['Peptide', 'Concentration_uM', 'Time_Window']
        merged = pd.merge(sub_t, sub_c, on=keys, how='outer', suffixes=('_Target', '_Control')).sort_values(keys)

        rows = []
        for _, r in merged.iterrows():
            pep_s = str(r['Peptide'])
            conc_f = round(float(r['Concentration_uM']), 2)
            tw = float(r['Time_Window'])
            t_perm = float(r['Perm_Mean_Target']) if pd.notna(r.get('Perm_Mean_Target')) else np.nan
            c_perm = float(r['Perm_Mean_Control']) if pd.notna(r.get('Perm_Mean_Control')) else np.nan
            t_v2l = float(r['V2L_Mean_Target']) if pd.notna(r.get('V2L_Mean_Target')) else np.nan
            c_v2l = float(r['V2L_Mean_Control']) if pd.notna(r.get('V2L_Mean_Control')) else np.nan

            d_perm = (t_perm - c_perm) if (np.isfinite(t_perm) and np.isfinite(c_perm)) else np.nan
            d_v2l = (t_v2l - c_v2l) if (np.isfinite(t_v2l) and np.isfinite(c_v2l)) else np.nan

            t_tlag, t_tmax, t_tmode = kin_map.get((pep_s, target_cell, conc_f), (np.nan, np.nan, "Auto"))
            c_tlag, c_tmax, c_tmode = kin_map.get((pep_s, control_cell, conc_f), (np.nan, np.nan, "Auto"))

            rows.append({
                'Peptide': pep_s,
                'Concentration_uM': conc_f,
                'Time_Bin_s': round(tw, 1),
                'Time_Interval': f"[{int(tw)}s, {int(tw + bin_sz)}s)",
                'Target_Cell': target_cell,
                'Target_Scope': str(r['Filter_Scope_Target']) if pd.notna(r.get('Filter_Scope_Target')) else "Global",
                'Target_V1L_Excl_Log10': round(float(r['V1L_Excl_Log10_Target']), 2) if pd.notna(r.get('V1L_Excl_Log10_Target')) else np.nan,
                'Target_V2L_Cutoff_Log10': round(float(r['V2L_Cutoff_Log10_Target']), 2) if pd.notna(r.get('V2L_Cutoff_Log10_Target')) else np.nan,
                'Target_Tlag_s': t_tlag,
                'Target_Tmax_s': t_tmax,
                'Target_Events_N': int(r['N_Events_Target']) if pd.notna(r.get('N_Events_Target')) else 0,
                'Target_Perm_Pct': round(t_perm, 2) if np.isfinite(t_perm) else np.nan,
                'Target_Perm_SEM': round(float(r['Perm_SEM_Target']), 2) if pd.notna(r.get('Perm_SEM_Target')) else 0.0,
                'Target_V2L_Pct': round(t_v2l, 2) if np.isfinite(t_v2l) else np.nan,
                'Target_V2L_SEM': round(float(r['V2L_SEM_Target']), 2) if pd.notna(r.get('V2L_SEM_Target')) else 0.0,
                'Target_V2R_Pct': round(float(r['V2R_Mean_Target']), 2) if pd.notna(r.get('V2R_Mean_Target')) else np.nan,
                'Target_FSC_Pct': round(float(r['FSC_Mean_Target']), 2) if pd.notna(r.get('FSC_Mean_Target')) else np.nan,
                'Control_Cell': control_cell,
                'Control_Scope': str(r['Filter_Scope_Control']) if pd.notna(r.get('Filter_Scope_Control')) else "Global",
                'Control_V1L_Excl_Log10': round(float(r['V1L_Excl_Log10_Control']), 2) if pd.notna(r.get('V1L_Excl_Log10_Control')) else np.nan,
                'Control_V2L_Cutoff_Log10': round(float(r['V2L_Cutoff_Log10_Control']), 2) if pd.notna(r.get('V2L_Cutoff_Log10_Control')) else np.nan,
                'Control_Tlag_s': c_tlag,
                'Control_Tmax_s': c_tmax,
                'Control_Events_N': int(r['N_Events_Control']) if pd.notna(r.get('N_Events_Control')) else 0,
                'Control_Perm_Pct': round(c_perm, 2) if np.isfinite(c_perm) else np.nan,
                'Control_Perm_SEM': round(float(r['Perm_SEM_Control']), 2) if pd.notna(r.get('Perm_SEM_Control')) else 0.0,
                'Control_V2L_Pct': round(c_v2l, 2) if np.isfinite(c_v2l) else np.nan,
                'Control_V2L_SEM': round(float(r['V2L_SEM_Control']), 2) if pd.notna(r.get('V2L_SEM_Control')) else 0.0,
                'Control_V2R_Pct': round(float(r['V2R_Mean_Control']), 2) if pd.notna(r.get('V2R_Mean_Control')) else np.nan,
                'Control_FSC_Pct': round(float(r['FSC_Mean_Control']), 2) if pd.notna(r.get('FSC_Mean_Control')) else np.nan,
                'Delta_Perm_Target_vs_Ctrl_Pct': round(d_perm, 2) if np.isfinite(d_perm) else np.nan,
                'Delta_V2L_Target_vs_Ctrl_Pct': round(d_v2l, 2) if np.isfinite(d_v2l) else np.nan
            })

        self.df_time_bin_overlap = pd.DataFrame(rows)

    def _evaluate_model_fit_for_group(self, pep, cell, grp, analyzer_name, control_cell, rng, n_boot=0):
        """
        Ajusta la función matemática especificada por analyzer_name sobre el grupo de muestras (pep, cell),
        calculando la curva de la función f(x), la ecuación numérica explícita evaluada y los coeficientes
        de determinación R² (sobre medias del grupo y sobre todas las réplicas individuales), R²_adj y RMSE.
        """
        agg = grp.groupby('Concentration_uM').agg(
            V0_Active=('V0_Active', 'mean'),
            V0_SEM=('V0_Active', 'sem'),
            Normalized_AUC_Pct=('Normalized_AUC_Pct', 'mean')
        ).reset_index()
        x = agg['Concentration_uM'].values.astype(float)
        y = agg['V0_Active'].values.astype(float)
        y_sem = agg['V0_SEM'].fillna(0.0).values.astype(float)
        raw_x = grp['Concentration_uM'].values.astype(float)
        raw_y = grp['V0_Active'].values.astype(float)
        max_v0 = float(np.max(y)) if len(y) > 0 else 0.0
        n_doses = len(x)

        if "Bi-Exponencial" in analyzer_name:
            max_d = float(np.max(x)) if len(x) > 0 else 1.0
            sub_full = self.df_full[
                (self.df_full['Peptide'] == pep) &
                (self.df_full['Cell_Line'] == cell) &
                (np.isclose(self.df_full['Concentration_uM'], max_d, atol=1e-2)) &
                (self.df_full['Time_Window'] >= 0)
            ]
            sub_t = sub_full.groupby('Time_Window').agg(
                p_m=('Permeabilization_Pct', 'mean'),
                p_se=('Permeabilization_Pct', 'sem')
            ).reset_index()
            t_vec = sub_t['Time_Window'].values.astype(float)
            p_vec = sub_t['p_m'].values.astype(float)
            p_sem = sub_t['p_se'].fillna(0.0).values.astype(float)
            raw_t = sub_full['Time_Window'].values.astype(float)
            raw_p = sub_full['Permeabilization_Pct'].values.astype(float)
            try:
                p0_bi = [max(1.0, min(15.0, float(np.max(p_vec)) * 0.2)), 0.04, max(2.0, float(np.max(p_vec)) * 0.8), 0.015, 75.0]
                b_bi = ([0.0, 0.0005, 0.0, 0.0002, 0.0], [50.0, 0.8, 100.0, 0.4, 540.0])
                popt_bi, _ = curve_fit(biexponential_time_model, t_vec, p_vec, p0=p0_bi, bounds=b_bi, maxfev=20000)
                p_ads, k_ads, p_lys, k_lys, tlag_bi = popt_bi
                t_grid = np.linspace(0.0, float(np.max(t_vec)) if len(t_vec) > 0 else 600.0, 100)
                p_smooth = biexponential_time_model(t_grid, *popt_bi)
                p_pred = biexponential_time_model(t_vec, *popt_bi)
                res_bi = p_vec - p_pred
                ss_tot_bi = float(np.sum((p_vec - np.mean(p_vec)) ** 2))
                r2_grp = max(0.0, 1.0 - float(np.sum(res_bi ** 2)) / ss_tot_bi) if ss_tot_bi > 1e-9 else 1.0
                n_pts_t = len(t_vec)
                r2_adj = max(0.0, 1.0 - (1.0 - r2_grp) * (n_pts_t - 1) / max(1, n_pts_t - 5)) if n_pts_t > 5 else r2_grp
                rmse_val = float(np.sqrt(np.mean(res_bi ** 2)))

                raw_pred = biexponential_time_model(raw_t, *popt_bi)
                ss_tot_raw = float(np.sum((raw_p - np.mean(raw_p)) ** 2))
                r2_rep = max(0.0, 1.0 - float(np.sum((raw_p - raw_pred) ** 2)) / ss_tot_raw) if ss_tot_raw > 1e-9 else r2_grp

                p_at_tlag = float(biexponential_time_model(np.array([tlag_bi]), *popt_bi)[0])
                eq_eval = f"P(t) = {p_ads:.1f}(1-e^-{k_ads:.3f}t) + {p_lys:.1f}(1-e^-{k_lys:.3f}·max(0,t-{tlag_bi:.0f}))"
                return {
                    'pep': pep, 'cell': cell, 'model_lbl': 'Bi-Exponencial (t)',
                    'x_pts': t_vec, 'y_pts': p_vec, 'y_sem': p_sem,
                    'raw_x': raw_t, 'raw_y': raw_p,
                    'x_fit': t_grid, 'y_fit': p_smooth,
                    'kd_marker': tlag_bi, 'y_half': p_at_tlag,
                    'eq_eval': eq_eval,
                    'r2_grp': r2_grp, 'r2_rep': r2_rep, 'r2_adj': r2_adj, 'rmse': rmse_val,
                    'r2_badge': f"[R²_grupo = {r2_grp:.4f} | R²_rép = {r2_rep:.4f}]",
                    'sat_rec': {
                        'Peptide': pep, 'Cell_Line': cell,
                        'Vmax_pct_s': f"k_lys={k_lys:.4f} s⁻¹",
                        'Kd_or_Half_uM': f"T_lag={tlag_bi:.1f} s",
                        'Cooperativity_Shape': f"P_ads={p_ads:.1f}%",
                        'Lytic_Efficiency': f"k_ads={k_ads:.3f} s⁻¹",
                        'R_Squared': f"{r2_grp:.4f}", 'R2_Adj': f"{r2_adj:.4f}",
                        'Model_Type': 'Bi-Exponencial (Dosis Máx)'
                    }
                }
            except Exception:
                pass

        x_grid = np.linspace(0.0, max(x) * 1.15 if len(x) > 0 else 5.0, 80)

        # Si la respuesta es prácticamente basal/inerte o tiene menos de 3 dosis, ajustar modelo basal/lineal o Hill suave
        if max_v0 < 0.006 or n_doses < 3:
            if len(x) >= 2 and np.std(x) > 1e-9:
                slope, intercept, r_val, _, _ = stats.linregress(x, y)
                y_pred = intercept + slope * x
                y_smooth = np.clip(intercept + slope * x_grid, 0.0, None)
                res = y - y_pred
                r2_grp = float(r_val ** 2)
                rmse_val = float(np.sqrt(np.mean(res ** 2)))
                raw_pred = intercept + slope * raw_x
                ss_raw = float(np.sum((raw_y - np.mean(raw_y)) ** 2))
                r2_rep = max(0.0, 1.0 - float(np.sum((raw_y - raw_pred) ** 2)) / ss_raw) if ss_raw > 1e-12 else r2_grp
                eq_eval = f"V(C) = {slope:.4f}·C + {intercept:.4f}  [Basal/Lineal]"
            else:
                mean_y = float(np.mean(y)) if len(y) > 0 else 0.0
                y_smooth = np.full_like(x_grid, mean_y)
                r2_grp, r2_rep, rmse_val = 1.0, 1.0, 0.0
                eq_eval = f"V(C) = {mean_y:.4f} (Constante Basal)"

            return {
                'pep': pep, 'cell': cell, 'model_lbl': 'Basal / Lineal',
                'x_pts': x, 'y_pts': y, 'y_sem': y_sem,
                'raw_x': raw_x, 'raw_y': raw_y,
                'x_fit': x_grid, 'y_fit': y_smooth,
                'kd_marker': None, 'y_half': None,
                'eq_eval': eq_eval,
                'r2_grp': r2_grp, 'r2_rep': r2_rep, 'r2_adj': r2_grp, 'rmse': rmse_val,
                'r2_badge': f"[R²_grupo = {r2_grp:.4f} | Control/Basal]",
                'sat_rec': {
                    'Peptide': pep, 'Cell_Line': cell,
                    'Vmax_pct_s': f"{max_v0:.4f}" if max_v0 >= 0.001 else "< 0.001",
                    'Kd_or_Half_uM': '> Dosis Máx', 'Cooperativity_Shape': '1.00',
                    'Lytic_Efficiency': '0.000',
                    'R_Squared': f"{r2_grp:.4f}", 'R2_Adj': f"{r2_grp:.4f}",
                    'Model_Type': 'Basal / Lineal (Inerte)'
                }
            }

        try:
            if "Michaelis-Menten" in analyzer_name:
                p0 = [max_v0 * 1.2, float(np.median(x))]
                bounds = ([max_v0 * 0.3, 0.01], [max_v0 * 4.0, max(x) * 4.0])
                popt, _ = curve_fit(michaelis_menten_model, x, y, p0=p0, bounds=bounds, maxfev=12000)
                vmax_f, kd_f = popt
                nh_f = 1.0
                y_pred = michaelis_menten_model(x, *popt)
                y_smooth = michaelis_menten_model(x_grid, *popt)
                raw_pred = michaelis_menten_model(raw_x, *popt)
                y_half_val = float(michaelis_menten_model(kd_f, *popt))
                n_params = 2
                model_lbl = "Michaelis-Menten"
                eq_eval = f"V(C) = ({vmax_f:.3f} · C) / ({kd_f:.2f} + C)"
            elif "Gompertz" in analyzer_name:
                p0 = [max_v0 * 1.1, float(np.median(x)), 1.5]
                bounds = ([max_v0 * 0.4, 0.01, 0.1], [max_v0 * 3.0, max(x) * 2.5, 10.0])
                popt, _ = curve_fit(gompertz_dose_model, x, y, p0=p0, bounds=bounds, maxfev=15000)
                vmax_f, kd_f, nh_f = popt
                y_pred = gompertz_dose_model(x, *popt)
                y_smooth = gompertz_dose_model(x_grid, *popt)
                raw_pred = gompertz_dose_model(raw_x, *popt)
                y_half_val = float(gompertz_dose_model(kd_f, *popt))
                n_params = 3
                model_lbl = "Gompertz Sigmoidal"
                eq_eval = f"V(C) = {vmax_f:.3f} · exp(-exp(-{nh_f:.2f}·(C - {kd_f:.2f})))"
            elif "Weibull" in analyzer_name:
                p0 = [max_v0 * 1.1, float(np.median(x)), 2.0]
                bounds = ([max_v0 * 0.4, 0.01, 0.4], [max_v0 * 3.0, max(x) * 2.5, 6.0])
                popt, _ = curve_fit(weibull_dose_model, x, y, p0=p0, bounds=bounds, maxfev=15000)
                vmax_f, kd_f, nh_f = popt
                y_pred = weibull_dose_model(x, *popt)
                y_smooth = weibull_dose_model(x_grid, *popt)
                raw_pred = weibull_dose_model(raw_x, *popt)
                y_half_val = float(weibull_dose_model(kd_f, *popt))
                n_params = 3
                model_lbl = "Weibull Stretched"
                eq_eval = f"V(C) = {vmax_f:.3f} · [1 - exp(-(C / {kd_f:.2f})^{nh_f:.2f})]"
            else:
                p0 = [max_v0 * 1.1, float(np.median(x)), 2.0]
                bounds = ([max_v0 * 0.4, 0.01, 0.5], [max_v0 * 3.0, max(x) * 2.5, 6.0])
                popt, _ = curve_fit(hill_model, x, y, p0=p0, bounds=bounds, maxfev=20000)
                vmax_f, kd_f, nh_f = popt
                y_pred = hill_model(x, *popt)
                y_smooth = hill_model(x_grid, *popt)
                raw_pred = hill_model(raw_x, *popt)
                y_half_val = float(hill_model(kd_f, *popt))
                n_params = 3
                model_lbl = "Allosteric Hill"
                eq_eval = f"V(C) = ({vmax_f:.3f} · C^{nh_f:.2f}) / ({kd_f:.2f}^{nh_f:.2f} + C^{nh_f:.2f})"

            res = y - y_pred
            ss_tot = float(np.sum((y - np.mean(y)) ** 2))
            r2_grp = max(0.0, 1.0 - (float(np.sum(res ** 2)) / ss_tot)) if ss_tot > 1e-12 else 1.0
            df_mod = n_doses - n_params
            r2_adj = max(0.0, 1.0 - (1.0 - r2_grp) * (n_doses - 1) / df_mod) if df_mod >= 1 else r2_grp
            rmse_val = float(np.sqrt(np.mean(res ** 2)))

            ss_tot_raw = float(np.sum((raw_y - np.mean(raw_y)) ** 2))
            r2_rep = max(0.0, 1.0 - (float(np.sum((raw_y - raw_pred) ** 2)) / ss_tot_raw)) if ss_tot_raw > 1e-12 else r2_grp

            boot_kd = []
            if "Hill" in model_lbl and n_boot > 0:
                for _ in range(n_boot):
                    y_b = [rng.choice(grp[grp['Concentration_uM'] == d]['V0_Active'].values, size=1)[0] for d in x]
                    try:
                        pb, _ = curve_fit(hill_model, x, y_b, p0=popt, bounds=bounds, maxfev=1200)
                        boot_kd.append(pb[1])
                    except Exception:
                        pass
            kd_str = f"{kd_f:.2f} [{np.percentile(boot_kd, 2.5):.2f}-{np.percentile(boot_kd, 97.5):.2f}]" if len(boot_kd) > 15 else f"{kd_f:.3f}"

            return {
                'pep': pep, 'cell': cell, 'model_lbl': model_lbl,
                'x_pts': x, 'y_pts': y, 'y_sem': y_sem,
                'raw_x': raw_x, 'raw_y': raw_y,
                'x_fit': x_grid, 'y_fit': y_smooth,
                'kd_marker': kd_f, 'y_half': y_half_val,
                'eq_eval': eq_eval,
                'r2_grp': r2_grp, 'r2_rep': r2_rep, 'r2_adj': r2_adj, 'rmse': rmse_val,
                'r2_badge': f"[R²_grupo = {r2_grp:.4f} | R²_rép = {r2_rep:.4f} | Kd = {kd_f:.2f}µM]",
                'sat_rec': {
                    'Peptide': pep, 'Cell_Line': cell,
                    'Vmax_pct_s': f"{vmax_f:.4f}",
                    'Kd_or_Half_uM': kd_str,
                    'Cooperativity_Shape': f"{nh_f:.2f}",
                    'Lytic_Efficiency': f"{(vmax_f / max(1e-6, kd_f)):.4f}",
                    'R_Squared': f"{r2_grp:.4f}",
                    'R2_Adj': f"{r2_adj:.4f}",
                    'Model_Type': model_lbl
                }
            }
        except Exception:
            slope, intercept, r_val, _, _ = stats.linregress(x, y)
            y_pred = intercept + slope * x
            y_smooth = np.clip(intercept + slope * x_grid, 0.0, None)
            r2_grp = float(r_val ** 2)
            rmse_val = float(np.sqrt(np.mean((y - y_pred) ** 2)))
            eq_eval = f"V(C) = {slope:.4f}·C + {intercept:.4f}"
            return {
                'pep': pep, 'cell': cell, 'model_lbl': 'Aprox. Lineal',
                'x_pts': x, 'y_pts': y, 'y_sem': y_sem,
                'raw_x': raw_x, 'raw_y': raw_y,
                'x_fit': x_grid, 'y_fit': y_smooth,
                'kd_marker': None, 'y_half': None,
                'eq_eval': eq_eval,
                'r2_grp': r2_grp, 'r2_rep': r2_grp, 'r2_adj': r2_grp, 'rmse': rmse_val,
                'r2_badge': f"[R²_grupo = {r2_grp:.4f}]",
                'sat_rec': {
                    'Peptide': pep, 'Cell_Line': cell,
                    'Vmax_pct_s': f"{max_v0:.4f}", 'Kd_or_Half_uM': 'N/A',
                    'Cooperativity_Shape': '1.00', 'Lytic_Efficiency': f"{slope:.4f}",
                    'R_Squared': f"{r2_grp:.4f}", 'R2_Adj': f"{r2_grp:.4f}",
                    'Model_Type': 'Linear Approx.'
                }
            }

    def _fit_selected_kinetic_analyzer(self, analyzer_name, control_cell, n_boot=200):
        saturation_records = []
        preview_curves = []
        eq_table_rows = []
        unique_peps = list(self.df_kinetics['Peptide'].unique())
        rng = np.random.default_rng(42)

        if "Michaelis-Menten" in analyzer_name:
            eq_str = "V_active(C) = (V_max · C) / (K_d + C)   [n_H = 1.0 fijo]"
            x_lbl, y_lbl = "Concentración (µM)", "V_active (%/s)"
        elif "Gompertz" in analyzer_name:
            eq_str = "V_active(C) = V_max · exp( -exp( -k_G · (C - K_d) ) )"
            x_lbl, y_lbl = "Concentración (µM)", "V_active (%/s)"
        elif "Weibull" in analyzer_name:
            eq_str = "V_active(C) = V_max · [ 1 - exp( -(C / K_d)^m ) ]"
            x_lbl, y_lbl = "Concentración (µM)", "V_active (%/s)"
        elif "Bi-Exponencial" in analyzer_name:
            eq_str = "P(t) = P_ads·(1 - e^(-k_ads·t)) + P_lys·(1 - e^(-k_lys·max(0, t - T_lag)))"
            x_lbl, y_lbl = "Tiempo post-adición (s)", "% Permeabilización (Dosis Máx)"
        else:
            eq_str = "V_active(C) = (V_max · C^nH) / (K_d^nH + C^nH)   [0.5 <= nH <= 6.0]"
            x_lbl, y_lbl = "Concentración (µM)", "V_active (%/s)"

        all_groups = sorted(list(self.df_kinetics.groupby(['Peptide', 'Cell_Line']).groups.keys()))
        group_labels = ["Todos los Grupos"] + [f"{p} | {c}" for p, c in all_groups]
        if hasattr(self, 'cmb_preview_group'):
            self.cmb_preview_group['values'] = group_labels
            if self.var_preview_group.get() not in group_labels:
                self.var_preview_group.set("Todos los Grupos")

        sel_grp_str = self.var_preview_group.get() if hasattr(self, 'var_preview_group') else "Todos los Grupos"
        prev_mode_str = self.var_preview_mode.get() if hasattr(self, 'var_preview_mode') else "Modelo Seleccionado (Por Grupo)"

        # 1. Ajustar el modelo seleccionado a cada grupo para la Tabla 6 y recolectar resultados
        group_fit_map = {}
        for idx_g, (pep, cell) in enumerate(all_groups):
            grp = self.df_kinetics[(self.df_kinetics['Peptide'] == pep) & (self.df_kinetics['Cell_Line'] == cell)]
            fit_res = self._evaluate_model_fit_for_group(pep, cell, grp, analyzer_name, control_cell, rng, n_boot=n_boot)
            group_fit_map[(pep, cell)] = fit_res
            saturation_records.append(fit_res['sat_rec'])

        self.df_saturation = pd.DataFrame(saturation_records)

        # 2. Construir las curvas y tabla del Panel de Previsualización según el modo elegido
        if "Comparar Todos los Modelos" in prev_mode_str and all_groups:
            # Si está en "Todos los Grupos", elegir el primer grupo de la línea celular diana (Target)
            target_c = self.ent_target.get().strip()
            if sel_grp_str != "Todos los Grupos" and "|" in sel_grp_str:
                pts = [s.strip() for s in sel_grp_str.split("|")]
                chosen_key = (pts[0], pts[1]) if (pts[0], pts[1]) in group_fit_map else all_groups[0]
            else:
                target_matches = [k for k in all_groups if k[1] == target_c]
                chosen_key = target_matches[0] if target_matches else all_groups[0]

            pep_c, cell_c = chosen_key
            grp_c = self.df_kinetics[(self.df_kinetics['Peptide'] == pep_c) & (self.df_kinetics['Cell_Line'] == cell_c)]
            compare_models = [
                ("Allosteric Hill (Regularizado nH <= 6.0)", "#004488", None),
                ("Michaelis-Menten Hiperbólico (nH = 1.0)", "#E66100", (5, 3)),
                ("Gompertz Sigmoidal Asimétrico", "#117733", (3, 2)),
                ("Weibull / Exponencial Estirado", "#BB5566", (6, 2, 2, 2)),
            ]
            x_lbl, y_lbl = "Concentración (µM)", "V_active (%/s)"
            eq_str = "Comparativa de Funciones: Hill vs Michaelis-Menten vs Gompertz vs Weibull"
            r2_vals_cmp = []
            for m_idx, (m_name, m_col, m_dash) in enumerate(compare_models):
                res_m = self._evaluate_model_fit_for_group(pep_c, cell_c, grp_c, m_name, control_cell, rng, n_boot=0)
                r2_vals_cmp.append(res_m['r2_grp'])
                preview_curves.append({
                    'label': f"{res_m['model_lbl']}",
                    'color': m_col,
                    'dash': m_dash,
                    'x_pts': res_m['x_pts'] if m_idx == 0 else [],
                    'y_pts': res_m['y_pts'] if m_idx == 0 else [],
                    'y_sem': res_m['y_sem'] if m_idx == 0 else [],
                    'raw_x': res_m['raw_x'] if m_idx == 0 else [],
                    'raw_y': res_m['raw_y'] if m_idx == 0 else [],
                    'x_fit': res_m['x_fit'],
                    'y_fit': res_m['y_fit'],
                    'kd_marker': res_m['kd_marker'] if m_idx == 0 else None,
                    'y_half': res_m['y_half'] if m_idx == 0 else None,
                    'eq_eval': res_m['eq_eval'],
                    'r2_badge': f"R²_grupo = {res_m['r2_grp']:.4f} | R²_rép = {res_m['r2_rep']:.4f}"
                })
                eq_table_rows.append((
                    f"{pep_c} | {cell_c}",
                    res_m['model_lbl'],
                    res_m['eq_eval'],
                    f"{res_m['r2_grp']:.4f}",
                    f"{res_m['r2_rep']:.4f}",
                    f"{res_m['r2_adj']:.4f}",
                    f"{res_m['rmse']:.4f}"
                ))
            best_idx = int(np.argmax(r2_vals_cmp)) if r2_vals_cmp else 0
            best_m_lbl = preview_curves[best_idx]['label'] if preview_curves else "Hill"
            best_r2 = max(r2_vals_cmp) if r2_vals_cmp else 0.0
            global_r2_banner = f"Grupo: {pep_c} ({cell_c}) | Mejor Ajuste: {best_m_lbl} (R² = {best_r2:.4f})"
            plot_title = f"Comparación Multi-Modelo sobre Grupo: {pep_c} ({cell_c})"
        else:
            # Modo estándar: graficar la función del modelo seleccionado para el grupo o todos los grupos
            r2_grp_list = []
            r2_rep_list = []
            for idx_g, (pep, cell) in enumerate(all_groups):
                grp_lbl = f"{pep} | {cell}"
                if sel_grp_str != "Todos los Grupos" and grp_lbl != sel_grp_str:
                    continue
                res_g = group_fit_map[(pep, cell)]
                col_g = CURVE_PALETTE[idx_g % len(CURVE_PALETTE)]
                dash_g = (4, 3) if cell == control_cell else None
                preview_curves.append({
                    'label': f"{pep} ({cell})",
                    'color': col_g,
                    'dash': dash_g,
                    'x_pts': res_g['x_pts'],
                    'y_pts': res_g['y_pts'],
                    'y_sem': res_g['y_sem'],
                    'raw_x': res_g['raw_x'],
                    'raw_y': res_g['raw_y'],
                    'x_fit': res_g['x_fit'],
                    'y_fit': res_g['y_fit'],
                    'kd_marker': res_g['kd_marker'],
                    'y_half': res_g['y_half'],
                    'eq_eval': res_g['eq_eval'],
                    'r2_badge': res_g['r2_badge']
                })
                r2_grp_list.append(res_g['r2_grp'])
                r2_rep_list.append(res_g['r2_rep'])
                eq_table_rows.append((
                    grp_lbl,
                    res_g['model_lbl'],
                    res_g['eq_eval'],
                    f"{res_g['r2_grp']:.4f}",
                    f"{res_g['r2_rep']:.4f}",
                    f"{res_g['r2_adj']:.4f}",
                    f"{res_g['rmse']:.4f}"
                ))

            mean_r2_grp = float(np.mean(r2_grp_list)) if r2_grp_list else 0.0
            mean_r2_rep = float(np.mean(r2_rep_list)) if r2_rep_list else 0.0
            global_r2_banner = (
                f"R² Medio del Grupo de Muestras = {mean_r2_grp:.4f}  |  "
                f"R² Global sobre Réplicas = {mean_r2_rep:.4f}  ({len(preview_curves)} grupos)"
            )
            plot_title = f"Función Ajustada: {analyzer_name}"

        if hasattr(self, 'tree_model_eq_r2'):
            for item in self.tree_model_eq_r2.get_children():
                self.tree_model_eq_r2.delete(item)
            for row_vals in eq_table_rows:
                self.tree_model_eq_r2.insert("", "end", values=row_vals)

        self.sat_preview.set_curves(
            preview_curves,
            title=plot_title,
            equation_str=eq_str,
            x_label=x_lbl,
            y_label=y_lbl,
            global_r2_str=global_r2_banner
        )
        if hasattr(self, 'lbl_model_preview_summary'):
            eq_lines = [f"• {r[0]} [{r[1]}]:  {r[2]}   ➔   R²_grupo = {r[3]} (R²_réplicas = {r[4]})" for r in eq_table_rows[:3]]
            sum_txt = f"{global_r2_banner}\n" + ("\n".join(eq_lines) if eq_lines else "Sin grupos activos.")
            self.lbl_model_preview_summary.config(text=sum_txt)

    # --------------------------------------------------------------------------
    # EXPORTACIÓN DE REPORTE FINAL AUDITABLE POR CATEGORÍAS CON TODAS LAS IMÁGENES
    # --------------------------------------------------------------------------
    def export_final_report_bundle(self):
        if self.df_full.empty or self.df_kinetics.empty:
            messagebox.showwarning("Sin Datos Procesados", "Ejecuta primero el Pre-Análisis antes de exportar el reporte.", parent=self.root)
            return

        out_dir = Path(self.ent_outdir.get().strip())
        out_dir.mkdir(parents=True, exist_ok=True)

        # Sincronizar todas las muestras y recalcular tablas con los Cut-Offs específicos configurados por el usuario
        self._sync_all_samples_with_active_user_cutoffs()
        self.recompute_all_models_and_tables()

        target_cell = self.ent_target.get().strip()
        control_cell = self.ent_control.get().strip()
        excl_v1l_log = float(self.var_log_excl_v1l.get())
        cutoff_log = float(self.var_log_cutoff.get())

        exported_files = []
        for name, df_out in [
            ("Category0_Preregistration_PreAnalysis_Parameters", self.df_preregistration),
            ("Category0B_User_Configured_Sample_CutOffs_Table", getattr(self, 'df_user_cutoffs_table', self.df_preregistration)),
            ("Category1A_Categorical_Effect_Milestones_T50_Tmax_Stats", self.df_categorical_summary),
            ("Category1B_Categorical_Proportions_By_Time_Window", self.df_cat_time_windows),
            ("Category1C_Categorical_Hypothesis_Tests_And_Isochrones", self.df_cat_hypothesis_tests),
            ("Category2_Table1_Early_Rate", self.df_v0_summary),
            ("Category2_Table2_Full_Rate", self.df_v0_600_summary),
            ("Category2_Table3_Active_Rate", self.df_v0_active_summary),
            ("Category2_Time_Bin_Overlapping_Target_vs_Control_Table", self.df_time_bin_overlap),
            ("Category3_Morphology_FSC_SSC", self.df_morphology_summary),
            ("Category4_Model_Parameters", self.df_saturation),
            ("Category5_Table4_Selectivity_Contrasts", self.df_table4.drop(columns=['raw_p'], errors='ignore')),
            ("Category5_Table5_Dynamic_Contrasts", self.df_table5.drop(columns=['raw_p'], errors='ignore')),
            ("Raw_Tube_Kinetics", self.df_kinetics)
        ]:
            csv_p = out_dir / f"{name}.csv"
            tsv_p = out_dir / f"{name}.tsv"
            try:
                df_out.to_csv(csv_p, index=False)
            except PermissionError:
                ts = datetime.now().strftime("%Y%m%d_%H%M%S")
                csv_p = out_dir / f"{name}_{ts}.csv"
                df_out.to_csv(csv_p, index=False)
            try:
                df_out.to_csv(tsv_p, sep="\t", index=False)
            except PermissionError:
                ts = datetime.now().strftime("%Y%m%d_%H%M%S")
                tsv_p = out_dir / f"{name}_{ts}.tsv"
                df_out.to_csv(tsv_p, sep="\t", index=False)
            exported_files.append(str(csv_p))

        # ----------------------------------------------------------------------
        # GENERACIÓN DE TODAS LAS FIGURAS EN FORMATO TIFF A 600 DPI Y BUFFERS RGB PARA EL PDF
        # ----------------------------------------------------------------------
        cat1_images = []
        cat2_images = []
        cat3_images = []
        cat4_images = []
        cat5_images = []

        # Resumen de rangos de Cut-Offs configurados por el usuario en todas las muestras
        if not self.df_kinetics.empty:
            v1l_u = sorted(self.df_kinetics['V1L_Excl_Log10'].round(2).unique())
            v2l_u = sorted(self.df_kinetics['V2L_Cutoff_Log10'].round(2).unique())
            v1l_range_str = f"{v1l_u[0]:.2f}" if len(v1l_u) == 1 else f"{v1l_u[0]:.2f}..{v1l_u[-1]:.2f}"
            v2l_range_str = f"{v2l_u[0]:.2f}" if len(v2l_u) == 1 else f"{v2l_u[0]:.2f}..{v2l_u[-1]:.2f}"
        else:
            v1l_range_str = f"{excl_v1l_log:.2f}"
            v2l_range_str = f"{cutoff_log:.2f}"

        # 0. Figura Comparativa Superpuesta de Célula Diana (target_cell) vs. Célula Control (control_cell) por Bin de Tiempo (TIFF 600 DPI)
        overlap_fig_path = out_dir / "Figure_Category2_Target_vs_Control_Overlapping_Kinetics_600DPI.tiff"
        overlap_rgb = render_highres_target_vs_control_overlap_png(
            overlap_fig_path, self.df_full, target_cell=target_cell, control_cell=control_cell, df_kinetics=self.df_kinetics
        )
        if overlap_rgb is not None:
            exported_files.append(str(overlap_fig_path))
            cat2_images.append((
                f"Figura II.0: Grafica Superpuesta por Bin de Tiempo — Celula Diana ({target_cell}, solido) vs. Celula Control ({control_cell}, discontinuo) [Cut V2L<{v2l_range_str}]",
                overlap_rgb
            ))
            cat5_images.append((
                f"Figura V.1: Comparacion Superpuesta de Selectividad por Bin Temporal ({target_cell} vs {control_cell} — Continua y V2L% con Cut-Offs Configurados)",
                overlap_rgb
            ))

        # 1. Figura de Histograma y Doble Cut-Off Log10(Fluorescencia) de las Muestras Seleccionadas (TIFF 600 DPI)
        sel_peps, sel_cells, sel_doses = self._get_selected_sample_filters()
        y_mode = "pct" if "Porcentaje" in self.var_hist_ymode.get() else "count"
        bins_edges, hist_series = self._build_histogram_series_data(sel_peps, sel_cells, sel_doses, y_mode=y_mode)
        hist_fig_path = out_dir / "Figure_Category1_Histogram_Cutoff_600DPI.tiff"
        hist_rgb = render_highres_histogram_png(
            hist_fig_path, bins_edges, hist_series, cutoff_log,
            f"Histograma de Fluorescencia: Exclusión V1L (<{v1l_range_str}) y Cut-Off V2L/V2R (<{v2l_range_str} log10)",
            "Eventos (%)" if y_mode == "pct" else "N° Eventos",
            excl_log=excl_v1l_log
        )
        exported_files.append(str(hist_fig_path))
        cat1_images.append((
            f"Figura I.1: Histograma Log10(Fluorescencia), Cut-Off 1 Exclusión V1L (<{v1l_range_str}) y Cut-Off 2 V2L/V2R (<{v2l_range_str} log10)",
            hist_rgb
        ))

        # 1B. Figura de Hitos Temporales (T_onset, T50 del Efecto, T_efecto_máx) y Balance de Proporciones V2L vs V2R (TIFF 600 DPI)
        mile_fig_path = out_dir / "Figure_Category1_Milestones_T50_Tmax_And_Proportions_600DPI.tiff"
        mile_rgb = render_highres_categorical_milestones_png(mile_fig_path, self.df_categorical_summary, target_cell=target_cell)
        if mile_rgb is not None:
            exported_files.append(str(mile_fig_path))
            cat1_images.append((
                f"Figura I.2: Hitos Temporales del Efecto (T_onset, T50 del Efecto, T_efecto_max) y Balance V2L/V2R con Cut-Offs Configurados ({target_cell})",
                mile_rgb
            ))

        # 1C. Figura de Mapa de Calor Espaciotemporal (% Eventos Afectados V2L por Compuesto, Dosis y Ventana de Tiempo, TIFF 600 DPI)
        hmap_fig_path = out_dir / "Figure_Category1_Spacetime_Heatmap_V2L_600DPI.tiff"
        hmap_rgb = render_highres_categorical_heatmap_png(hmap_fig_path, self.df_cat_time_windows, target_cell=target_cell)
        if hmap_rgb is not None:
            exported_files.append(str(hmap_fig_path))
            cat1_images.append((
                f"Figura I.3: Mapa de Calor Espaciotemporal de Proporcion de Eventos Afectados (% V2L) por Ventana de Tiempo y Cut-Off ({target_cell})",
                hmap_rgb
            ))

        # 2. Exportar e incrustar Figuras de Cinética Continua, Categórica (V2L%) y Morfológica (FSC) por Compuesto (TIFF 600 DPI)
        unique_peps = list(self.df_full['Peptide'].unique())
        for idx_p, pep in enumerate(unique_peps):
            sub = self.df_full[self.df_full['Peptide'] == pep].copy()
            unique_doses = sorted(sub['Concentration_uM'].unique())
            ramp = DOSE_RAMPS.get(pep, GENERIC_RAMPS[idx_p % len(GENERIC_RAMPS)])
            safe_pep = re.sub(r'[^a-zA-Z0-9_-]', '_', str(pep))

            series_cont = []
            series_cat_v2l = []
            series_fsc = []
            vlines_export = [(0.0, '#333333', True)]

            for (cell, conc), grp in sub.groupby(['Cell_Line', 'Concentration_uM']):
                d_idx = unique_doses.index(conc) if conc in unique_doses else 0
                col_hex = ramp[min(d_idx, len(ramp) - 1)]
                agg = grp.groupby('Time_Window').agg(
                    m_perm=('Permeabilization_Pct', 'mean'), se_perm=('Permeabilization_Pct', 'sem'),
                    m_v2l=('Affected_V2L_Pct', 'mean'), se_v2l=('Affected_V2L_Pct', 'sem'),
                    m_fsc=('Norm_FSC_Pct', 'mean'), se_fsc=('Norm_FSC_Pct', 'sem')
                ).reset_index()

                # Obtener los Cut-Offs específicos configurados por el usuario para esta muestra (V1L, V2L, T_lag, T_max)
                s_v2l = float(grp['Applied_V2L_Cutoff_Log10'].iloc[0]) if 'Applied_V2L_Cutoff_Log10' in grp else cutoff_log
                s_tlag, s_tmax, s_tmode = 60.0, 300.0, "Auto"
                if not self.df_kinetics.empty:
                    k_sub = self.df_kinetics[
                        (self.df_kinetics['Peptide'] == pep) &
                        (self.df_kinetics['Cell_Line'] == cell) &
                        (np.isclose(self.df_kinetics['Concentration_uM'], conc, atol=1e-2))
                    ]
                    if not k_sub.empty:
                        s_tlag = float(k_sub['T_lag'].mean())
                        s_tmax = float(k_sub['T_max'].mean())
                        s_tmode = str(k_sub['Tlag_Mode'].iloc[0])

                lbl_cont_curve = f"{conc:.2f} uM ({cell} | Tlag={s_tlag:.0f}s, Tmax={s_tmax:.0f}s [{s_tmode}])"
                lbl_cat_curve = f"{conc:.2f} uM ({cell} | Cut<{s_v2l:.2f} | Tlag={s_tlag:.0f}s, Tmax={s_tmax:.0f}s)"
                lbl_fsc_curve = f"{conc:.2f} uM ({cell})"

                series_cont.append({
                    'x': agg['Time_Window'].values, 'y': agg['m_perm'].values,
                    'se': agg['se_perm'].fillna(0).values, 'color': col_hex,
                    'dashed': (cell == control_cell), 'label': lbl_cont_curve
                })
                series_cat_v2l.append({
                    'x': agg['Time_Window'].values, 'y': agg['m_v2l'].values,
                    'se': agg['se_v2l'].fillna(0).values, 'color': col_hex,
                    'dashed': (cell == control_cell), 'label': lbl_cat_curve
                })
                series_fsc.append({
                    'x': agg['Time_Window'].values, 'y': agg['m_fsc'].values,
                    'se': agg['se_fsc'].fillna(0).values, 'color': col_hex,
                    'dashed': (cell == control_cell), 'label': lbl_fsc_curve
                })

                # Trazar líneas verticales de T_lag y T_max configurados (para célula diana o curvas con modo Manual)
                if cell == target_cell or s_tmode == "Manual":
                    if 0.0 < s_tlag < 600.0:
                        vlines_export.append((float(s_tlag), col_hex, True))
                    if 0.0 < s_tmax < 600.0:
                        vlines_export.append((float(s_tmax), col_hex, False))

            # Guardar e incrustar Cinética Continua en TIFF a 600 DPI
            kin_cont_path = out_dir / f"Figure_Category2_Kinetics_Continuous_{safe_pep}_600DPI.tiff"
            rgb_cont = render_highres_plot_png(
                kin_cont_path, series_cont, vlines_export,
                f"Cinetica Continua (Mediana FITC): {pep} [con T_lag y T_max Configurados]", "TIEMPO (S)", "PERMEABILIZACION (%)",
                (-60.0, 600.0), (-5.0, 105.0)
            )
            exported_files.append(str(kin_cont_path))
            cat2_images.append((f"Figura II.{idx_p+1}A: Cinetica de Permeabilizacion Continua ({pep}) [con T_lag y T_max Configurados]", rgb_cont))

            # Guardar e incrustar Cinética Categórica (Proporción de Afectados V2L %) en TIFF a 600 DPI
            kin_cat_path = out_dir / f"Figure_Category1_Kinetics_Categorical_V2L_{safe_pep}_600DPI.tiff"
            rgb_cat = render_highres_plot_png(
                kin_cat_path, series_cat_v2l, vlines_export,
                f"Proporcion de Eventos Afectados (V2L %): {pep} [Cut-Offs Configurados]", "TIEMPO (S)", "AFECTADOS V2L (%)",
                (-60.0, 600.0), (-5.0, 105.0)
            )
            exported_files.append(str(kin_cat_path))
            cat1_images.append((f"Figura I.{idx_p+4}: Evolucion Temporal de Proporcion de Eventos Afectados V2L (%) — {pep} [Cut-Offs Configurados]", rgb_cat))

            # Guardar e incrustar Dinámica Morfológica FSC (% Basal) en TIFF a 600 DPI
            fsc_vals_all = np.concatenate([s['y'] for s in series_fsc]) if series_fsc else np.array([100.0])
            fsc_ylim = (min(75.0, float(np.min(fsc_vals_all)) * 0.9), max(145.0, float(np.max(fsc_vals_all)) * 1.08))
            morph_path = out_dir / f"Figure_Category3_Morphology_FSC_{safe_pep}_600DPI.tiff"
            rgb_morph = render_highres_plot_png(
                morph_path, series_fsc, [(0.0, '#333333', True)],
                f"Dinamica Morfologica FSC (% Basal): {pep}", "TIEMPO (S)", "FSC NORMALIZADO (%)",
                (-60.0, 600.0), fsc_ylim
            )
            exported_files.append(str(morph_path))
            cat3_images.append((f"Figura III.{idx_p+1}: Tumefaccion Osmotica y Tamano Celular FSC (% Basal) — {pep}", rgb_morph))

        # 3. Exportar e incrustar Figura Global de Saturación en TIFF a 600 DPI
        fig_sat_path = out_dir / "Figure_Category4_Rate_Saturation_Active_600DPI.tiff"
        sat_series = []
        max_conc = max(self.df_kinetics['Concentration_uM'].max(), 1.0)
        max_vact = max(self.df_kinetics['V0_Active'].max(), 0.1)
        for idx_p, (pep, grp_p) in enumerate(self.df_kinetics.groupby('Peptide')):
            col_hex = MAIN_COLORS[idx_p % len(MAIN_COLORS)]
            for cell, grp_c in grp_p.groupby('Cell_Line'):
                agg = grp_c.groupby('Concentration_uM').agg(m=('V0_Active', 'mean'), se=('V0_Active', 'sem')).reset_index()
                sat_series.append({
                    'x': agg['Concentration_uM'].values, 'y': agg['m'].values,
                    'se': agg['se'].fillna(0).values, 'color': col_hex,
                    'dashed': (cell == control_cell), 'label': f"{pep} ({cell})"
                })
        rgb_sat = render_highres_plot_png(
            fig_sat_path, sat_series, [],
            "Curvas de Saturacion Farmacodinamica (V_active)", "CONCENTRACION (UM)", "V_ACTIVE (%/S)",
            (0.0, max_conc * 1.1), (0.0, max_vact * 1.2)
        )
        exported_files.append(str(fig_sat_path))
        cat4_images.append(("Figura IV.1: Curvas de Saturacion Alosterica sobre la Tasa Litica Activa (V_active vs Dosis)", rgb_sat))

        # ----------------------------------------------------------------------
        # 4. ESTRUCTURAR EL REPORTE PDF AUDITABLE POR CATEGORÍAS DE ANÁLISIS
        # ----------------------------------------------------------------------
        pdf_filename = out_dir / "The_Ultimate_Biostatistical_Audit_Report.pdf"

        t_prereg_r = [[
            r['Peptide'], r['Cell_Line'], f"{r['Concentration_uM']:.2f}",
            f"{r['Bin_s']}s", str(r['Filter_Scope']),
            f"<{r['V1L_Excl_Log10']:.2f}", f"<{r['V2L_Cutoff_Log10']:.2f}",
            f"{r['T_lag_s']:.0f}s", f"{r['T_max_s']:.0f}s",
            str(r['Tlag_Mode']),
            f"{int(r['Raw_Events']):,}", f"{int(r['V1L_Omitted']):,}",
            f"{int(r['V1R_Analyzed']):,} ({r['Retained_Pct']:.1f}%)"
        ] for _, r in self.df_preregistration.iterrows()]

        t_user_cuts_r = [[
            str(r['Sample_ID'])[:26], str(r['Filter_Scope']),
            f"<{r['V1L_Excl_Log10']:.2f}", f"<{r['V2L_Cutoff_Log10']:.2f}",
            f">{r['FSC_Min']:.0f}", f">{r['SSC_Min']:.0f}",
            f"{r['Bin_s']}s", f"{r['T_lag_s']:.0f}s", f"{r['T_max_s']:.0f}s",
            str(r['Tlag_Mode']), f"{r['Max_V2L_Pct']:.1f}%", f"{r['Max_Perm_Pct']:.1f}%"
        ] for _, r in getattr(self, 'df_user_cutoffs_table', pd.DataFrame()).iterrows()]

        t_cat_milestones_r = [[
            r['Peptide'], r['Cell_Line'], f"{r['Concentration_uM']:.2f}",
            f"<{r['V2L_Cutoff_Log10']:.2f}",
            f"{r['Min_V2L_Pct']:.1f}%", f"{r['Max_V2L_Pct']:.1f}%", f"{r['Min_Unaffected_V2R_Pct']:.1f}%",
            f"{r['T_Stat_Onset_s']:.0f}s", f"{r['T50_Effect_s']:.0f}s", f"{r['LT50_s']:.0f}s",
            f"{r['T_Effect_Max_s']:.0f}s", f"{r['Cohen_h_Max']:.2f}", str(r['Affectation_State'])[:18]
        ] for _, r in self.df_categorical_summary.iterrows()]

        t_cat_mech_r = [[
            r['Peptide'], r['Cell_Line'], f"{r['Concentration_uM']:.2f}",
            f"{r['AUC_V2L_Pct']:.1f}%", f"{r['V_Cat_Active']:.4f}",
            f"{r['Span_T10_T90_s']:.0f}s", f"{r['Max_Entropy_Bits']:.2f}",
            f"{r['Bimodality_BC']:.3f}", f"{r['Hazard_Rate_s']:.4f}",
            f"{r['Mode_Delta']:+.1f}", str(r['Mechanism_Class'])[:20]
        ] for _, r in self.df_categorical_summary.iterrows()]

        # Subconjunto representativo de ventanas de tiempo clave (-60, 0, 60, 120, 180, 300, 450, 585s) para el PDF
        key_tw_set = {-60.0, 0.0, 60.0, 120.0, 180.0, 300.0, 450.0, 585.0}
        if not self.df_cat_time_windows.empty:
            sub_tw_pdf = self.df_cat_time_windows[
                (self.df_cat_time_windows['Cell_Line'] == target_cell) &
                (self.df_cat_time_windows['Time_Window'].isin(key_tw_set))
            ]
            if sub_tw_pdf.empty:
                sub_tw_pdf = self.df_cat_time_windows.iloc[::4]
            t_cat_win_r = [[
                r['Peptide'], f"{r['Concentration_uM']:.2f}", f"<{r['V2L_Cutoff_Log10']:.2f}", f"{r['Time_Window']:.0f}s",
                f"{int(r['Affected_V2L_Events'])}/{int(r['Unaffected_V2R_Events'])}",
                f"{r['Mean_Affected_V2L_Pct']:.1f}+-{r['SEM_V2L_Pct']:.1f}%",
                f"{r['Mean_Unaffected_V2R_Pct']:.1f}%",
                r['Wilson_95CI_V2L'], f"{r['Ratio_V2L_to_V2R']:.2f}",
                f"{r['Cohen_h_vs_Basal']:.2f}", f"{r['q_FDR']:.4f}", str(r['Window_Status'])[:16]
            ] for _, r in sub_tw_pdf.iterrows()]
        else:
            t_cat_win_r = []

        t_cat_tests_r = [[
            str(r['Analysis_Type'])[:22], str(r['Peptide']), str(r['Comparison_Scope'])[:22],
            str(r['Metric_Focus'])[:22], str(r['Stat_Name'])[:20], str(r['Effect_Size'])[:16],
            str(r['p_value']), str(r['q_FDR']), str(r['Decision'])[:22]
        ] for _, r in self.df_cat_hypothesis_tests.iterrows()] if not self.df_cat_hypothesis_tests.empty else []

        t1_r = [[r['Peptide'], r['Cell_Line'], f"{r['Concentration_uM']:.2f}", f"<{r['V2L_Cutoff_Log10']:.2f}", f"{r['Molecules_Per_Cell']:.2e}", f"{r['Mean_V0_Early']:.4f}", f"{r['SEM_V0_Early']:.4f}"] for _, r in self.df_v0_summary.iterrows()]
        t2_r = [[r['Peptide'], r['Cell_Line'], f"{r['Concentration_uM']:.2f}", f"<{r['V2L_Cutoff_Log10']:.2f}", f"{r['Molecules_Per_Cell']:.2e}", f"{r['Mean_V0_600']:.4f}", f"{r['SEM_V0_600']:.4f}"] for _, r in self.df_v0_600_summary.iterrows()]
        t3_r = [[r['Peptide'], r['Cell_Line'], f"{r['Concentration_uM']:.2f}", f"<{r['V2L_Cutoff_Log10']:.2f}", f"{r['T_lag_s']:.0f}-{r['T_max_s']:.0f}s ({r['Tlag_Mode']})", f"{r['Molecules_Per_Cell']:.2e}", f"{r['Mean_V0_Active']:.4f}", f"{r['SEM_V0_Active']:.4f}"] for _, r in self.df_v0_active_summary.iterrows()]

        if not self.df_time_bin_overlap.empty:
            t_timebin_r = [[
                str(r['Peptide']),
                f"{r['Concentration_uM']:.2f}",
                str(r['Time_Interval']),
                f"<{r['Target_V2L_Cutoff_Log10']:.2f}" if pd.notna(r['Target_V2L_Cutoff_Log10']) else "N/A",
                f"{int(r['Target_Events_N'])}",
                f"{r['Target_Perm_Pct']:.1f}+-{r['Target_Perm_SEM']:.1f}%" if pd.notna(r['Target_Perm_Pct']) else "N/A",
                f"{r['Target_V2L_Pct']:.1f}+-{r['Target_V2L_SEM']:.1f}%" if pd.notna(r['Target_V2L_Pct']) else "N/A",
                f"<{r['Control_V2L_Cutoff_Log10']:.2f}" if pd.notna(r['Control_V2L_Cutoff_Log10']) else "N/A",
                f"{int(r['Control_Events_N'])}",
                f"{r['Control_Perm_Pct']:.1f}+-{r['Control_Perm_SEM']:.1f}%" if pd.notna(r['Control_Perm_Pct']) else "N/A",
                f"{r['Control_V2L_Pct']:.1f}+-{r['Control_V2L_SEM']:.1f}%" if pd.notna(r['Control_V2L_Pct']) else "N/A",
                f"{r['Delta_Perm_Target_vs_Ctrl_Pct']:+.1f}%" if pd.notna(r['Delta_Perm_Target_vs_Ctrl_Pct']) else "N/A",
                f"{r['Delta_V2L_Target_vs_Ctrl_Pct']:+.1f}%" if pd.notna(r['Delta_V2L_Target_vs_Ctrl_Pct']) else "N/A"
            ] for _, r in self.df_time_bin_overlap.iterrows()]
        else:
            t_timebin_r = []

        t_morph_r = [[
            r['Peptide'], r['Cell_Line'], f"{r['Concentration_uM']:.2f}",
            f"{r['End_Norm_FSC_Pct']:.1f}%", f"{0.0 if np.isnan(r['SEM_FSC']) else r['SEM_FSC']:.2f}",
            f"{r['End_Norm_SSC_Pct']:.1f}%", f"{0.0 if np.isnan(r['SEM_SSC']) else r['SEM_SSC']:.2f}"
        ] for _, r in self.df_morphology_summary.iterrows()]

        t4_r = [[
            r['Peptide'], f"{r['Dose_uM']:.2f}", f"{r['THP1_Tlag']:.0f}", f"{r['THP1_Tmax']:.0f}",
            f"{r['THP1_Vact']:.4f}", f"{r['THP1_AUC']:.1f}", f"{r['PBMC_AUC']:.1f}",
            f"{r['Delta_AUC']:.1f}", f"{r['Odds_Ratio_V2L']:.1f}x", r['CI_95'], r['p_value'], r['q_Adj']
        ] for _, r in self.df_table4.iterrows()] if not self.df_table4.empty else []

        t5_r = [[r['Peptide'], f"{r['Dose_uM']:.2f}", f"{r['T_lag']:.1f}", f"{r['T_max']:.1f}", f"{r['V_active']:.4f}", f"{r['Mean_AUC']:.2f}", f"{r['Delta_vs_Ref']:.2f}", r['CI_95'], r['p_value'], r['q_Adj']] for _, r in self.df_table5.iterrows()] if not self.df_table5.empty else []
        t6_r = [[r['Peptide'], r['Cell_Line'], r['Vmax_pct_s'], r['Kd_or_Half_uM'], r['Cooperativity_Shape'], r['Lytic_Efficiency'], r['R_Squared'], r['R2_Adj'], r['Model_Type']] for _, r in self.df_saturation.iterrows()]

        categories_report = [
            {
                'category_title': "CATEGORIA 0: DESCRIPCION DEL PRERREGISTRO Y CUT-OFFS ESPECIFICOS CONFIGURADOS POR EL USUARIO",
                'description': (
                    f"Gating Morfologico y Cut-Offs especificos por muestra: Exclusion V1L (<{v1l_range_str} log10), "
                    f"Discriminacion V2L/V2R (<{v2l_range_str} log10) y Ventana Activa [T_lag, T_max]."
                ),
                'tables': [
                    ("Tabla 0.1: Registro Auditable de Parametros de Pre-Analisis, Cut-Offs Aplicados y Retencion de Eventos (V1R vs V1L)",
                     ["Compuesto", "Celula", "Dosis", "Bin", "Alcance", "Excl V1L", "Cut V2L", "T_lag", "T_max", "Modo T", "Crudos", "Omit V1L", "V1R Analizados"],
                     t_prereg_r),
                    ("Tabla 0.2: Matriz Interactiva de Cut-Offs Especificos Configurados por el Usuario por Muestra (Panel 2)",
                     ["Muestra", "Alcance", "Excl V1L", "Cut V2L", "FSC Min", "SSC Min", "Bin", "T_lag", "T_max", "Modo T", "Max V2L%", "Max Perm%"],
                     t_user_cuts_r)
                ],
                'images': []
            },
            {
                'category_title': f"CATEGORIA I: PANEL CATEGORICO DE PROPORCIONES (AFECTADOS V2L vs NO AFECTADOS V2R) Y VENTANAS DE TIEMPO",
                'description': (
                    f"Proporciones V2L (<{v2l_range_str} log10) vs V2R (>={v2l_range_str}), % Min/Max, Tiempo a Efecto Maximo (T_max_eff), "
                    f"T50 del Efecto (relativo), LT50 (absoluto), Cochran-Armitage, Isocronas EC50(t) y pruebas exactas de afectacion."
                ),
                'tables': [
                    ("Tabla I.1: Hitos del Efecto Categorico con Cut-Off Configurado (% Afectados Min/Max, % No Afectados Min, T_onset, T50, LT50, T_efecto_max)",
                     ["Compuesto", "Celula", "Dosis", "Cut V2L", "Min V2L", "Max V2L", "Min V2R", "T_onset", "T50 Efec", "LT50 Pob", "T_ef_max", "h Cohen", "Decision Afectacion"],
                     t_cat_milestones_r),
                    ("Tabla I.2: Dinamica Mecanistica Categorica (AUC V2L, Velocidad V_cat, Span T10-T90, Entropia H, Bimodalidad BC de Sarle y Modo)",
                     ["Compuesto", "Celula", "Dosis", "AUC V2L", "V_cat", "Span 10-90", "Entropia H", "BC Sarle", "Hazard k", "Delta Modo", "Mecanismo"],
                     t_cat_mech_r),
                    (f"Tabla I.3: Proporciones entre Eventos Afectados (V2L) y No Afectados (V2R) por Ventana de Tiempo ({target_cell}, con Cut-Off Configurado)",
                     ["Compuesto", "Dosis", "Cut V2L", "Tiempo", "N V2L/V2R", "% Afect V2L", "% Sano V2R", "IC95% Wilson", "Odds V2L", "h Cohen", "q-FDR", "Estado Ventana"],
                     t_cat_win_r),
                    ("Tabla I.4: Pruebas Estadisticas para Definir Afectacion (Tendencia Cochran-Armitage, Isocronas EC50(t) y Contrastes Exactos)",
                     ["Tipo Prueba", "Compuesto", "Condicion / Tiempo", "Proporcion / EC50", "Estadistico", "Tam Efecto", "p-val", "q-FDR", "Conclusion"],
                     t_cat_tests_r)
                ],
                'images': cat1_images
            },
            {
                'category_title': "CATEGORIA II: CINETICA CONTINUA DE PERMEABILIZACION Y TASAS DE FASE ACTIVA SEGMENTADA",
                'description': "Extraccion de latencia pre-estacionaria (T_lag), inicio de meseta (T_max) y pendiente litica activa (V_active) segun los Cut-Offs configurados.",
                'tables': [
                    ("Tabla II.1: Tasa de Permeabilizacion Temprana (V0: 0-120 s) y Cut-Off Configurado",
                     ["Compuesto", "Celula", "Dosis (uM)", "Cut V2L", "Moleculas/Cel", "V0 0-120s (%/s)", "SEM"], t1_r),
                    ("Tabla II.2: Tasa de Ventana Completa (V0: 0-600 s) y Cut-Off Configurado",
                     ["Compuesto", "Celula", "Dosis (uM)", "Cut V2L", "Moleculas/Cel", "V0 0-600s (%/s)", "SEM"], t2_r),
                    ("Tabla II.3: Tasa de Fase Litica Activa Real (V_active: T_lag a T_max Configurados)",
                     ["Compuesto", "Celula", "Dosis (uM)", "Cut V2L", "Ventana [Tlag-Tmax]", "Moleculas/Cel", "V_active (%/s)", "SEM"], t3_r),
                    (f"Tabla II.4: Tabla de Datos por Bin de Tiempo — Superposicion Celula Diana ({target_cell}) vs. Celula Control ({control_cell}) con Cut-Offs",
                     ["Compuesto", "Dosis", "Bin Tiempo", f"Cut {target_cell}", f"N {target_cell}", f"Perm {target_cell}", f"V2L {target_cell}", f"Cut {control_cell}", f"N {control_cell}", f"Perm {control_cell}", f"V2L {control_cell}", "Delta Perm", "Delta V2L"],
                     t_timebin_r)
                ],
                'images': cat2_images
            },
            {
                'category_title': "CATEGORIA III: DINAMICA MORFOLOGICA CELULAR (DISPERSION FRONTAL FSC Y LATERAL SSC)",
                'description': "Monitoreo simultaneo de tumefaccion osmotica celular (FSC-A) y cambios de complejidad/granularidad citoplasmatica (SSC-A).",
                'tables': [
                    ("Tabla III.1: Cambio Relativo en Dispersion Frontal (FSC % Basal) y Lateral (SSC % Basal) al Final del Ensayo",
                     ["Compuesto", "Celula", "Dosis (uM)", "FSC Final (%)", "SEM FSC", "SSC Final (%)", "SEM SSC"], t_morph_r)
                ],
                'images': cat3_images
            },
            {
                'category_title': f"CATEGORIA IV: MODELADO FARMACODINAMICO Y COOPERATIVIDAD ALOSTERICA ({self.var_kinetic_analyzer.get()})",
                'description': "Ajuste no lineal de tasas activas y estimacion de parametros de afinidad aparente (Kd), capacidad maxima (Vmax) y cooperatividad (nH).",
                'tables': [
                    ("Tabla IV.1: Parametros Farmacodinamicos del Analizador Seleccionado e Intervalos Bootstrap 95%",
                     ["Compuesto", "Celula", "Vmax / Rate", "Kd [95% CI]", "nH / Forma", "Eficiencia", "R2", "R2 Adj", "Modelo"], t6_r)
                ],
                'images': cat4_images
            },
            {
                'category_title': f"CATEGORIA V: INFERENCIA BIOESTADISTICA, SELECTIVIDAD Y ODDS RATIO ({self.var_fdr_method.get()})",
                'description': f"Contrastes de selectividad ({target_cell} vs {control_cell}), Odds Ratio (OR) de eventos V2L y potencia frente al estado basal.",
                'tables': [
                    (f"Tabla V.1: Selectividad Celular ({target_cell} vs {control_cell}), Delta AUC y Odds Ratio Categorico (V2L)",
                     ["Compuesto", "Dosis", "T_lag", "T_max", "V_act", f"{target_cell}", f"{control_cell}", "Delta", "OR V2L", "95% CI", "p-val", "q-Adj"], t4_r),
                    (f"Tabla V.2: Contrastes Dinamicos de Potencia ({self.var_stat_mode.get()})",
                     ["Compuesto", "Dosis", "T_lag", "T_max", "V_act", "Mean AUC", "Delta", "95% CI", "p-val", "q-Adj"], t5_r)
                ],
                'images': cat5_images
            }
        ]

        subtitle_str = (
            f"Autor/Lab: {self.ent_rep_author.get().strip()} | Target: {target_cell} | "
            f"Control: {control_cell} | Canal: {self.var_perm_channel.get()} | V1L<{v1l_range_str} | V2L<{v2l_range_str} log10 | TIFF 600 DPI"
        )
        saved_pdf_path = generate_native_pdf_report(pdf_filename, self.ent_rep_title.get().strip(), subtitle_str, categories_report)
        exported_files.append(str(saved_pdf_path))

        self.txt_report_summary.delete("1.0", tk.END)
        insights_txt = self.txt_cat_insights.get("1.0", tk.END).strip() if hasattr(self, 'txt_cat_insights') else ""
        summary_lines = [
            "================================================================================",
            f" REPORTE AUDITABLE POR CATEGORÍAS CON CUT-OFFS DEL USUARIO E IMÁGENES TIFF 600 DPI ({datetime.now().strftime('%Y-%m-%d %H:%M:%S')})",
            "================================================================================",
            f"• Directorio de Exportación: {out_dir.resolve()}",
            f"• Formato de Imágenes Exportadas: TIFF 600 DPI (con compresión sin pérdida LZW)",
            f"• Canal de Permeabilidad Usado: {self.var_perm_channel.get()}",
            f"• Cut-Off 1 Exclusión V1L Configurado: < {v1l_range_str} log10",
            f"• Cut-Off 2 Discriminación V2L/V2R Configurado: < {v2l_range_str} log10",
            f"• Alcance Activo del Prefiltro: {self.var_cutoff_scope.get()}",
            f"• Analizador Cinético/Farmacodinámico: {self.var_kinetic_analyzer.get()}",
            f"• Modo de Comparación Estadística: {self.var_stat_mode.get()} ({self.var_fdr_method.get()})",
            "",
            f"--- ARCHIVOS E IMÁGENES TIFF 600 DPI EXPORTADOS ({len(exported_files)}) ---"
        ] + [f"  [OK] {fp}" for fp in exported_files] + [
            "",
            "--- CATEGORÍA 0.A: TABLA DE PRERREGISTRO Y PARÁMETROS DE PRE-ANÁLISIS ---",
            self.df_preregistration.to_string(index=False) if not self.df_preregistration.empty else "Sin datos",
            "",
            "--- CATEGORÍA 0.B: TABLA DE CUT-OFFS ESPECÍFICOS CONFIGURADOS POR EL USUARIO POR MUESTRA ---",
            self.df_user_cutoffs_table.to_string(index=False) if hasattr(self, 'df_user_cutoffs_table') and not self.df_user_cutoffs_table.empty else "Sin datos",
            "",
            "--- CATEGORÍA I.A: HITOS DEL EFECTO CATEGÓRICO (CON CUT-OFFS CONFIGURADOS POR EL USUARIO) ---",
            self.df_categorical_summary.to_string(index=False) if not self.df_categorical_summary.empty else "Sin datos",
            "",
            "--- CATEGORÍA I.C: PRUEBAS ESTADÍSTICAS DE AFECTACIÓN, COCHRAN-ARMITAGE E ISOCRONAS EC50(t) ---",
            self.df_cat_hypothesis_tests.to_string(index=False) if not self.df_cat_hypothesis_tests.empty else "Sin datos",
            "",
            insights_txt,
            "",
            "--- CATEGORÍA IV: PARÁMETROS DEL MODELO FARMACODINÁMICO (TABLA 6) ---",
            self.df_saturation.to_string(index=False) if not self.df_saturation.empty else "Sin datos",
            "",
            f"--- CATEGORÍA II.4 / V: TABLA DE DATOS POR BIN DE TIEMPO ({target_cell} vs {control_cell}) CON CUT-OFFS ---",
            self.df_time_bin_overlap.to_string(index=False) if not self.df_time_bin_overlap.empty else "Sin datos"
        ]
        self.txt_report_summary.insert(tk.END, "\n".join(summary_lines))
        self.lbl_global_status.config(text=f"Estado: Reporte PDF + TIFF 600 DPI con Cut-Offs del Usuario exportado en {out_dir.name}")

        messagebox.showinfo(
            "Exportación Completada",
            f"Se ha generado el Reporte PDF estructurado por Categorías de Análisis junto con todas las imágenes en formato TIFF a 600 DPI y las tablas CSV/TSV integrando los Cut-Offs específicos configurados por el usuario en:\n\n{out_dir.resolve()}",
            parent=self.root
        )

    # --------------------------------------------------------------------------
    # GESTIÓN DE SESIÓN Y CONFIGURACIÓN (GUARDAR / EXPORTAR / CARGAR .JSON)
    # --------------------------------------------------------------------------
    def _collect_session_state_dict(self):
        """Serializa el estado completo de la sesión de trabajo (configuración global, muestras, cut-offs y ventanas)."""
        serialized_thresholds = {}
        for (pep, cell, conc), val in self.curve_thresholds.items():
            key_str = f"{pep}||{cell}||{float(conc):.4f}"
            serialized_thresholds[key_str] = {
                'mode': str(val.get('mode', 'auto')),
                't_lag': float(val.get('t_lag', 60.0)),
                't_max': float(val.get('t_max', 300.0))
            }

        serialized_overrides = {}
        for fp, ov in self.sample_filter_overrides.items():
            f_max_v = ov.get('fsc_max')
            s_max_v = ov.get('ssc_max')
            serialized_overrides[str(fp)] = {
                'bin_size': int(ov.get('bin_size', 15)),
                'fsc_min': float(ov.get('fsc_min', 50000.0)),
                'fsc_max': float(f_max_v) if (f_max_v is not None and np.isfinite(f_max_v)) else None,
                'ssc_min': float(ov.get('ssc_min', 10000.0)),
                'ssc_max': float(s_max_v) if (s_max_v is not None and np.isfinite(s_max_v)) else None,
                'log_excl_v1l': float(ov.get('log_excl_v1l', 3.0)),
                'log_cutoff': float(ov.get('log_cutoff', 4.85)),
                'scope_label': str(ov.get('scope_label', 'Global'))
            }

        sel_peps = [self.lst_compounds.get(i) for i in self.lst_compounds.curselection()] if hasattr(self, 'lst_compounds') else []
        sel_cells = [self.lst_cells.get(i) for i in self.lst_cells.curselection()] if hasattr(self, 'lst_cells') else []
        sel_doses = [self.lst_doses.get(i) for i in self.lst_doses.curselection()] if hasattr(self, 'lst_doses') else []

        return {
            'session_version': '2.0',
            'saved_at': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'global_parameters': {
                'target_cell': self.ent_target.get().strip(),
                'control_cell': self.ent_control.get().strip(),
                'cell_density': self.ent_density.get().strip(),
                'perm_channel': self.var_perm_channel.get(),
                'fsc_min': self.ent_fsc_min.get().strip(),
                'fsc_max': self.ent_fsc_max.get().strip(),
                'ssc_min': self.ent_ssc_min.get().strip(),
                'ssc_max': self.ent_ssc_max.get().strip(),
                'bin_size': self.ent_bin_size.get().strip(),
                'log_excl_v1l': float(self.var_log_excl_v1l.get()),
                'log_cutoff': float(self.var_log_cutoff.get()),
                'cutoff_scope': self.var_cutoff_scope.get(),
                'hist_ymode': self.var_hist_ymode.get(),
                'perm_display_mode': self.var_perm_display_mode.get(),
                'scatter_norm': bool(self.var_scatter_norm.get()),
                'active_curve': self.cmb_active_curve.get().strip(),
                'cat_windows': self.ent_cat_windows.get().strip() if hasattr(self, 'ent_cat_windows') else '0-60, 60-180, 180-360, 360-600',
                'cat_vis_mode': self.var_cat_vis_mode.get() if hasattr(self, 'var_cat_vis_mode') else '',
                'kinetic_analyzer': self.var_kinetic_analyzer.get(),
                'stat_mode': self.var_stat_mode.get(),
                'fdr_method': self.var_fdr_method.get(),
                'bootstrap_n': self.ent_bootstrap_n.get().strip(),
                'report_title': self.ent_rep_title.get().strip() if hasattr(self, 'ent_rep_title') else '',
                'report_author': self.ent_rep_author.get().strip() if hasattr(self, 'ent_rep_author') else '',
                'outdir': self.ent_outdir.get().strip() if hasattr(self, 'ent_outdir') else ''
            },
            'selected_filters': {
                'compounds': sel_peps,
                'cells': sel_cells,
                'doses': sel_doses
            },
            'file_records': [dict(r) for r in self.file_records],
            'sample_filter_overrides': serialized_overrides,
            'curve_thresholds': serialized_thresholds
        }

    def _apply_session_state_dict(self, state_dict, source_label="Sesión"):
        """Restaura el estado completo de una sesión guardada y recalcula todas las vistas y tablas."""
        g = state_dict.get('global_parameters', {})

        def _set_entry(ent_widget, val_str):
            if ent_widget is not None and val_str is not None:
                ent_widget.delete(0, tk.END)
                ent_widget.insert(0, str(val_str))

        _set_entry(self.ent_target, g.get('target_cell', 'THP1'))
        _set_entry(self.ent_control, g.get('control_cell', 'Erythrocytes'))
        _set_entry(self.ent_density, g.get('cell_density', '100000'))
        _set_entry(self.ent_fsc_min, g.get('fsc_min', '50000'))
        _set_entry(self.ent_fsc_max, g.get('fsc_max', 'Sin límite'))
        _set_entry(self.ent_ssc_min, g.get('ssc_min', '10000'))
        _set_entry(self.ent_ssc_max, g.get('ssc_max', 'Sin límite'))
        _set_entry(self.ent_bin_size, g.get('bin_size', '15'))
        _set_entry(self.ent_bootstrap_n, g.get('bootstrap_n', '200'))

        if hasattr(self, 'ent_cat_windows') and g.get('cat_windows'):
            _set_entry(self.ent_cat_windows, g.get('cat_windows'))
        if hasattr(self, 'ent_rep_title') and g.get('report_title'):
            _set_entry(self.ent_rep_title, g.get('report_title'))
        if hasattr(self, 'ent_rep_author') and g.get('report_author'):
            _set_entry(self.ent_rep_author, g.get('report_author'))
        if hasattr(self, 'ent_outdir') and g.get('outdir'):
            _set_entry(self.ent_outdir, g.get('outdir'))

        if g.get('perm_channel'):
            self.var_perm_channel.set(g['perm_channel'])
        if 'log_excl_v1l' in g:
            self.var_log_excl_v1l.set(float(g['log_excl_v1l']))
        if 'log_cutoff' in g:
            self.var_log_cutoff.set(float(g['log_cutoff']))
        if g.get('cutoff_scope'):
            self.var_cutoff_scope.set(g['cutoff_scope'])
        if g.get('hist_ymode'):
            self.var_hist_ymode.set(g['hist_ymode'])
        if g.get('perm_display_mode'):
            self.var_perm_display_mode.set(g['perm_display_mode'])
        if 'scatter_norm' in g:
            self.var_scatter_norm.set(bool(g['scatter_norm']))
        if g.get('kinetic_analyzer'):
            self.var_kinetic_analyzer.set(g['kinetic_analyzer'])
        if g.get('stat_mode'):
            self.var_stat_mode.set(g['stat_mode'])
        if g.get('fdr_method'):
            self.var_fdr_method.set(g['fdr_method'])
        if hasattr(self, 'var_cat_vis_mode') and g.get('cat_vis_mode'):
            self.var_cat_vis_mode.set(g['cat_vis_mode'])

        # Restaurar registros de archivos (.fcs) si están incluidos en la sesión
        saved_records = state_dict.get('file_records', [])
        if saved_records:
            self.clear_all()
            for meta in saved_records:
                fp_p = Path(meta['Filepath'])
                if not fp_p.exists() and "demo_fcs_data" in str(fp_p):
                    generate_demo_fcs_dataset(fp_p.parent.parent.parent)
                if fp_p.exists():
                    self.file_records.append(dict(meta))
                    self.tree.insert(
                        "", "end", iid=meta['Filepath'],
                        values=(
                            meta['Filename'], meta['Compound'], meta['Cell_Line'],
                            meta['Sample_Role'], f"{float(meta['Concentration_uM']):.2f}",
                            meta['Replicate'], f"{float(meta.get('Exclusion_s', 0.0)):.1f}"
                        )
                    )

        # Restaurar umbrales de T_lag / T_max por curva
        self.curve_thresholds.clear()
        for k_str, th_val in state_dict.get('curve_thresholds', {}).items():
            parts = k_str.split("||")
            if len(parts) == 3:
                self.curve_thresholds[(parts[0], parts[1], float(parts[2]))] = {
                    'mode': str(th_val.get('mode', 'auto')),
                    't_lag': float(th_val.get('t_lag', 60.0)),
                    't_max': float(th_val.get('t_max', 300.0))
                }

        # Si hay archivos cargados, extraer caché de eventos y aplicar overrides de cada muestra
        if self.file_records:
            saved_overrides = state_dict.get('sample_filter_overrides', {})
            perm_chan = self.var_perm_channel.get().strip() or "FITC (Por defecto)"
            default_bin = int(float(self.ent_bin_size.get().strip() or 15))
            default_fsc_min = float(self.ent_fsc_min.get().strip() or 50000)
            default_fsc_max = self._parse_bound_value(self.ent_fsc_max.get())
            default_ssc_min = float(self.ent_ssc_min.get().strip() or 10000)
            default_ssc_max = self._parse_bound_value(self.ent_ssc_max.get())
            default_v1l = float(self.var_log_excl_v1l.get())
            default_v2l = float(self.var_log_cutoff.get())

            self.tube_event_cache.clear()
            self.sample_filter_overrides.clear()
            records = []
            tot_raw, tot_v1l, tot_ret = 0, 0, 0
            used_channel_name = perm_chan

            for meta in self.file_records:
                fp = meta['Filepath']
                is_cc = meta['Sample_Role'].startswith('CC')
                ov = saved_overrides.get(fp, {
                    'bin_size': default_bin,
                    'fsc_min': default_fsc_min, 'fsc_max': default_fsc_max,
                    'ssc_min': default_ssc_min, 'ssc_max': default_ssc_max,
                    'log_excl_v1l': default_v1l, 'log_cutoff': default_v2l,
                    'scope_label': 'Global'
                })
                f_max_val = float(ov['fsc_max']) if ov.get('fsc_max') is not None else default_fsc_max
                s_max_val = float(ov['ssc_max']) if ov.get('ssc_max') is not None else default_ssc_max
                ov_clean = dict(ov)
                ov_clean['fsc_max'] = f_max_val
                ov_clean['ssc_max'] = s_max_val
                try:
                    df_t, ev_cache = extract_tube_kinetics(
                        fp, is_cc=is_cc, bin_size=int(ov_clean.get('bin_size', default_bin)),
                        fsc_min=float(ov_clean.get('fsc_min', default_fsc_min)),
                        fsc_max=f_max_val,
                        ssc_min=float(ov_clean.get('ssc_min', default_ssc_min)),
                        ssc_max=s_max_val,
                        perm_channel=perm_chan,
                        log_excl_v1l=float(ov_clean.get('log_excl_v1l', default_v1l)),
                        log_cutoff=float(ov_clean.get('log_cutoff', default_v2l))
                    )
                    df_t['Applied_Scope'] = str(ov_clean.get('scope_label', 'Global'))
                    df_t['Is_CC'] = is_cc
                    self.tube_event_cache[fp] = ev_cache
                    self.sample_filter_overrides[fp] = ov_clean
                    used_channel_name = str(df_t['Perm_Channel_Used'].iloc[0])
                    tot_raw += int(df_t['Total_Raw_Events'].iloc[0])
                    tot_v1l += int(df_t['Excluded_V1L_Events'].iloc[0])
                    tot_ret += int(df_t['Retained_Events'].iloc[0])
                    records.append(df_t)
                except Exception as ex:
                    print(f"[!] Aviso restaurando {meta['Filename']}: {ex}")

            if records:
                cells_per_ml = float(self.ent_density.get().strip() or 1e5)
                pct_ret = (tot_ret / max(1, tot_raw)) * 100.0
                self.lbl_gate_stats.config(
                    text=f"Canal: {used_channel_name} | V1R Analizados: {tot_ret:,} ({pct_ret:.1f}%) | V1L Omitidos: {tot_v1l:,}"
                )
                self._rebuild_normalized_full_dataframe(records, cells_per_ml)
                self._populate_preanalysis_selectors()

                if g.get('active_curve') and g['active_curve'] in self.cmb_active_curve['values']:
                    self.cmb_active_curve.set(g['active_curve'])

                self.recompute_kinetics_dataframe()
                self.on_active_curve_selected()
                self.recompute_all_models_and_tables()
                self.lbl_global_status.config(text=f"Estado: {source_label} restaurada ({len(self.file_records)} archivos)")

    def save_session_to_file(self, fpath=None, silent=False):
        """Permite guardar y exportar la sesión completa de configuración y trabajo en un archivo .json."""
        out_dir = Path(self.ent_outdir.get().strip() or DEFAULT_OUTPUT_DIR)
        out_dir.mkdir(parents=True, exist_ok=True)
        if not fpath:
            default_name = f"sesion_citometria_{datetime.now().strftime('%Y%m%d_%H%M')}.json"
            fpath = filedialog.asksaveasfilename(
                parent=self.root,
                title="Guardar / Exportar Sesión de Configuración (.json)",
                initialdir=str(out_dir),
                initialfile=default_name,
                defaultextension=".json",
                filetypes=[("Archivos de Sesión JSON", "*.json"), ("Todos los archivos", "*.*")]
            )
        if not fpath:
            return None
        state_dict = self._collect_session_state_dict()
        Path(fpath).write_text(json.dumps(state_dict, indent=2, ensure_ascii=False), encoding="utf-8")
        self.lbl_global_status.config(text=f"Estado: Sesión exportada en {Path(fpath).name}")
        if not silent:
            messagebox.showinfo(
                "Sesión Guardada",
                f"La sesión de trabajo (parámetros globales, archivos .fcs, Cut-Offs por muestra y ventanas T_lag/T_max) se ha exportado exitosamente en:\n\n{fpath}",
                parent=self.root
            )
        return Path(fpath)

    def load_session_from_file(self, fpath=None, silent=False):
        """Permite cargar e importar una sesión de configuración .json para continuar trabajando."""
        out_dir = Path(self.ent_outdir.get().strip() or DEFAULT_OUTPUT_DIR)
        if not fpath:
            init_dir = str(out_dir) if out_dir.exists() else str(Path.cwd())
            fpath = filedialog.askopenfilename(
                parent=self.root,
                title="Cargar / Importar Sesión de Configuración (.json)",
                initialdir=init_dir,
                filetypes=[("Archivos de Sesión JSON", "*.json"), ("Todos los archivos", "*.*")]
            )
        if not fpath or not Path(fpath).exists():
            return False
        try:
            state_dict = json.loads(Path(fpath).read_text(encoding="utf-8"))
            self._apply_session_state_dict(state_dict, source_label=Path(fpath).name)
            if not silent:
                messagebox.showinfo(
                    "Sesión Cargada",
                    f"Se ha restaurado exitosamente la sesión de trabajo desde:\n\n{fpath}\n\n"
                    f"• Muestras cargadas: {len(self.file_records)}\n"
                    f"• Cut-Offs personalizados restaurados: {len(self.sample_filter_overrides)}\n"
                    f"• Ventanas T_lag/T_max restauradas: {len(self.curve_thresholds)}",
                    parent=self.root
                )
            return True
        except Exception as exc:
            if not silent:
                messagebox.showerror("Error al Cargar Sesión", f"No se pudo importar el archivo de sesión:\n{exc}", parent=self.root)
            return False

    def quick_save_session(self):
        quick_path = Path.cwd() / "last_session_config.json"
        self.save_session_to_file(fpath=quick_path, silent=False)

    def quick_load_session(self):
        quick_path = Path.cwd() / "last_session_config.json"
        if not quick_path.exists():
            messagebox.showwarning(
                "Sin Guardado Rápido",
                f"No se encontró '{quick_path.name}' en el directorio de trabajo actual.\nUsa primero '💾 Guardar Sesión'.",
                parent=self.root
            )
            return
        self.load_session_from_file(fpath=quick_path, silent=False)

    # --------------------------------------------------------------------------
    # MENÚ DE AYUDA: EXPLICACIÓN DETALLADA DE MODELOS MATEMÁTICOS Y CÁLCULOS
    # --------------------------------------------------------------------------
    def _get_mathematical_models_reference_sections(self):
        """Retorna las secciones estructuradas con todos los modelos matemáticos y cálculos estadísticos."""
        sec1_models = """================================================================================
 1. MODELOS FARMACODINÁMICOS Y CINÉTICOS NO LINEALES (PANEL 4)
================================================================================
El software ajusta mediante mínimos cuadrados no lineales (Trust-Region Reflective,
scipy.optimize.curve_fit) la tasa activa de permeabilización V_active (%/s) en función
de la concentración molar del compuesto C (µM), o bien la trayectoria temporal P(t):

1.1 MODELO ALOSTÉRICO DE HILL REGULARIZADO (Por defecto)
--------------------------------------------------------------------------------
  Ecuación:
      V_active(C) = ( V_max · C^n_H ) / ( K_d^n_H + C^n_H )

  Parámetros y Significado Físico:
  • V_max (%/s): Capacidad catalítica o velocidad máxima de lisis a saturación.
  • K_d (µM): Concentración de semi-saturación aparente (EC50 cinética donde V = 0.5·V_max).
  • n_H (adimensional, regularizado en [0.5, 6.0]): Coeficiente de cooperatividad de Hill.
      - n_H ≈ 1.0: Unión independiente tipo Langmuir (sin cooperatividad).
      - n_H > 1.5: Cooperatividad alostérica positiva (oligomerización de múltiples
        monómeros de péptido para ensamblar un poro toroidal o en barril funcional).
      - n_H < 1.0: Cooperatividad negativa o heterogeneidad de sitios de unión en membrana.
  • Eficiencia Lítica = V_max / K_d (% · s⁻¹ · µM⁻¹): Índice global de potencia catalítica.

1.2 MODELO HIPERBÓLICO DE MICHAELIS-MENTEN (n_H = 1.0 fijo)
--------------------------------------------------------------------------------
  Ecuación:
      V_active(C) = ( V_max · C ) / ( K_d + C )

  • Modela una cinética de saturación hiperbólica de primer orden sin requerimiento
    de oligomerización cooperativa en la bicapa lipídica.

1.3 MODELO SIGMOIDAL ASIMÉTRICO DE GOMPERTZ
--------------------------------------------------------------------------------
  Ecuación:
      V_active(C) = V_max · exp( -exp( -k_G · (C - K_d) ) )

  • Captura respuestas dosis-efecto asimétricas donde el inicio de la nucleación
    del poro ocurre abruptamente superada una concentración crítica (punto de inflexión K_d),
    seguido de una aproximación gradual hacia la asíntota V_max.

1.4 MODELO DE WEIBULL / EXPONENCIAL ESTIRADO (FRACTAL)
--------------------------------------------------------------------------------
  Ecuación:
      V_active(C) = V_max · [ 1 - exp( -( C / K_d )^m ) ]

  • Derivado de la teoría de percolación y cinética de falla de membranas.
  • El exponente de forma m cuantifica la distribución de umbrales de resistencia
    mecánica de la población celular ante la tensión interfacial inducida por el péptido.

1.5 MODELO BI-EXPONENCIAL SECUENCIAL EN TIEMPO (A Dosis Máxima)
--------------------------------------------------------------------------------
  Ecuación:
      P(t) = P_ads · (1 - exp(-k_ads · t)) + P_lys · (1 - exp(-k_lys · max(0, t - T_lag)))

  • Desacopla explícitamente dos regímenes temporales sobre la curva P(t):
      1) Fase rápida pre-lítica de adsorción electrostática / perturbación superficial
         con amplitud P_ads (%) y constante de velocidad k_ads (s⁻¹).
      2) Fase secundaria de colapso de permeabilidad / apertura de poros tras el tiempo
         de latencia crítica T_lag (s), con amplitud P_lys (%) y tasa k_lys (s⁻¹).
"""

        sec2_kinetics = """================================================================================
 2. GATING MORFOLÓGICO, SUBPOBLACIONES (V1L/V2L/V2R) Y EXTRACCIÓN CINÉTICA
================================================================================

2.1 JERARQUÍA DE FILTRADO Y SUBPOBLACIONES EN EL HISTOGRAMA LOG10(FITC)
--------------------------------------------------------------------------------
Cada evento celular i registrado en el archivo .fcs posee dispersión frontal (FSC_i),
dispersión lateral (SSC_i), fluorescencia de Calceína (FITC_i) y tiempo (t_i):

  1) Gate Morfológico (Células Íntegras en Tamaño/Granularidad):
         FSC_min < FSC_i <= FSC_max    y    SSC_min < SSC_i <= SSC_max

  2) Cut-Off 1 — Exclusión de Detritos / Células Sin Carga (V1L vs V1R):
         • V1L (Omitidos):  log10(FITC_i) < CutOff_V1L   (células muertas de origen o sin calceína)
         • V1R (Analizados): log10(FITC_i) >= CutOff_V1L  (población viable total N_V1R)

  3) Cut-Off 2 — Clasificación Categórica de Afectación (V2L vs V2R dentro de V1R):
         • V2L (Eventos Afectados / Permeabilizados):
               CutOff_V1L <= log10(FITC_i) < CutOff_V2L
         • V2R (Eventos No Afectados / Retención Intacta):
               log10(FITC_i) >= CutOff_V2L

  Calibración Automática desde el Control Basal CC (-60 a 0 s):
         CutOff_V2L_auto = Percentil_05( log10(FITC_CC) ) ≈ Media_CC - 1.645 · SD_CC
  Garantiza una especificidad basal >= 95% (<= 5% de falsos positivos en t < 0 s).

2.2 NORMALIZACIÓN CONTINUA POR VENTANA TEMPORAL (BINNING)
--------------------------------------------------------------------------------
El eje temporal se divide en ventanas de ancho Δt (por defecto 15 s). Para cada ventana t_k:
  • Retención_%(t_k) = [ Mediana_FITC(t_k) / Mediana_FITC_Basal_CC ] × 100
  • Permeabilización_Continua_%(t_k) = clip( 100 - Retención_%(t_k), 0, 100 )
  • Proporción Categórica Afectada:
        %V2L(t_k) = [ N_V2L(t_k) / N_V1R(t_k) ] × 100
        %V2R(t_k) = 100 - %V2L(t_k)

2.3 EXTRACCIÓN AUTOMÁTICA Y MANUAL DE T_lag, T_max Y V_active
--------------------------------------------------------------------------------
  • T_lag (Tiempo de Latencia Pre-Estacionaria, s):
    Instante donde la señal supera de forma sostenida el umbral de ruido basal
    (Umbral = max(5.0%, Media_basal + 3·SD_basal)) e inicia la fase de pendiente máxima.
  • T_max (Inicio de Meseta de Saturación, s):
    Primer instante posterior a T_lag donde se alcanza el 90% de la permeabilización
    máxima del tubo (o el fin de la ventana de máxima pendiente).
  • V_active (%/s):
    Pendiente de regresión lineal ordinaria calculada exclusivamente dentro de la
    ventana de fase lítica activa [T_lag, T_max]. Elimina el sesgo de dilución temporal
    que sufren las tasas de ventana fija (V0_0-120s y V0_0-600s) cuando T_lag > 120 s.
  • AUC Normalizada (%):
        AUC_norm = [ Integral_Trapezoidal_{0}^{600} P(t) dt / (600 s × 100%) ] × 100
  • Estequiometría Interfacial (Moléculas de Péptido / Célula):
        Mol/Célula = ( Concentración_µM × 10⁻⁹ mol/mL × N_Avogadro ) / ( Células/mL )
"""

        sec3_categorical = """================================================================================
 3. HITOS TEMPORALES (T50 EFECTO vs LT50), ENTROPÍA, HAZARD Y BIMODALIDAD
================================================================================

3.1 DISTINCIÓN ENTRE T50 DEL EFECTO Y LT50 POBLACIONAL
--------------------------------------------------------------------------------
  • % Afectados Mínimo (V2L_min) y Máximo (V2L_max):
    Proporción mínima y máxima de eventos en la región V2L durante la ventana post-adición,
    junto con el tiempo exacto para lograr el efecto máximo (T_efecto_máx).
  • T50 del Efecto (s):
    Tiempo requerido para alcanzar el 50% de la amplitud neta propia de esa curva:
        Meta_50_Efecto = V2L_basal + 0.50 · (V2L_max - V2L_basal)
    (Siempre está definido para cualquier compuesto activo, incluso si V2L_max < 50%).
  • LT50 Poblacional (s):
    Tiempo requerido para que exactamente el 50% absoluto de toda la población celular
    viable cruce hacia la región afectada (%V2L(t) >= 50%).
  • Rango Dinámico de Transición (ΔT_{10→90}) e Índice de Sincronía (Sync_Index):
        ΔT_{10→90} = T90_efecto - T10_efecto
        Sync_Index = (V2L_max - V2L_basal) / ΔT_{10→90}   (%/s)
    Valores altos indican un colapso poblacional altamente sincronizado (tipo avalancha).

3.2 ENTROPÍA DE SHANNON H(t) Y TASA DE RIESGO INSTANTÁNEO (HAZARD RATE λ)
--------------------------------------------------------------------------------
  • Entropía Binaria de Shannon H(t) [bits]:
        H(t) = -p(t)·log2(p(t)) - (1 - p(t))·log2(1 - p(t)),   donde p(t) = %V2L(t)/100
    Mide la heterogeneidad binaria de la población celular: H = 0 bits cuando todas las
    células están sanas (p=0) o todas lisadas (p=1), y alcanza su máximo H = 1.000 bit
    exactamente en el punto de máxima coexistencia bimodal (p = 0.50, en t = LT50).
  • Tasa de Riesgo Instantáneo de Conversión (Nelson-Aalen Hazard Rate, λ en s⁻¹):
        S(t) = %V2R(t) / 100   (Fracción de supervivencia intacta)
        Λ(t) = -ln( S(t) )     (Riesgo acumulado de Nelson-Aalen)
        λ_media = Pendiente de regresión lineal de Λ(t) vs t en la fase activa.

3.3 COEFICIENTE DE BIMODALIDAD DE SARLE (BC) Y CLASIFICACIÓN DEL MECANISMO
--------------------------------------------------------------------------------
  • Coeficiente de Bimodalidad de Sarle sobre log10(FITC):
        BC = ( Asimetría_Skewness² + 1 ) / Curtosis_Pearson
  • Índice de Divergencia Continua-Categórica:
        Δ_modo = | Permeabilización_Continua_% - %V2L |
  • Clasificación Automática del Mecanismo de Membrana:
      - "Todo-o-Nada (Poro Estable / Bimodal)" si BC > 0.555 o Δ_modo > 15%:
        Coexisten células 100% vacías y células 100% llenas (típico de poros estables).
      - "Relevo Graduado (Fuga Homogénea / Detergente)" si BC <= 0.555 y Δ_modo <= 15%:
        Toda la población pierde fluorescencia de manera uniforme y gradual.
"""

        sec4_stats = """================================================================================
 4. PRUEBAS ESTADÍSTICAS PARA DEFINIR LA AFECTACIÓN Y SELECTIVIDAD
================================================================================

4.1 PRUEBA EXACTA DE FISHER, CHI-CUADRADO (χ²) Y TAMAÑO DEL EFECTO h DE COHEN
--------------------------------------------------------------------------------
En cada ventana de tiempo y en el pico de efecto se construye la tabla de contingencia 2×2
de conteos celulares reales:
                  | Eventos Afectados (V2L) | Eventos No Afectados (V2R) |
    Tratamiento   |         a = V2L_trat    |          b = V2R_trat      |
    Basal / Ref   |         c = V2L_basal   |          d = V2R_basal     |

  • Estadístico χ² de Pearson (con corrección de continuidad de Yates) y Prueba Exacta
    de Fisher bilateral (p_Fisher).
  • Tamaño del Efecto para Proporciones h de Cohen (independiente del tamaño muestral N):
        h = 2 · arcsin( sqrt(p_trat) ) - 2 · arcsin( sqrt(p_basal) )
    Interpretación biofísica de |h|:
      - |h| < 0.20 : Efecto trivial / Ruido instrumental
      - 0.20 <= |h| < 0.50 : Afectación Leve / Incipiente
      - 0.50 <= |h| < 0.80 : Afectación Moderada
      - |h| >= 0.80 : Afectación Severa / Lítica

4.2 TIEMPO DE INICIO ESTADÍSTICO (T_onset_stat) Y REGLA DE DECISIÓN DE AFECTACIÓN
--------------------------------------------------------------------------------
  • T_onset_stat (s): Primera ventana temporal post-adición donde la proporción de
    afectados supera simultáneamente significancia estadística (p_Fisher < 0.01) y
    relevancia biológica (h de Cohen >= 0.20).
  • Clasificación Categórica de Afectación:
      - NO AFECTADO (Inerte): p_Fisher >= 0.01 o h_max < 0.20
      - AFECTACIÓN LEVE (Incipiente): p_Fisher < 0.01 y 0.20 <= h_max < 0.50
      - AFECTACIÓN MODERADA: p_Fisher < 0.01 y 0.50 <= h_max < 0.80
      - AFECTACIÓN SEVERA (Lítica): p_Fisher < 0.01 y (h_max >= 0.80 o V2L_max >= 50%)

4.3 INTERVALOS DE CONFIANZA DE WILSON 95%, ODDS RATIO Y TENDENCIA COCHRAN-ARMITAGE
--------------------------------------------------------------------------------
  • IC 95% de Wilson Score para la proporción p = V2L / N:
        IC_95% = [ p + z²/(2N) ± z · sqrt( p(1-p)/N + z²/(4N²) ) ] / (1 + z²/N),  z = 1.96
    Evita límites negativos o superiores a 100% cuando p es cercano a 0% o 100%.
  • Razón de Momios (Odds Ratio, OR) con corrección de Haldane-Anscombe (+0.5):
        OR = [ (a + 0.5) · (d + 0.5) ] / [ (b + 0.5) · (c + 0.5) ]
  • Prueba de Tendencia Dosis-Respuesta de Cochran-Armitage (Z_trend):
    Evalúa en cada ventana temporal si la proporción de eventos afectados V2L crece
    monotónicamente con la concentración del compuesto (C_1 < C_2 < ... < C_k).

4.4 INFERENCIA PARAMÉTRICA ENTRE RÉPLICAS Y CORRECCIÓN POR MULTIPLICIDAD (FDR)
--------------------------------------------------------------------------------
  • Prueba t de dos muestras de Welch (varianzas desiguales, grados de libertad de
    Welch-Satterthwaite) para contrastar AUC Normalizada y V_active entre líneas
    celulares (Target vs Control), entre compuestos o entre dosis sucesivas.
  • Corrección de Valores p por Comparaciones Múltiples:
      - Benjamini-Hochberg (FDR, control de tasa de falsos descubrimientos q-value).
      - Holm-Bonferroni (control fuerte de tasa de error familiar FWER).
"""
        return [
            ("📐 1. Modelos Farmacodinámicos y Cinéticos", sec1_models),
            ("🔬 2. Gating (V1L/V1R), Cut-Off V2L/V2R y Cinética", sec2_kinetics),
            ("⏱️ 3. Hitos (T50 vs LT50), Entropía, Hazard y Bimodalidad", sec3_categorical),
            ("📊 4. Pruebas Estadísticas de Afectación y Selectividad", sec4_stats),
        ]

    def open_mathematical_models_help_dialog(self, initial_tab=0):
        """Abre la ventana interactiva de Ayuda explicando todos los modelos matemáticos y cálculos."""
        dlg = tk.Toplevel(self.root)
        dlg.title("❓ Manual de Modelos Matemáticos, Ecuaciones Biofísicas y Cálculos Estadísticos")
        dlg.geometry("980x680")
        dlg.transient(self.root)

        hdr = tk.Frame(dlg, bg="#002B49", pady=6, padx=12)
        hdr.pack(fill="x")
        tk.Label(
            hdr,
            text="📐 REFERENCIA MATEMÁTICA Y BIOESTADÍSTICA DE LA ESTACIÓN DE TRABAJO DE CITOMETRÍA",
            bg="#002B49", fg="white", font=("DejaVu Sans", 10, "bold")
        ).pack(side="left")
        tk.Button(
            hdr, text="📄 Exportar Manual en Markdown (.md)",
            bg="#117864", fg="white", font=("DejaVu Sans", 8, "bold"),
            command=self.export_markdown_guide_file
        ).pack(side="right", padx=4)

        nb = ttk.Notebook(dlg)
        nb.pack(fill="both", expand=True, padx=10, pady=8)

        sections = self._get_mathematical_models_reference_sections()
        for tab_title, content_txt in sections:
            frm = ttk.Frame(nb, padding=6)
            nb.add(frm, text=f" {tab_title} ")
            txt = tk.Text(frm, wrap="word", font=("Consolas", 9), bg="#FCFDFE", fg="#1C2833", padx=10, pady=8)
            vsb = ttk.Scrollbar(frm, orient="vertical", command=txt.yview)
            txt.configure(yscrollcommand=vsb.set)
            txt.pack(side="left", fill="both", expand=True)
            vsb.pack(side="right", fill="y")
            txt.insert("1.0", content_txt)
            txt.configure(state="disabled")

        if 0 <= initial_tab < len(sections):
            nb.select(initial_tab)

        return dlg

    # --------------------------------------------------------------------------
    # GUÍA DE USO DEL SOFTWARE (INTERACTIVA + EXPORTACIÓN A MARKDOWN)
    # --------------------------------------------------------------------------
    def _get_user_guide_sections(self):
        """Retorna las secciones paso a paso de la Guía de Uso del Software."""
        g1 = """================================================================================
 PASO 1: CARGA DE DATOS (.FCS), CONFIGURACIÓN DE CANALES Y METADATOS (PESTAÑA 1)
================================================================================
1. Nomenclatura Automática o Edición Manual:
   • Puedes cargar archivos estándar de citometría (.fcs versiones 2.0, 3.0 o 3.1) haciendo
     clic en '📂 Subir Archivos .fcs'.
   • Si deseas explorar todas las funciones de inmediato con datos realistas validados,
     haz clic en '🧪 Generar y Cargar Set Demo (.fcs)'.
   • El software autodetecta del nombre del archivo:
       - Rol de muestra: Control Basal sin péptido ('CC', -60 a 0 s) vs Experimental ('EXP', 0 a 600 s).
       - Nombre del Compuesto / Péptido (ej. Pep1, Pep2).
       - Línea Celular (ej. THP1 vs Erythrocytes).
       - Concentración en µM (ej. 0.62, 1.25, 2.50, 5.00, 10.00 µM) y Réplica (R1, R2, R3).
   • Si tus archivos tienen otros nombres, selecciona una o varias filas en la tabla y usa
     '✏️ Editar Metadatos de Seleccionados' para asignar Compuesto, Línea, Dosis o Rol en lote.

2. Configuración Experimental Global:
   • Define la 'Línea Celular Diana (Target)' (ej. THP1) y la 'Línea Control Sano' (ej. Erythrocytes).
   • Ajusta la 'Densidad Celular (cél/mL)' para calcular la relación estequiométrica exacta
     de Moléculas de Péptido por Célula.
   • Selecciona el 'Canal de Permeabilidad' (FITC, FITC-A, FL1-A, etc.) y pulsa el botón verde:
     '🚀 EJECUTAR PRE-ANÁLISIS E IR A GRÁFICAS DINÁMICAS'.
"""

        g2 = """================================================================================
 PASO 2: PRE-ANÁLISIS DINÁMICO Y TABLA INTERACTIVA DE CUT-OFFS (PESTAÑA 2)
================================================================================
El Panel 2 está dividido en 4 cajas de control a la izquierda y 4 lienzos interactivos a la derecha:

1. Caja 1 — Prefiltros Morfológicos (FSC/SSC) y Doble Cut-Off en el Histograma:
   • Línea Morada Punteada (Cut-Off 1: Exclusión V1L):
     Todo evento con log10(FITC) inferior a este umbral se clasifica como V1L (detrito o célula
     no cargada) y se excluye del denominador para no sesgar el análisis.
   • Línea Roja Discontinua (Cut-Off 2: Discriminación V2L vs V2R):
     Divide las células viables analizadas (V1R) en:
       - V2L (Izquierda, sombreado rojo): Eventos AFECTADOS / Permeabilizados.
       - V2R (Derecha): Eventos NO AFECTADOS / Intactos.
   • Puedes arrastrar ambas líneas directamente con el ratón sobre el Histograma superior,
     o usar el botón 'Auto CC' para calibrar el Cut-Off V2L al percentil 5% del Control Basal.
   • Puedes aplicar los umbrales de manera '🌐 Global (a TODAS las muestras)' o de forma
     '🎯 Discriminada (solo a las muestras seleccionadas)'.

2. Caja 2 y Caja 3 — Selección de Curvas y Ventanas Arrastrables T_lag / T_max:
   • Filtra qué Compuestos, Líneas Celulares y Dosis deseas superponer en tiempo real.
   • Arrastra con el ratón la línea azul punteada (T_lag: inicio de fase lítica) y la línea
     roja punto-raya (T_max: inicio de meseta) directamente sobre la Gráfica 1 de Permeabilización,
     o usa los deslizadores. La pendiente activa V_active se recalcula al instante.

3. Caja 4 — Tabla Interactiva de Cut-Offs por Muestra (¡NUEVO!):
   • Muestra en formato tabular cada muestra experimental con su Alcance (Global, Discriminado
     o Individual), V1L_Log, CutOff V2L_Log, FSC_Min, SSC_Min, T_lag, T_max, Modo y %V2L Máx.
   • CLIC SIMPLE sobre una fila: Selecciona automáticamente esa curva en los gráficos,
     sincroniza el histograma y actualiza los deslizadores.
   • DOBLE CLIC sobre cualquier celda numérica (V1L_Log, V2L_Log, FSC_Min, SSC_Min, T_lag o T_max):
     Abre un editor en línea directamente sobre la celda. Al pulsar Enter, actualiza el Cut-Off
     o la ventana temporal específicamente para esa muestra ('Individual') y recalcula todo.
   • También puedes pulsar '✏️ Editar Muestra Sel.' para abrir el editor modal completo o
     '📋 Exportar Cut-Offs (.csv)' para guardar la tabla de umbrales.
"""

        g3 = """================================================================================
 PASO 3: PANEL CATEGÓRICO DE AFECTADOS vs NO AFECTADOS Y VENTANAS (PESTAÑA 3)
================================================================================
1. Configuración de Ventanas de Tiempo Personalizadas:
   • Ingresa rangos de tiempo separados por comas en 'Ventanas de Tiempo (s)'
     (por defecto: '0-60, 60-180, 180-360, 360-600') y haz clic en '🔄 Actualizar Panel Categórico'.

2. Cuatro Modos de Visualización Interactiva en el Lienzo Superior:
   • 1) Barras Apiladas: % Afectados (V2L rojo) vs % No Afectados (V2R verde) en cada ventana
        de tiempo para la muestra seleccionada, indicando Odds Ratio (OR) y p de Fisher.
   • 2) Cinética Multi-Dosis: Curvas de evolución de % Afectados (V2L) para todas las dosis
        con marcadores automáticos de T50 del Efecto (diamante) y Tiempo de Efecto Máximo (estrella).
   • 3) Mapa de Calor (Heatmap): Matriz Compuesto × Dosis × Ventana de Tiempo con código de color
        de severidad e indicación de significancia estadística (*, **, ***).
   • 4) Isocronas Dosis-Respuesta: Curvas de % Afectados vs Concentración para cada ventana
        temporal con estimación de EC50(t) dependiente del tiempo.

3. Cuatro Sub-Pestañas de Tablas e Interpretación Biofísica:
   • Sub-pestaña 1: Hitos del Efecto (% Afectados Mín/Máx, T_efecto_máx, T50 del efecto, LT50,
     ΔT10-90, h de Cohen, OR, RR, p_Fisher y clasificación de Afectación).
   • Sub-pestaña 2: Proporciones detalladas por Ventana de Tiempo con IC 95% de Wilson.
   • Sub-pestaña 3: Pruebas de Tendencia de Cochran-Armitage e Isocronas EC50(t).
   • Sub-pestaña 4: Interpretación Biofísica Automática y sugerencias mecanísticas.
"""

        g4 = """================================================================================
 PASO 4: MODELOS CINÉTICOS, SESIONES (.JSON) Y REPORTE FINAL 600 DPI (PESTAÑAS 4 Y 5)
================================================================================
1. Pestaña 4 — Selección de Modelos Matemáticos y Vista Previa Interactiva:
   • Elige entre los 5 analizadores matemáticos:
       - Allosteric Hill (Regularizado nH <= 6.0)
       - Michaelis-Menten Hiperbólico (nH = 1.0)
       - Gompertz Sigmoidal Asimétrico
       - Weibull / Exponencial Estirado
       - Bi-Exponencial Secuencial en Tiempo (k_ads + k_lysis)
   • A la derecha verás la 'Vista Previa Interactiva del Modelo Matemático Usado', que muestra
     en tiempo real la ecuación algebraica del modelo activo, los puntos experimentales con
     barras de error SEM, las curvas suavizadas ajustadas, las líneas guía de K_d (o T_lag)
     y el resumen de parámetros (V_max, K_d, n_H, R²).
   • Revisa las Tablas 0 a 6 a la izquierda (Prerregistro, V0 temprana, V0 completa, V_active,
     Resumen Categórico, Selectividad Celular, Contrastes Dinámicos y Parámetros del Modelo).

2. Guardar / Exportar y Cargar Sesión de Trabajo (.json):
   • En cualquier momento puedes usar el menú superior '💾 Sesión y Configuración' o los botones
     de la barra superior ('💾 Guardar Sesión' [Ctrl+S] / '📂 Cargar Sesión' [Ctrl+O]).
   • El archivo .json guarda absolutamente toda tu configuración: parámetros globales, lista de
     archivos .fcs y sus metadatos, Cut-Offs individuales de cada muestra en la Tabla del Panel 2,
     y ventanas T_lag/T_max manuales o automáticas, permitiéndote retomar el trabajo exactamente
     donde lo dejaste o compartir la sesión auditable con otro investigador.

3. Pestaña 5 — Generación del Reporte Final Auditable (PDF + PNG 600 DPI + CSV/TSV):
   • Define el título del reporte, autor/laboratorio y carpeta de salida.
   • Haz clic en '💾 GENERAR Y EXPORTAR REPORTE FINAL EN PDF E IMÁGENES 600 DPI'.
   • El software exporta automáticamente todas las tablas en CSV/TSV, las figuras de alta
     resolución a 600 DPI (incluyendo la lámina de 4 paneles del Análisis Categórico) y el
     Reporte PDF completo estructurado por Categorías de Análisis (Categorías 0 a V).
"""
        return [
            ("📁 Paso 1: Carga de Datos y Metadatos", g1),
            ("🔍 Paso 2: Pre-Análisis y Tabla de Cut-Offs", g2),
            ("🧬 Paso 3: Panel Categórico y Ventanas", g3),
            ("📊 Paso 4: Modelos, Sesiones (.json) y Reporte", g4),
        ]

    def open_user_guide_dialog(self):
        """Abre la ventana interactiva con la Guía Paso a Paso de Uso del Software."""
        dlg = tk.Toplevel(self.root)
        dlg.title("📖 Guía Paso a Paso de Uso del Software de Cinética de Citometría")
        dlg.geometry("960x660")
        dlg.transient(self.root)

        hdr = tk.Frame(dlg, bg="#002B49", pady=6, padx=12)
        hdr.pack(fill="x")
        tk.Label(
            hdr,
            text="📖 GUÍA PRÁCTICA DE USO DEL SOFTWARE — FLUJO DE TRABAJO COMPLETO (PESTAÑAS 1 A 5)",
            bg="#002B49", fg="white", font=("DejaVu Sans", 10, "bold")
        ).pack(side="left")
        tk.Button(
            hdr, text="📄 Exportar Guía Completa en Markdown (.md)",
            bg="#117864", fg="white", font=("DejaVu Sans", 8, "bold"),
            command=self.export_markdown_guide_file
        ).pack(side="right", padx=4)

        nb = ttk.Notebook(dlg)
        nb.pack(fill="both", expand=True, padx=10, pady=8)

        for tab_title, content_txt in self._get_user_guide_sections():
            frm = ttk.Frame(nb, padding=6)
            nb.add(frm, text=f" {tab_title} ")
            txt = tk.Text(frm, wrap="word", font=("Consolas", 9), bg="#FCFDFE", fg="#1C2833", padx=10, pady=8)
            vsb = ttk.Scrollbar(frm, orient="vertical", command=txt.yview)
            txt.configure(yscrollcommand=vsb.set)
            txt.pack(side="left", fill="both", expand=True)
            vsb.pack(side="right", fill="y")
            txt.insert("1.0", content_txt)
            txt.configure(state="disabled")

        return dlg

    def export_markdown_guide_file(self, fpath=None, silent=False):
        """Exporta un documento Markdown completo con la Guía de Uso y el Manual de Modelos Matemáticos."""
        default_md = Path.cwd() / "GUIA_DE_USO_Y_MODELOS_MATEMATICOS.md"
        if not fpath:
            fpath = filedialog.asksaveasfilename(
                parent=self.root,
                title="Exportar Guía de Uso y Manual Matemático (.md)",
                initialdir=str(Path.cwd()),
                initialfile=default_md.name,
                defaultextension=".md",
                filetypes=[("Markdown files", "*.md"), ("All files", "*.*")]
            )
        if not fpath:
            return None

        lines = [
            "# 🧬 Guía de Uso y Manual de Modelos Matemáticos — Estación de Trabajo de Cinética de Citometría de Flujo (`.fcs`)",
            "",
            f"*Generado automáticamente el {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*",
            "",
            "---",
            "",
            "## Parte I: Guía Paso a Paso de Uso del Software",
            ""
        ]
        for title, body in self._get_user_guide_sections():
            lines.append(f"### {title}")
            lines.append("```text")
            lines.append(body.strip())
            lines.append("```")
            lines.append("")

        lines.append("---")
        lines.append("")
        lines.append("## Parte II: Modelos Matemáticos, Ecuaciones Biofísicas y Cálculos Estadísticos")
        lines.append("")
        for title, body in self._get_mathematical_models_reference_sections():
            lines.append(f"### {title}")
            lines.append("```text")
            lines.append(body.strip())
            lines.append("```")
            lines.append("")

        Path(fpath).write_text("\n".join(lines), encoding="utf-8")
        if not silent:
            messagebox.showinfo(
                "Guía Exportada",
                f"La Guía de Uso y Manual de Modelos Matemáticos se ha guardado en:\n\n{fpath}",
                parent=self.root
            )
        return Path(fpath)


