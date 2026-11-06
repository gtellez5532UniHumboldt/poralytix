"""
Pruebas Unitarias del Sistema de Internacionalización (i18n) y Selección de Idioma.
Verifica traducciones completas ES/EN, conmutación dinámica en todos los paneles,
submenús/sub-pestañas, encabezados de tablas, comboboxes y persistencia en sesión.
"""
import unittest
import os
import sys
import tkinter as tk
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from cytoworkstation.i18n import (
    t, get_current_language, set_current_language,
    SUPPORTED_LANGUAGES, TRANSLATIONS
)
from cytoworkstation.gui.app import FlowCytometryKineticWorkstation


class TestI18nSystem(unittest.TestCase):
    def setUp(self):
        set_current_language('es')

    def tearDown(self):
        set_current_language('es')

    def test_01_supported_languages(self):
        self.assertIn('es', SUPPORTED_LANGUAGES)
        self.assertIn('en', SUPPORTED_LANGUAGES)
        self.assertEqual(get_current_language(), 'es')

    def test_02_set_language(self):
        set_current_language('en')
        self.assertEqual(get_current_language(), 'en')
        self.assertEqual(t('status_ready'), "Status: Ready (Native Engine Active)")

        set_current_language('es')
        self.assertEqual(get_current_language(), 'es')
        self.assertEqual(t('status_ready'), "Estado: Listo (Motor Nativo Activo)")

    def test_03_translations_symmetry(self):
        """Verifica que todas las claves en español tengan su contraparte en inglés y viceversa."""
        es_keys = set(TRANSLATIONS['es'].keys())
        en_keys = set(TRANSLATIONS['en'].keys())

        missing_in_en = es_keys - en_keys
        missing_in_es = en_keys - es_keys

        self.assertEqual(missing_in_en, set(), f"Claves faltantes en inglés: {missing_in_en}")
        self.assertEqual(missing_in_es, set(), f"Claves faltantes en español: {missing_in_es}")
        self.assertGreaterEqual(len(es_keys), 230)

    def test_04_translation_formatting(self):
        text_es = t('editor_title', lang='es').format("Sample_1")
        self.assertIn("Sample_1", text_es)
        self.assertIn("Editar", text_es)

        text_en = t('editor_title', lang='en').format("Sample_1")
        self.assertIn("Sample_1", text_en)
        self.assertIn("Edit", text_en)

    def test_05_fallback_behavior(self):
        self.assertEqual(t("non_existent_key_12345"), "non_existent_key_12345")

    def test_06_gui_language_switching(self):
        """Verifica que la interfaz gráfica cambie títulos, botones, pestañas y tablas dinámicamente."""
        root = tk.Tk()
        root.withdraw()
        try:
            app = FlowCytometryKineticWorkstation(root)

            # Estado inicial en Español
            self.assertEqual(app.current_language, 'es')
            self.assertIn("Estación de Trabajo", root.title())
            self.assertIn("1. Carga de Datos", app.notebook.tab(app.tab_ingest, "text"))
            self.assertEqual(app.tree.heading("Filename", "text"), "Archivo .fcs")

            # Conmutar a Inglés
            app.set_language('en')
            self.assertEqual(app.current_language, 'en')
            self.assertIn("Cell Permeabilization", root.title())
            self.assertIn("1. Data Ingestion", app.notebook.tab(app.tab_ingest, "text"))
            self.assertEqual(app.tree.heading("Filename", "text"), ".fcs File")
            self.assertEqual(app.btn_top_save.cget("text"), "💾 Save Session")

            # Conmutar de vuelta a Español
            app.set_language('es')
            self.assertEqual(app.current_language, 'es')
            self.assertIn("Estación de Trabajo", root.title())
            self.assertIn("1. Carga de Datos", app.notebook.tab(app.tab_ingest, "text"))
            self.assertEqual(app.tree.heading("Filename", "text"), "Archivo .fcs")
            self.assertEqual(app.btn_top_save.cget("text"), "💾 Guardar Sesión")

        finally:
            root.destroy()

    def test_07_all_panels_and_submenus_language_switching(self):
        """Verifica exhaustivamente que los contenidos de los paneles 2, 3, 4, 5 y submenús se traduzcan."""
        root = tk.Tk()
        root.withdraw()
        try:
            app = FlowCytometryKineticWorkstation(root)

            # --- Conmutar a INGLÉS ---
            app.set_language('en')

            # 1. Pestañas principales
            self.assertEqual(app.notebook.tab(app.tab_preanalysis, "text"), t("tab_2_title", lang='en'))
            self.assertEqual(app.notebook.tab(app.tab_catpanel, "text"), t("tab_3_title", lang='en'))
            self.assertEqual(app.notebook.tab(app.tab_prereport, "text"), t("tab_4_title", lang='en'))
            self.assertEqual(app.notebook.tab(app.tab_finalreport, "text"), t("tab_5_title", lang='en'))

            # 2. Panel 2 (Pre-análisis)
            self.assertEqual(app.lblfrm_tab2_gate.cget("text"), t("tab2_gate_title", lang='en'))
            self.assertEqual(app.lbl_tab2_perm_chan.cget("text"), t("tab2_perm_channel", lang='en'))
            self.assertEqual(app.lbl_tab2_fsc_min.cget("text"), t("tab2_fsc_min", lang='en'))
            self.assertEqual(app.btn_tab2_apply_all.cget("text"), t("tab2_btn_apply_all", lang='en'))
            self.assertEqual(app.btn_tab2_apply_sel.cget("text"), t("tab2_btn_apply_sel", lang='en'))
            self.assertEqual(app.lblfrm_tab2_sel.cget("text"), t("tab2_sel_title", lang='en'))
            self.assertEqual(app.lbl_tab2_compounds.cget("text"), t("tab2_compounds", lang='en'))
            self.assertEqual(app.lbl_tab2_cells.cget("text"), t("tab2_cell_lines", lang='en'))
            self.assertEqual(app.lbl_tab2_doses.cget("text"), t("tab2_doses", lang='en'))
            self.assertEqual(app.btn_tab2_sel_all.cget("text"), t("tab2_btn_sel_all", lang='en'))
            self.assertEqual(app.btn_tab2_refresh_plots.cget("text"), t("tab2_btn_refresh", lang='en'))
            self.assertEqual(app.lblfrm_tab2_thresh.cget("text"), t("tab2_thresh_title", lang='en'))
            self.assertEqual(app.lbl_tab2_active_curve.cget("text"), t("tab2_active_curve", lang='en'))
            self.assertEqual(app.rad_tab2_thresh_auto.cget("text"), t("tab2_radio_auto", lang='en'))
            self.assertEqual(app.btn_tab2_auto_all.cget("text"), t("tab2_btn_auto_all", lang='en'))
            self.assertEqual(app.lblfrm_tab2_cut_tbl.cget("text"), t("tab2_cut_tbl_title", lang='en'))
            self.assertEqual(app.btn_tab2_edit_sample.cget("text"), t("tab2_btn_edit_sample", lang='en'))
            self.assertEqual(app.btn_tab2_export_cutoffs.cget("text"), t("tab2_btn_export_csv", lang='en'))
            self.assertEqual(app.btn_tab2_goto_cat.cget("text"), t("tab2_btn_go_cat", lang='en'))
            # Encabezados de tabla de cut-offs
            self.assertEqual(app.tree_sample_cutoffs.heading("Muestra", "text"), t("tab2_col_sample", lang='en'))
            self.assertEqual(app.tree_sample_cutoffs.heading("Scope", "text"), t("tab2_col_scope", lang='en'))

            # 3. Panel 3 (Panel Categórico y sus 4 sub-pestañas / submenús)
            self.assertEqual(app.lblfrm_tab3_ctrl.cget("text"), t("tab3_ctrl_title", lang='en'))
            self.assertEqual(app.lbl_tab3_pep.cget("text"), t("tab3_compound", lang='en'))
            self.assertEqual(app.lbl_tab3_cell.cget("text"), t("tab3_cell_line", lang='en'))
            self.assertEqual(app.btn_tab3_recompute.cget("text"), t("tab3_btn_recalc", lang='en'))
            self.assertEqual(app.btn_tab3_goto_final.cget("text"), t("tab3_btn_export_report", lang='en'))
            # 4 Sub-pestañas del Notebook categórico
            self.assertEqual(app.cat_nb.tab(0, "text"), t("tab3_subtab_milestones", lang='en'))
            self.assertEqual(app.cat_nb.tab(1, "text"), t("tab3_subtab_windows", lang='en'))
            self.assertEqual(app.cat_nb.tab(2, "text"), t("tab3_subtab_tests", lang='en'))
            self.assertEqual(app.cat_nb.tab(3, "text"), t("tab3_subtab_insights", lang='en'))

            # 4. Panel 4 (Pre-Reporte, Sub-notebook de 9 tablas y Previsualización)
            self.assertEqual(app.lblfrm_tab4_ctrl.cget("text"), t("tab4_ctrl_title", lang='en'))
            self.assertEqual(app.lbl_tab4_analyzer.cget("text"), t("tab4_analyzer", lang='en'))
            self.assertEqual(app.lbl_tab4_stat_mode.cget("text"), t("tab4_stat_mode", lang='en'))
            self.assertEqual(app.btn_tab4_recompute.cget("text"), t("tab4_btn_recalc", lang='en'))
            self.assertEqual(app.btn_tab4_explain.cget("text"), t("tab4_btn_explain", lang='en'))
            # 9 Sub-pestañas del tbl_notebook
            self.assertEqual(app.tbl_notebook.tab(0, "text"), t("tab4_subtab_t0", lang='en'))
            self.assertEqual(app.tbl_notebook.tab(1, "text"), t("tab4_subtab_t1", lang='en'))
            self.assertEqual(app.tbl_notebook.tab(2, "text"), t("tab4_subtab_t2", lang='en'))
            self.assertEqual(app.tbl_notebook.tab(3, "text"), t("tab4_subtab_t3", lang='en'))
            self.assertEqual(app.tbl_notebook.tab(4, "text"), t("tab4_subtab_tcat", lang='en'))
            self.assertEqual(app.tbl_notebook.tab(5, "text"), t("tab4_subtab_t4", lang='en'))
            self.assertEqual(app.tbl_notebook.tab(6, "text"), t("tab4_subtab_t5", lang='en'))
            self.assertEqual(app.tbl_notebook.tab(7, "text"), t("tab4_subtab_t6", lang='en'))
            self.assertEqual(app.tbl_notebook.tab(8, "text"), t("tab4_subtab_ttimebin", lang='en'))
            # Panel de previsualización y ecuaciones
            self.assertEqual(app.lblfrm_tab4_sat.cget("text"), t("tab4_sat_title", lang='en'))
            self.assertEqual(app.lbl_tab4_prev_group.cget("text"), t("tab4_sample_group", lang='en'))
            self.assertEqual(app.lbl_tab4_prev_mode.cget("text"), t("tab4_view", lang='en'))
            self.assertEqual(app.lblfrm_tab4_eq.cget("text"), t("tab4_eq_title", lang='en'))
            self.assertEqual(app.tree_model_eq_r2.heading("Grupo", "text"), t("tab4_col_group", lang='en'))
            self.assertEqual(app.tree_model_eq_r2.heading("Modelo", "text"), t("tab4_col_model", lang='en'))
            self.assertEqual(app.tree_model_eq_r2.heading("Ecuacion_Evaluada", "text"), t("tab4_col_eq", lang='en'))

            # 5. Panel 5 (Reporte Final)
            self.assertEqual(app.lblfrm_tab5_top.cget("text"), t("tab5_top_title", lang='en'))
            self.assertEqual(app.lbl_tab5_outdir.cget("text"), t("tab5_outdir", lang='en'))
            self.assertEqual(app.btn_tab5_browse.cget("text"), t("tab5_btn_browse", lang='en'))
            self.assertEqual(app.lbl_tab5_rep_title.cget("text"), t("tab5_rep_title", lang='en'))
            self.assertEqual(app.btn_tab5_export_bundle.cget("text"), t("tab5_btn_export_bundle", lang='en'))
            self.assertEqual(app.lblfrm_tab5_preview.cget("text"), t("tab5_preview_title", lang='en'))

            # --- Conmutar de vuelta a ESPAÑOL ---
            app.set_language('es')

            # 1. Pestañas principales
            self.assertEqual(app.notebook.tab(app.tab_preanalysis, "text"), t("tab_2_title", lang='es'))
            self.assertEqual(app.notebook.tab(app.tab_catpanel, "text"), t("tab_3_title", lang='es'))
            self.assertEqual(app.notebook.tab(app.tab_prereport, "text"), t("tab_4_title", lang='es'))
            self.assertEqual(app.notebook.tab(app.tab_finalreport, "text"), t("tab_5_title", lang='es'))

            # 2. Panel 2
            self.assertEqual(app.lblfrm_tab2_gate.cget("text"), t("tab2_gate_title", lang='es'))
            self.assertEqual(app.lbl_tab2_perm_chan.cget("text"), t("tab2_perm_channel", lang='es'))
            self.assertEqual(app.lbl_tab2_fsc_min.cget("text"), t("tab2_fsc_min", lang='es'))
            self.assertEqual(app.btn_tab2_apply_all.cget("text"), t("tab2_btn_apply_all", lang='es'))
            self.assertEqual(app.tree_sample_cutoffs.heading("Muestra", "text"), t("tab2_col_sample", lang='es'))

            # 3. Panel 3 y sus submenús
            self.assertEqual(app.lblfrm_tab3_ctrl.cget("text"), t("tab3_ctrl_title", lang='es'))
            self.assertEqual(app.cat_nb.tab(0, "text"), t("tab3_subtab_milestones", lang='es'))
            self.assertEqual(app.cat_nb.tab(1, "text"), t("tab3_subtab_windows", lang='es'))
            self.assertEqual(app.cat_nb.tab(2, "text"), t("tab3_subtab_tests", lang='es'))
            self.assertEqual(app.cat_nb.tab(3, "text"), t("tab3_subtab_insights", lang='es'))

            # 4. Panel 4 y sus submenús
            self.assertEqual(app.lblfrm_tab4_ctrl.cget("text"), t("tab4_ctrl_title", lang='es'))
            self.assertEqual(app.tbl_notebook.tab(0, "text"), t("tab4_subtab_t0", lang='es'))
            self.assertEqual(app.tbl_notebook.tab(1, "text"), t("tab4_subtab_t1", lang='es'))
            self.assertEqual(app.tbl_notebook.tab(4, "text"), t("tab4_subtab_tcat", lang='es'))
            self.assertEqual(app.tbl_notebook.tab(8, "text"), t("tab4_subtab_ttimebin", lang='es'))
            self.assertEqual(app.tree_model_eq_r2.heading("Grupo", "text"), t("tab4_col_group", lang='es'))
            self.assertEqual(app.tree_model_eq_r2.heading("Modelo", "text"), t("tab4_col_model", lang='es'))

            # 5. Panel 5
            self.assertEqual(app.lblfrm_tab5_top.cget("text"), t("tab5_top_title", lang='es'))
            self.assertEqual(app.lbl_tab5_outdir.cget("text"), t("tab5_outdir", lang='es'))
            self.assertEqual(app.btn_tab5_browse.cget("text"), t("tab5_btn_browse", lang='es'))
            self.assertEqual(app.btn_tab5_export_bundle.cget("text"), t("tab5_btn_export_bundle", lang='es'))

        finally:
            root.destroy()

    def test_08_session_persistence_of_language(self):
        """Verifica que el idioma se serialice y restaure correctamente en sesiones .json."""
        root = tk.Tk()
        root.withdraw()
        try:
            app = FlowCytometryKineticWorkstation(root)

            # Cambiar a inglés y serializar estado
            app.set_language('en')
            state = app._collect_session_state_dict()
            self.assertEqual(state.get('language'), 'en')

            # Reinicializar en español y restaurar estado
            app.set_language('es')
            self.assertEqual(app.current_language, 'es')

            app._apply_session_state_dict(state)
            self.assertEqual(app.current_language, 'en')
            self.assertIn("Cell Permeabilization", root.title())

        finally:
            root.destroy()


if __name__ == "__main__":
    unittest.main()
