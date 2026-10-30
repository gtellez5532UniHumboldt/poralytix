"""
Pruebas Unitarias del Sistema de Internacionalización (i18n) y Selección de Idioma.
Verifica traducciones completas ES/EN, conmutación dinámica y persistencia en sesión.
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

    def test_04_translation_formatting(self):
        text_es = t('editor_title', lang='es').format("Sample_1")
        self.assertIn("Sample_1", text_es)
        self.assertIn("Editar", text_es)

        text_en = t('editor_title', lang='en').format("Sample_1")
        self.assertIn("Sample_1", text_en)
        self.assertIn("Edit", text_en)

    def test_05_fallback_behavior(self):
        # Clave inexistente devuelve la misma clave
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

    def test_07_session_persistence_of_language(self):
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
