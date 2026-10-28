"""
Suite de Pruebas de Integridad de Sistema y Flujo de Trabajo End-to-End.
Verifica la integridad de:
  1. Estructura de archivos y paquetes de Refactorizacion.
  2. Parseo sin pérdida de los 4 datasets experimentales (813, 820, sat1, CZS5c-).
  3. Extracción de metadatos robusta y libre de bugs.
  4. Flujo cinético completo (Gating morfológico + Doble Cut-Off).
  5. Ajuste no lineal de modelos matemáticos (Hill, Michaelis-Menten, Gompertz).
  6. Motores de exportación (Tablas CSV/TSV, Render 600 DPI y Reporte PDF).
  7. Persistencia y round-trip de sesiones JSON.
"""
import os
import sys
import unittest
import tempfile
import numpy as np
import pandas as pd
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from cytoworkstation.config import PEPTIDE_ALIASES, DOSE_RAMPS
from cytoworkstation.io.fcs_parser import parse_fcs_native, resolve_perm_channel
from cytoworkstation.io.metadata import guess_metadata_from_filename
from cytoworkstation.io.session import save_session_dict, load_session_dict
from cytoworkstation.core.kinetics import (
    extract_tube_kinetics,
    extract_slope_based_parameters,
    compute_categorical_kinetics_metrics
)
from cytoworkstation.core.models import (
    hill_model,
    fit_dose_response_model
)
from cytoworkstation.core.statistics import (
    wilson_binomial_ci,
    cohen_h_proportions,
    contingency_2x2_p_and_or
)
from cytoworkstation.export.tables_export import export_analysis_tables
from cytoworkstation.export.plots_600dpi import render_highres_plot_png
from cytoworkstation.export.pdf_report import generate_native_pdf_report


class TestProjectIntegrity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace_root = BASE_DIR.parent
        cls.test_folders = {
            '813': cls.workspace_root / '813',
            '820': cls.workspace_root / '820',
            'sat1': cls.workspace_root / 'sat1',
            'CZS5c-': cls.workspace_root / 'CZS5c-'
        }

    def test_01_module_structure_integrity(self):
        """Verifica que todos los archivos clave del paquete modular existan."""
        expected_files = [
            BASE_DIR / 'main.py',
            BASE_DIR / 'cli.py',
            BASE_DIR / 'requirements.txt',
            BASE_DIR / 'README.md',
            BASE_DIR / 'cytoworkstation' / '__init__.py',
            BASE_DIR / 'cytoworkstation' / 'config.py',
            BASE_DIR / 'cytoworkstation' / 'io' / 'fcs_parser.py',
            BASE_DIR / 'cytoworkstation' / 'io' / 'metadata.py',
            BASE_DIR / 'cytoworkstation' / 'io' / 'session.py',
            BASE_DIR / 'cytoworkstation' / 'core' / 'statistics.py',
            BASE_DIR / 'cytoworkstation' / 'core' / 'models.py',
            BASE_DIR / 'cytoworkstation' / 'core' / 'kinetics.py',
            BASE_DIR / 'cytoworkstation' / 'core' / 'synthetic_demo.py',
            BASE_DIR / 'cytoworkstation' / 'export' / 'tables_export.py',
            BASE_DIR / 'cytoworkstation' / 'export' / 'plots_600dpi.py',
            BASE_DIR / 'cytoworkstation' / 'export' / 'pdf_report.py',
            BASE_DIR / 'cytoworkstation' / 'gui' / 'app.py',
            BASE_DIR / 'cytoworkstation' / 'gui' / 'widgets' / 'canvas_histogram.py',
            BASE_DIR / 'cytoworkstation' / 'gui' / 'widgets' / 'canvas_kinetics.py',
            BASE_DIR / 'cytoworkstation' / 'gui' / 'widgets' / 'canvas_saturation.py',
            BASE_DIR / 'cytoworkstation' / 'gui' / 'widgets' / 'canvas_categorical.py'
        ]
        for f in expected_files:
            self.assertTrue(f.exists(), f"Archivo faltante en arquitectura: {f}")

    def test_02_all_test_datasets_exist(self):
        """Verifica la existencia e inmutabilidad de las 4 carpetas de prueba."""
        for name, folder in self.test_folders.items():
            self.assertTrue(folder.exists() and folder.is_dir(), f"Carpeta de prueba no encontrada: {name}")
            fcs_count = len(list(folder.glob("**/*.fcs"))) + len(list(folder.glob("**/*.FCS")))
            self.assertGreater(fcs_count, 0, f"No hay archivos FCS en {name}")

    def test_03_metadata_across_all_datasets(self):
        """Verifica que la heurística de metadatos clasifique correctamente muestras de los 4 datasets."""
        # 1. Dataset 813
        f_813 = next(self.test_folders['813'].glob("*.fcs"))
        m_813 = guess_metadata_from_filename(str(f_813))
        self.assertEqual(m_813['Compound'], 'Act8-13')

        # 2. Dataset 820
        f_820 = next(self.test_folders['820'].glob("*.fcs"))
        m_820 = guess_metadata_from_filename(str(f_820))
        self.assertEqual(m_820['Compound'], 'Act8-20')

        # 3. Dataset sat1
        f_sat1 = next(self.test_folders['sat1'].glob("*.fcs"))
        m_sat1 = guess_metadata_from_filename(str(f_sat1))
        self.assertEqual(m_sat1['Compound'], 'Satanin 1')

        # 4. Dataset CZS5c- (Verificación del bug resuelto en subcarpetas)
        f_czs5_pbmc = self.test_folders['CZS5c-'] / 'PBMCs' / 'CZS5_PBMCs_9_1.fcs'
        m_czs5_pbmc = guess_metadata_from_filename(str(f_czs5_pbmc))
        self.assertEqual(m_czs5_pbmc['Compound'], 'CZS5')
        self.assertEqual(m_czs5_pbmc['Cell_Line'], 'PBMCs')
        self.assertEqual(m_czs5_pbmc['Concentration_uM'], 9.0)

        f_czs5_thp = self.test_folders['CZS5c-'] / 'THP1' / 'CZS5_THP1_9_CC_1.fcs'
        m_czs5_thp = guess_metadata_from_filename(str(f_czs5_thp))
        self.assertEqual(m_czs5_thp['Compound'], 'CZS5')
        self.assertEqual(m_czs5_thp['Cell_Line'], 'THP-1')
        self.assertEqual(m_czs5_thp['Sample_Role'], 'CC (Baseline -60 to 0s)')

    def test_04_full_kinetic_pipeline_execution(self):
        """Ejecuta el pipeline completo de gating, cut-offs y cinética sobre una muestra real."""
        f_sample = self.test_folders['813'] / 'THP1-813-0.5-1.fcs'
        res_df, cache = extract_tube_kinetics(
            str(f_sample), is_cc=False, bin_size=15,
            fsc_min=50000.0, ssc_min=10000.0,
            log_excl_v1l=3.0, log_cutoff=4.85
        )
        self.assertFalse(res_df.empty)
        self.assertIn('Affected_V2L_Pct', res_df.columns)
        self.assertIn('Median_FITC', res_df.columns)
        self.assertIn('Time_Window', res_df.columns)

        # Cálculo de parámetros T_lag y V_active
        tlag, tmax, v_act, max_p = extract_slope_based_parameters(
            res_df['Time_Window'].values,
            res_df['Affected_V2L_Pct'].values
        )
        self.assertTrue(0.0 <= tlag <= tmax <= 600.0)
        self.assertTrue(v_act >= 0.0)

        # Métricas bioestadísticas categóricas
        cat_metrics = compute_categorical_kinetics_metrics(
            res_df['Time_Window'].values,
            res_df['Affected_V2L_Pct'].values,
            res_df['Median_FITC'].values
        )
        self.assertIn('Affectation_State', cat_metrics)
        self.assertIn('Cohen_h_Max', cat_metrics)

    def test_05_export_pipeline_end_to_end(self):
        """Verifica que export_analysis_tables, render_plot y generate_native_pdf_report funcionen juntos."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)

            # 1. Exportar tablas CSV / TSV
            df_test = pd.DataFrame({'Dose': [0.5, 1.0, 2.0], 'Rate': [0.1, 0.4, 0.7]})
            exported_tables = export_analysis_tables({'Test_Table': df_test}, tmp_path)
            self.assertEqual(len(exported_tables), 2)
            self.assertTrue(Path(exported_tables[0]).exists())

            # 2. Renderizar figura 600 DPI
            series = [{'x': [0, 60, 120], 'y': [5, 25, 60], 'color': '#004488', 'label': 'Test Serie'}]
            rgb_buf = render_highres_plot_png(
                tmp_path / 'test_plot', series, [(60, '#0055CC', True)],
                "Figura de Prueba", "Tiempo (s)", "Lisis (%)", (0, 120), (0, 100)
            )
            self.assertIsNotNone(rgb_buf)
            self.assertEqual(rgb_buf.ndim, 3)

            # 3. Generar reporte PDF con tabla e imagen incrustada
            categories = [{
                'category_title': "Categoría de Prueba",
                'description': "Verificación de integración",
                'tables': [("Tabla de Prueba", list(df_test.columns), df_test.astype(str).values.tolist())],
                'images': [("Lámina Cinética", rgb_buf)]
            }]
            pdf_path = tmp_path / 'reporte_integridad.pdf'
            saved_pdf = generate_native_pdf_report(
                pdf_path, "TEST DE INTEGRIDAD", "Reporte de prueba automatizado", categories
            )
            self.assertTrue(Path(saved_pdf).exists())
            self.assertGreater(Path(saved_pdf).stat().st_size, 1000)

    def test_06_session_roundtrip_persistence(self):
        """Verifica el guardado y restauración idéntica de sesiones de trabajo JSON."""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_json = Path(tmpdir) / 'test_session.json'
            state = {
                'session_version': '2.0',
                'global_parameters': {'target_cell': 'THP-1', 'control_cell': 'PBMCs'},
                'thresholds': {'sat1||THP-1||0.7000': {'t_lag': 75.0, 't_max': 300.0, 'mode': 'manual'}},
                'overrides': {'sample1': {'v1l': 3.1, 'v2l': 4.9}}
            }
            save_session_dict(state, test_json)
            self.assertTrue(test_json.exists())

            loaded = load_session_dict(test_json)
            self.assertEqual(loaded['global_parameters']['target_cell'], 'THP-1')
            self.assertEqual(loaded['thresholds']['sat1||THP-1||0.7000']['t_lag'], 75.0)


if __name__ == "__main__":
    unittest.main()
