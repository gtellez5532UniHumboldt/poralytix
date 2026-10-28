"""
Pruebas de lectura e integridad del parser binario de archivos .fcs.
"""
import unittest
from pathlib import Path
from cytoworkstation.io.fcs_parser import parse_fcs_native, resolve_perm_channel


class TestFCSParser(unittest.TestCase):
    def test_parse_real_fcs_813(self):
        fpath = r"D:\Usuario\Documents\Analisis_Citometria\analizador\813\THP1-813-0.5-1.fcs"
        meta, df = parse_fcs_native(fpath)
        self.assertGreater(len(df), 1000)
        self.assertIn('$TOT', meta)
        self.assertIn('FSC-A', df.columns)
        self.assertIn('SSC-A', df.columns)

        perm_col = resolve_perm_channel(df.columns)
        self.assertIsNotNone(perm_col)
        self.assertIn('FL1', perm_col.upper())

    def test_parse_real_fcs_czs5(self):
        fpath = r"D:\Usuario\Documents\Analisis_Citometria\analizador\CZS5c-\THP1\CZS5_THP1_9_1.fcs"
        meta, df = parse_fcs_native(fpath)
        self.assertGreater(len(df), 1000)
        self.assertIn('Time', df.columns)


if __name__ == "__main__":
    unittest.main()
