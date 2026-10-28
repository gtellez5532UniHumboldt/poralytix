"""
Pruebas unitarias para la extracción y deducción de metadatos experimentales.
Valida la corrección del bug en subcarpetas y el reconocimiento del nuevo compuesto CZS5.
"""
import unittest
from pathlib import Path
from cytoworkstation.io.metadata import guess_metadata_from_filename


class TestMetadataExtraction(unittest.TestCase):
    def test_czs5_nested_pbmcs(self):
        fpath = r"D:\Usuario\Documents\Analisis_Citometria\analizador\CZS5c-\PBMCs\CZS5_PBMCs_9_1.fcs"
        meta = guess_metadata_from_filename(fpath)
        self.assertEqual(meta['Compound'], "CZS5")
        self.assertEqual(meta['Cell_Line'], "PBMCs")
        self.assertEqual(meta['Concentration_uM'], 9.0)
        self.assertEqual(meta['Sample_Role'], "Experimental (0 to 600s)")
        self.assertEqual(meta['Replicate'], "R1")

    def test_czs5_nested_thp1_cc(self):
        fpath = r"D:\Usuario\Documents\Analisis_Citometria\analizador\CZS5c-\THP1\CZS5_THP1_9_CC_2.fcs"
        meta = guess_metadata_from_filename(fpath)
        self.assertEqual(meta['Compound'], "CZS5")
        self.assertEqual(meta['Cell_Line'], "THP-1")
        self.assertEqual(meta['Concentration_uM'], 9.0)
        self.assertEqual(meta['Sample_Role'], "CC (Baseline -60 to 0s)")
        self.assertEqual(meta['Replicate'], "R2")

    def test_813_standard(self):
        fpath = r"D:\Usuario\Documents\Analisis_Citometria\analizador\813\THP1-813-0.5-1.fcs"
        meta = guess_metadata_from_filename(fpath)
        self.assertEqual(meta['Compound'], "Act8-13")
        self.assertEqual(meta['Cell_Line'], "THP-1")
        self.assertEqual(meta['Concentration_uM'], 0.5)
        self.assertEqual(meta['Sample_Role'], "Experimental (0 to 600s)")
        self.assertEqual(meta['Replicate'], "R1")

    def test_sat1_standard(self):
        fpath = r"D:\Usuario\Documents\Analisis_Citometria\analizador\sat1\PBMCs-CC-sat1-5.7-3.fcs"
        meta = guess_metadata_from_filename(fpath)
        self.assertEqual(meta['Compound'], "Satanin 1")
        self.assertEqual(meta['Cell_Line'], "PBMCs")
        self.assertEqual(meta['Concentration_uM'], 5.7)
        self.assertEqual(meta['Sample_Role'], "CC (Baseline -60 to 0s)")
        self.assertEqual(meta['Replicate'], "R3")

    def test_820_standard(self):
        fpath = r"D:\Usuario\Documents\Analisis_Citometria\analizador\820\THP1-820-0.54-1.fcs"
        meta = guess_metadata_from_filename(fpath)
        self.assertEqual(meta['Compound'], "Act8-20")
        self.assertEqual(meta['Cell_Line'], "THP-1")
        self.assertEqual(meta['Concentration_uM'], 0.54)
        self.assertEqual(meta['Sample_Role'], "Experimental (0 to 600s)")
        self.assertEqual(meta['Replicate'], "R1")


if __name__ == "__main__":
    unittest.main()
