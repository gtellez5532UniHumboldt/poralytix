"""
Pruebas unitarias de modelos no lineales de ajuste dosis-respuesta.
"""
import unittest
import numpy as np
from cytoworkstation.core.models import (
    hill_model,
    michaelis_menten_model,
    gompertz_dose_model,
    fit_dose_response_model
)


class TestModels(unittest.TestCase):
    def test_hill_model_eval(self):
        c = np.array([0.0, 1.0, 2.0, 10.0])
        res = hill_model(c, vmax=1.0, kd=2.0, nh=2.0)
        self.assertAlmostEqual(res[0], 0.0, places=4)
        self.assertAlmostEqual(res[2], 0.5, places=4)  # en C = Kd, V = 0.5*Vmax

    def test_michaelis_menten_eval(self):
        c = np.array([0.0, 2.0])
        res = michaelis_menten_model(c, vmax=1.0, km=2.0)
        self.assertAlmostEqual(res[0], 0.0, places=4)
        self.assertAlmostEqual(res[1], 0.5, places=4)

    def test_fit_hill_curve(self):
        doses = np.array([0.5, 1.0, 2.0, 4.0, 8.0])
        # Datos simulados con Hill Vmax=0.8, Kd=2.0, nH=2.5
        true_y = hill_model(doses, 0.8, 2.0, 2.5)
        fit_res = fit_dose_response_model(doses, true_y, analyzer_name="Allosteric Hill")
        self.assertIsNotNone(fit_res)
        self.assertGreater(fit_res['r2_grp'], 0.98)
        self.assertAlmostEqual(fit_res['vmax'], 0.8, delta=0.05)
        self.assertAlmostEqual(fit_res['kd'], 2.0, delta=0.1)


if __name__ == "__main__":
    unittest.main()
