"""
Modelos Farmacodinámicos y Cinéticos de Permeabilización de Membrana.
Ajustes no lineales con scipy.optimize.curve_fit:
  - Hill Alostérico Regularizado (nH <= 6.0)
  - Michaelis-Menten Hiperbólico (nH = 1.0)
  - Gompertz Sigmoidal Asimétrico
  - Weibull / Exponencial Estirado (Fractal)
  - Bi-Exponencial Secuencial en Tiempo (Fase Adsorción + Fase Lisis)
"""
import numpy as np
from scipy.optimize import curve_fit


def hill_model(c, vmax, kd, nh):
    """
    Modelo alostérico de Hill regularizado:
    V(C) = (Vmax * C^nH) / (Kd^nH + C^nH)
    """
    c_safe = np.maximum(0.0, np.asarray(c, dtype=float))
    return (vmax * (c_safe ** nh)) / ((kd ** nh) + (c_safe ** nh) + 1e-12)


def michaelis_menten_model(c, vmax, km):
    """
    Modelo de Michaelis-Menten (Hill con nH = 1.0 fijo):
    V(C) = (Vmax * C) / (Km + C)
    """
    c_safe = np.maximum(0.0, np.asarray(c, dtype=float))
    return (vmax * c_safe) / (km + c_safe + 1e-12)


def gompertz_dose_model(c, vmax, c_infl, k_g):
    """
    Modelo sigmoidal asimétrico de Gompertz:
    V(C) = Vmax * exp(-exp(-k_G * (C - c_infl)))
    """
    c_safe = np.maximum(0.0, np.asarray(c, dtype=float))
    return vmax * np.exp(-np.exp(-k_g * (c_safe - c_infl)))


def weibull_dose_model(c, vmax, lam, k_w):
    """
    Modelo de Weibull / exponencial estirado:
    V(C) = Vmax * (1 - exp(-(C / lambda)^k_w))
    """
    c_safe = np.maximum(0.0, np.asarray(c, dtype=float))
    return vmax * (1.0 - np.exp(-((c_safe / (lam + 1e-12)) ** k_w)))


def biexponential_time_model(t, p_ads, k_ads, p_lysis, k_lysis, t_lag):
    """
    Modelo bi-exponencial secuencial en el tiempo:
    P(t) = P_ads * (1 - exp(-k_ads * t)) + P_lysis * (1 - exp(-k_lysis * max(0, t - t_lag)))
    """
    t_arr = np.maximum(0.0, np.asarray(t, dtype=float))
    term_ads = p_ads * (1.0 - np.exp(-k_ads * t_arr))
    t_active = np.maximum(0.0, t_arr - t_lag)
    term_lysis = p_lysis * (1.0 - np.exp(-k_lysis * t_active))
    return term_ads + term_lysis


def fit_dose_response_model(x_doses, y_rates, analyzer_name="Allosteric Hill", raw_x=None, raw_y=None):
    """
    Ajusta la curva de velocidad frente a la dosis según el modelo seleccionado.
    
    Parameters
    ----------
    x_doses : array-like
        Concentraciones ensayadas (µM)
    y_rates : array-like
        Tasas medias observadas (%/s o V_active)
    analyzer_name : str
        Nombre del modelo ('Allosteric Hill', 'Michaelis-Menten', 'Gompertz', 'Weibull')
    raw_x, raw_y : array-like, optional
        Puntos experimentales de todas las réplicas individuales
        
    Returns
    -------
    dict
        Parámetros óptimos, predicciones, métricas R² y R² ajustado
    """
    x = np.asarray(x_doses, dtype=float)
    y = np.asarray(y_rates, dtype=float)
    n_doses = len(x)
    if n_doses < 2:
        return None

    if raw_x is None:
        raw_x = x
    if raw_y is None:
        raw_y = y

    max_v0 = max(0.01, float(np.max(y)))
    max_c = max(1.0, float(np.max(x)))
    x_grid = np.linspace(0.0, max_c * 1.15, 120)

    if "Michaelis-Menten" in analyzer_name:
        p0 = [max_v0 * 1.1, float(np.median(x))]
        bounds = ([max_v0 * 0.4, 0.01], [max_v0 * 3.0, max_c * 3.0])
        popt, _ = curve_fit(michaelis_menten_model, x, y, p0=p0, bounds=bounds, maxfev=15000)
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
        bounds = ([max_v0 * 0.4, 0.01, 0.1], [max_v0 * 3.0, max_c * 2.5, 10.0])
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
        bounds = ([max_v0 * 0.4, 0.01, 0.5], [max_v0 * 3.0, max_c * 3.0, 8.0])
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
        # Allosteric Hill por defecto (Regularizado nH en [0.5, 6.0])
        p0 = [max_v0 * 1.1, float(np.median(x)), 2.0]
        bounds = ([max_v0 * 0.4, 0.01, 0.5], [max_v0 * 3.0, max_c * 2.5, 6.0])
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
    r2_rep = max(0.0, 1.0 - float(np.sum((raw_y - raw_pred) ** 2)) / ss_tot_raw) if ss_tot_raw > 1e-12 else r2_grp

    return {
        'model_lbl': model_lbl,
        'vmax': float(vmax_f),
        'kd': float(kd_f),
        'nh': float(nh_f),
        'popt': popt,
        'eq_eval': eq_eval,
        'x_fit': x_grid,
        'y_fit': y_smooth,
        'y_half': y_half_val,
        'r2_grp': float(r2_grp),
        'r2_rep': float(r2_rep),
        'r2_adj': float(r2_adj),
        'rmse': rmse_val
    }
