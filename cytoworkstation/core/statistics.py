"""
Módulo Bioestadístico: Intervalos de Wilson, Cochran-Armitage, Tablas de Contingencia 2x2,
Odds Ratio, Riesgo Relativo, Tamaño del efecto h de Cohen y ajuste de p-valores.
"""
import math
import numpy as np
import scipy.stats as stats


def wilson_binomial_ci(k_aff, n_tot, z=1.95996):
    """
    Calcula el intervalo de confianza binomial exacto de Wilson al 95% (en porcentaje 0-100%).
    
    Parameters
    ----------
    k_aff : float or int
        Número de eventos afectados (positivos)
    n_tot : float or int
        Total de eventos analizados en el denominador
    z : float
        Puntaje Z normal (1.95996 para 95% de confianza)
        
    Returns
    -------
    tuple (lo, hi)
        Límites inferior y superior en porcentaje [0, 100]
    """
    if n_tot <= 0:
        return 0.0, 0.0
    p = float(np.clip(k_aff / float(n_tot), 0.0, 1.0))
    denom = 1.0 + (z ** 2) / n_tot
    center = (p + (z ** 2) / (2.0 * n_tot)) / denom
    half = (z * math.sqrt((p * (1.0 - p) + (z ** 2) / (4.0 * n_tot)) / n_tot)) / denom
    lo = max(0.0, (center - half) * 100.0)
    hi = min(100.0, (center + half) * 100.0)
    return lo, hi


def cohen_h_proportions(p_exp_pct, p_ref_pct):
    """
    Calcula el tamaño del efecto h de Cohen para la diferencia entre dos proporciones:
      h = 2 * arcsin(sqrt(p_exp)) - 2 * arcsin(sqrt(p_ref))
      
    Interpretación bioestadística:
      |h| < 0.20: Sin afectación / Ruido basal
      0.20 <= h < 0.50: Afectación Leve
      0.50 <= h < 0.80: Afectación Moderada
      h >= 0.80: Afectación Severa / Lisis Masiva
    """
    p1 = float(np.clip(p_exp_pct / 100.0, 0.0, 1.0))
    p0 = float(np.clip(p_ref_pct / 100.0, 0.0, 1.0))
    return float(2.0 * math.asin(math.sqrt(p1)) - 2.0 * math.asin(math.sqrt(p0)))


def contingency_2x2_p_and_or(aff_exp, unaff_exp, aff_ref, unaff_ref):
    """
    Ejecuta prueba estadística sobre tabla de contingencia 2x2 [Afectados V2L vs No Afectados V2R]:
      - Chi-cuadrado de Pearson con corrección de Yates (o Exacta de Fisher cuando N < 200 o celdas < 5).
      - Odds Ratio (OR) con corrección de Haldane-Anscombe (+0.5) e IC 95% de Woolf.
      - Riesgo Relativo (RR = P_exp / P_ref).
      
    Returns
    -------
    tuple (chi2_stat, p_val, or_val, or_lo, or_hi, rr_val)
    """
    a = max(0, int(round(aff_exp)))
    b = max(0, int(round(unaff_exp)))
    c = max(0, int(round(aff_ref)))
    d = max(0, int(round(unaff_ref)))

    # Odds Ratio y su IC 95% con corrección Haldane-Anscombe
    a_c, b_c, c_c, d_c = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_val = float((a_c * d_c) / (b_c * c_c))
    se_ln_or = math.sqrt(1.0 / a_c + 1.0 / b_c + 1.0 / c_c + 1.0 / d_c)
    or_lo = math.exp(math.log(or_val) - 1.96 * se_ln_or)
    or_hi = math.exp(math.log(or_val) + 1.96 * se_ln_or)

    p_exp = a_c / (a_c + b_c)
    p_ref = c_c / (c_c + d_c)
    rr_val = float(p_exp / max(1e-6, p_ref))

    n_total = a + b + c + d
    if n_total == 0 or (a + b) == 0 or (c + d) == 0:
        return 0.0, 1.0, or_val, or_lo, or_hi, rr_val

    table = np.array([[a, b], [c, d]], dtype=np.int64)
    try:
        if n_total < 200 or np.min(table) < 5:
            _, p_val = stats.fisher_exact(table, alternative='two-sided')
            z_eq = stats.norm.isf(max(1e-16, min(0.5, p_val / 2.0)))
            chi2_stat = float(z_eq ** 2)
        else:
            chi2_stat, p_val, _, _ = stats.chi2_contingency(table, correction=True)
            chi2_stat = float(chi2_stat)
            p_val = float(p_val)
    except Exception:
        chi2_stat, p_val = 0.0, 1.0

    return chi2_stat, float(np.clip(p_val, 0.0, 1.0)), or_val, or_lo, or_hi, rr_val


def cochran_armitage_trend_test(doses, aff_counts, total_counts):
    """
    Prueba de tendencia dosis-respuesta de Cochran-Armitage para proporciones binomiales
    (Afectados V2L vs No Afectados V2R a través de concentraciones crecientes).
    
    Returns
    -------
    tuple (z_stat, p_val, slope)
    """
    w = np.asarray(doses, dtype=float)
    r = np.asarray(aff_counts, dtype=float)
    n = np.asarray(total_counts, dtype=float)
    valid = n > 0
    if np.sum(valid) < 2:
        return 0.0, 1.0, 0.0
    w, r, n = w[valid], r[valid], n[valid]
    N = float(np.sum(n))
    R = float(np.sum(r))
    if N <= 0 or R <= 0 or R >= N:
        return 0.0, 1.0, 0.0

    p_bar = R / N
    w_bar = float(np.sum(n * w) / N)
    num = float(np.sum(r * (w - w_bar)))
    denom_sq = p_bar * (1.0 - p_bar) * float(np.sum(n * ((w - w_bar) ** 2)))
    if denom_sq <= 1e-12:
        return 0.0, 1.0, 0.0
    z_stat = num / math.sqrt(denom_sq)
    p_val = float(2.0 * stats.norm.sf(abs(z_stat)))
    props_pct = (r / n) * 100.0
    slope = float(stats.linregress(w, props_pct).slope) if len(w) >= 2 and np.std(w) > 0 else 0.0
    return float(z_stat), float(np.clip(p_val, 0.0, 1.0)), slope


def interpolate_crossing_time(t_arr, y_arr, target_val, start_time=0.0):
    """
    Calcula por interpolación lineal entre ventanas temporales el instante exacto (s) en que
    la curva cruza por primera vez el umbral target_val (para t >= start_time).
    """
    mask = t_arr >= start_time
    t_sub = np.asarray(t_arr[mask], dtype=float)
    y_sub = np.asarray(y_arr[mask], dtype=float)
    if len(t_sub) == 0:
        return 600.0
    hits = np.where(y_sub >= target_val)[0]
    if len(hits) == 0:
        return 600.0
    idx = int(hits[0])
    if idx == 0:
        return float(t_sub[0])
    t0, t1 = float(t_sub[idx - 1]), float(t_sub[idx])
    y0, y1 = float(y_sub[idx - 1]), float(y_sub[idx])
    if y1 <= y0 + 1e-9:
        return t1
    frac = np.clip((target_val - y0) / (y1 - y0), 0.0, 1.0)
    return float(t0 + frac * (t1 - t0))


def classify_affectation_status(p_adj, cohen_h, max_v2l_pct):
    """Define estadísticamente y biofísicamente el grado de afectación celular."""
    if p_adj >= 0.05 or cohen_h < 0.20 or max_v2l_pct < 5.0:
        return "No Afectado / Basal (NS)"
    if cohen_h < 0.50 or max_v2l_pct < 20.0:
        return "Afectación Leve Significativa (*)"
    if cohen_h < 0.80 or max_v2l_pct < 50.0:
        return "Afectación Moderada (**)"
    return "Afectación Severa / Lisis Masiva (***)"


def adjust_pvalues(p_vals, method="Benjamini-Hochberg (FDR)"):
    """
    Ajusta p-valores por comparaciones múltiples:
    - Benjamini-Hochberg (FDR)
    - Holm-Bonferroni
    - Sin ajuste
    """
    p_arr = np.asarray(p_vals, dtype=float)
    n = len(p_arr)
    if n == 0:
        return np.array([], dtype=float)
    if "Sin ajuste" in method:
        return np.clip(p_arr, 0.0, 1.0)
    if "Holm" in method:
        order = np.argsort(p_arr)
        q_sorted = np.maximum.accumulate(p_arr[order] * (n - np.arange(n)))
        q_final = np.empty(n, float)
        q_final[order] = np.clip(q_sorted, 0.0, 1.0)
        return q_final
    # Benjamini-Hochberg por defecto
    order = np.argsort(p_arr)
    ranks = np.empty(n, int)
    ranks[order] = np.arange(1, n + 1)
    q_vals = p_arr * n / ranks
    q_vals_sorted = np.minimum.accumulate(q_vals[order][::-1])[::-1]
    q_final = np.empty(n, float)
    q_final[order] = np.clip(q_vals_sorted, 0.0, 1.0)
    return q_final
