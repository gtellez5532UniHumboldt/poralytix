"""
Motor de Extracción Cinética, Gating Morfológico, Doble Cut-Off y Ventanas Temporales.
Calcula series temporales de permeabilización, V0 temprana, V0 completa y V_active.
"""
import os
import numpy as np
import pandas as pd
import scipy.stats as stats
try:
    from scipy.integrate import trapezoid
except ImportError:
    from scipy.integrate import trapz as trapezoid

from ..io.fcs_parser import parse_fcs_native, resolve_perm_channel
from .statistics import (
    wilson_binomial_ci,
    cohen_h_proportions,
    contingency_2x2_p_and_or,
    interpolate_crossing_time,
    classify_affectation_status
)


def extract_tube_kinetics(
    fpath, is_cc, bin_size=15,
    fsc_min=50000.0, fsc_max=np.inf,
    ssc_min=10000.0, ssc_max=np.inf,
    perm_channel="FITC (Por defecto)",
    log_excl_v1l=3.00,
    log_cutoff=4.85
):
    """
    Filtra eventos viables por FSC, SSC y el Cut-Off 1 de Exclusión (V1L), calcula sobre las células
    retenidas (V1R) las métricas continuas y clasifica eventos en Afectados (V2L) vs No Afectados (V2R).
    
    Returns
    -------
    tuple (pd.DataFrame, dict)
        (DataFrame binned por ventanas de tiempo, caché de arrays numéricos para re-análisis interactivo)
    """
    meta, data = parse_fcs_native(fpath)
    cols = {str(c).upper().strip(): c for c in data.columns}

    fitc_col = resolve_perm_channel(data.columns, perm_channel)
    fsc_col = next((cols[c] for c in ['FSC-A', 'FSC-H', 'FSC'] if c in cols), None)
    ssc_col = next((cols[c] for c in ['SSC-A', 'SSC-H', 'SSC'] if c in cols), None)
    time_col = next((cols[c] for c in ['TIME', 'TIME_MS', 'T'] if c in cols), None)

    if not fitc_col or not fsc_col or not ssc_col or not time_col:
        raise KeyError(
            f"Canales requeridos ausentes en {os.path.basename(fpath)}. "
            f"Buscado: '{perm_channel}'. Disponibles: {list(data.columns)}"
        )

    total_raw = len(data)
    timestep = float(meta.get('$TIMESTEP', 0.00002))
    t0 = data[time_col].min()
    time_s_all = (data[time_col] - t0) * timestep
    if time_s_all.max() < 1.0 and (data[time_col].max() - t0) > 10.0:
        time_s_all = data[time_col] - t0

    data = data.copy()
    if is_cc:
        t_valid = time_s_all <= 60.0
        data = data.loc[t_valid].copy()
        data['Time_s'] = time_s_all.loc[t_valid] - 60.0
        data = data[data['Time_s'] >= -60.0].copy()
    else:
        t_valid = time_s_all <= 600.0
        data = data.loc[t_valid].copy()
        data['Time_s'] = time_s_all.loc[t_valid]

    data['Log10_FITC'] = np.log10(np.maximum(1.0, data[fitc_col].values))

    event_cache = {
        'time_s': data['Time_s'].values.astype(np.float32),
        'fsc': data[fsc_col].values.astype(np.float32),
        'ssc': data[ssc_col].values.astype(np.float32),
        'fitc_raw': data[fitc_col].values.astype(np.float32),
        'log10_fitc': data['Log10_FITC'].values.astype(np.float32),
        'total_raw': int(total_raw),
        'perm_col': str(fitc_col),
        'is_cc': bool(is_cc)
    }

    res = apply_gate_and_cutoffs_to_cache(
        fpath, event_cache, bin_size=bin_size,
        fsc_min=fsc_min, fsc_max=fsc_max,
        ssc_min=ssc_min, ssc_max=ssc_max,
        log_excl_v1l=log_excl_v1l, log_cutoff=log_cutoff,
        scope_label="Global"
    )
    return res, event_cache


def apply_gate_and_cutoffs_to_cache(
    fpath, ev_cache, bin_size=15,
    fsc_min=50000.0, fsc_max=np.inf,
    ssc_min=10000.0, ssc_max=np.inf,
    log_excl_v1l=3.00, log_cutoff=4.85,
    scope_label="Global"
):
    """
    Aplica gating morfológico y doble corte de fluorescencia sobre la caché de eventos en memoria.
    """
    fsc_arr = ev_cache['fsc']
    ssc_arr = ev_cache['ssc']
    log_arr = ev_cache['log10_fitc']
    fitc_arr = ev_cache['fitc_raw']
    time_s = ev_cache['time_s']
    total_raw = ev_cache['total_raw']
    fitc_col = ev_cache['perm_col']

    morph_mask = (fsc_arr > fsc_min) & (ssc_arr > ssc_min)
    if np.isfinite(fsc_max):
        morph_mask = morph_mask & (fsc_arr <= fsc_max)
    if np.isfinite(ssc_max):
        morph_mask = morph_mask & (ssc_arr <= ssc_max)

    morph_count = int(np.sum(morph_mask))
    if morph_count == 0:
        raise ValueError(f"Sin eventos tras filtro FSC/SSC en {os.path.basename(fpath)}")

    if log_excl_v1l > 0.0:
        v1r_mask = morph_mask & (log_arr >= log_excl_v1l)
    else:
        v1r_mask = morph_mask

    v1r_count = int(np.sum(v1r_mask))
    excluded_v1l_count = max(0, morph_count - v1r_count)
    if v1r_count == 0:
        raise ValueError(
            f"Todos los eventos fueron excluidos por el Cut-Off V1L (< {log_excl_v1l:.2f} log10) en {os.path.basename(fpath)}"
        )

    tw_arr = (time_s[v1r_mask] // bin_size) * bin_size
    log_v1r = log_arr[v1r_mask]
    fitc_v1r = fitc_arr[v1r_mask]
    fsc_v1r = fsc_arr[v1r_mask]
    ssc_v1r = ssc_arr[v1r_mask]

    is_aff_bool = (log_v1r < log_cutoff)
    aff_v2l = is_aff_bool.astype(np.float64) * 100.0
    aff_cnt = is_aff_bool.astype(np.int64)

    if len(log_v1r) >= 10:
        sk_val = float(stats.skew(log_v1r))
        ku_val = float(stats.kurtosis(log_v1r, fisher=False))
        sarle_bc = float((sk_val ** 2 + 1.0) / max(1e-6, ku_val))
    else:
        sarle_bc = 0.33

    df_ev = pd.DataFrame({
        'Time_Window': tw_arr,
        'FITC': fitc_v1r,
        'Log10_FITC': log_v1r,
        'FSC': fsc_v1r,
        'SSC': ssc_v1r,
        'Is_Affected_V2L': aff_v2l,
        'Affected_Count': aff_cnt
    })

    res = df_ev.groupby('Time_Window').agg(
        Median_FITC=('FITC', 'median'),
        Mean_FSC=('FSC', 'mean'),
        Mean_SSC=('SSC', 'mean'),
        Events=('FITC', 'count'),
        Affected_V2L_Events=('Affected_Count', 'sum'),
        Affected_V2L_Pct=('Is_Affected_V2L', 'mean')
    ).reset_index()

    res['Affected_V2L_Events'] = res['Affected_V2L_Events'].astype(int)
    res['Unaffected_V2R_Events'] = (res['Events'] - res['Affected_V2L_Events']).astype(int)
    res['Unaffected_V2R_Pct'] = 100.0 - res['Affected_V2L_Pct']
    res['Odds_V2L_V2R'] = (res['Affected_V2L_Events'] + 0.5) / (res['Unaffected_V2R_Events'] + 0.5)
    res['Sarle_Bimodality_BC'] = float(sarle_bc)
    res['Filepath'] = str(fpath)
    res['Perm_Channel_Used'] = str(fitc_col)
    res['Total_Raw_Events'] = int(total_raw)
    res['Morph_Events'] = int(morph_count)
    res['Excluded_V1L_Events'] = int(excluded_v1l_count)
    res['Retained_Events'] = int(v1r_count)
    res['Applied_FSC_Min'] = float(fsc_min)
    res['Applied_FSC_Max'] = float(fsc_max)
    res['Applied_SSC_Min'] = float(ssc_min)
    res['Applied_SSC_Max'] = float(ssc_max)
    res['Applied_V1L_Excl_Log10'] = float(log_excl_v1l)
    res['Applied_V2L_Cutoff_Log10'] = float(log_cutoff)
    res['Applied_Scope'] = str(scope_label)
    return res


def extract_slope_based_parameters(
    t_arr_full, p_arr_full, is_control_cell=False, start_time=0.0,
    manual_tlag=None, manual_tmax=None, v2l_arr_full=None
):
    """
    Calcula T_lag (inicio de fase activa) y T_max (meseta de lisis) y la tasa activa V_active (%/s).
    """
    t_arr_full = np.asarray(t_arr_full, dtype=float)
    p_arr_full = np.asarray(p_arr_full, dtype=float)

    post_mask = t_arr_full >= 0.0
    if not np.any(post_mask):
        return 60.0, 300.0, 0.0, 0.0

    t_post = t_arr_full[post_mask]
    p_post = p_arr_full[post_mask]
    max_p = float(np.max(p_post)) if len(p_post) > 0 else 0.0

    y_sig = p_post.copy()
    if v2l_arr_full is not None and len(v2l_arr_full) == len(t_arr_full):
        v2l_post = np.asarray(v2l_arr_full, dtype=float)[post_mask]
        v2l_base = float(np.median(v2l_post[:min(3, len(v2l_post))])) if len(v2l_post) > 0 else 0.0
        v2l_net = np.clip(v2l_post - v2l_base, 0.0, 100.0)
        p_base_init = float(np.median(p_post[:min(3, len(p_post))])) if len(p_post) > 0 else 0.0
        p_span = max_p - p_base_init
        v2l_span = float(np.max(v2l_net)) if len(v2l_net) > 0 else 0.0
        if p_span < 5.0 and v2l_span >= p_span:
            y_sig = v2l_net
        elif v2l_span > 1.25 * p_span:
            y_sig = 0.5 * p_post + 0.5 * v2l_net

    if len(y_sig) >= 3:
        y_smooth = y_sig.copy()
        y_smooth[1:-1] = 0.25 * y_sig[:-2] + 0.50 * y_sig[1:-1] + 0.25 * y_sig[2:]
    else:
        y_smooth = y_sig.copy()

    dt_bin = float(np.median(np.diff(t_post))) if len(t_post) >= 2 else 15.0
    dt_bin = max(5.0, dt_bin)
    t_end_max = float(t_post[-1]) if len(t_post) > 0 else 600.0
    eff_start = float(np.clip(start_time, 0.0, max(0.0, t_end_max - 4.0 * dt_bin)))

    if manual_tlag is not None and manual_tmax is not None:
        t_lag = float(np.clip(manual_tlag, 0.0, max(0.0, t_end_max - dt_bin)))
        t_max = float(np.clip(manual_tmax, t_lag + dt_bin, t_end_max))
    else:
        valid_mask = t_post >= eff_start
        t_val = t_post[valid_mask] if np.any(valid_mask) else t_post
        y_val = y_smooth[valid_mask] if np.any(valid_mask) else y_smooth

        n_pts = len(t_val)
        if n_pts < 3:
            t_lag, t_max = eff_start, min(t_end_max, eff_start + 180.0)
        else:
            y_base = float(np.percentile(y_val[:min(4, n_pts)], 25.0))
            peak_idx_rel = int(np.argmax(y_val))
            y_peak = float(y_val[peak_idx_rel])
            net_amp = y_peak - y_base

            win_w = min(4, max(2, n_pts // 5))
            best_slope = -1e9
            best_idx = 0
            search_limit = max(win_w, peak_idx_rel + 1)
            for i in range(0, min(n_pts - win_w + 1, search_limit)):
                tw = t_val[i:i + win_w]
                yw = y_val[i:i + win_w]
                if np.std(tw) > 1e-6:
                    sl, _, _, _, _ = stats.linregress(tw, yw)
                    if sl > best_slope:
                        best_slope = float(sl)
                        best_idx = i

            if net_amp < 1.8 or best_slope <= 0.002:
                if best_slope > 0.0 and n_pts >= 4:
                    t_lag_raw = float(t_val[best_idx])
                    t_max_raw = min(t_end_max, t_lag_raw + max(120.0, 6.0 * dt_bin))
                else:
                    t_lag_raw = max(eff_start, 60.0)
                    t_max_raw = min(t_end_max, max(t_lag_raw + 120.0, 300.0))
                t_lag = float(np.clip(round(t_lag_raw / dt_bin) * dt_bin, eff_start, max(eff_start, t_end_max - 2.0 * dt_bin)))
                t_max = float(np.clip(round(t_max_raw / dt_bin) * dt_bin, t_lag + 2.0 * dt_bin, t_end_max))
            else:
                mid_i = min(n_pts - 1, best_idx + win_w // 2)
                t_infl = float(t_val[mid_i])
                y_infl = float(y_val[mid_i])
                t_lag_tan = t_infl - (y_infl - y_base) / max(1e-5, best_slope)
                t_max_tan = t_infl + (y_peak - y_infl) / max(1e-5, best_slope)

                thresh_12 = y_base + 0.12 * net_amp
                thresh_90 = y_base + 0.90 * net_amp

                idx_12 = np.where(y_val >= thresh_12)[0]
                t_12 = float(t_val[idx_12[0]]) if len(idx_12) > 0 else eff_start

                if eff_start <= t_lag_tan <= t_infl:
                    t_lag_cand = 0.60 * t_lag_tan + 0.40 * t_12
                else:
                    t_lag_cand = t_12
                t_lag_cand = np.clip(t_lag_cand, eff_start, max(eff_start, float(t_val[peak_idx_rel]) - 2.0 * dt_bin))
                t_lag = float(t_val[int(np.argmin(np.abs(t_val - t_lag_cand)))])

                idx_90 = np.where((t_val > t_lag) & (y_val >= thresh_90))[0]
                t_90 = float(t_val[idx_90[0]]) if len(idx_90) > 0 else float(t_val[peak_idx_rel])
                if t_max_tan > t_lag + dt_bin:
                    t_max_cand = 0.50 * min(t_end_max, t_max_tan) + 0.50 * t_90
                else:
                    t_max_cand = t_90
                t_max_cand = np.clip(t_max_cand, t_lag + 2.0 * dt_bin, t_end_max)
                t_max = float(t_val[int(np.argmin(np.abs(t_val - t_max_cand)))])
                if t_max <= t_lag:
                    t_max = min(t_end_max, t_lag + max(30.0, 2.0 * dt_bin))

    idx_phase = np.where((t_post >= t_lag) & (t_post <= t_max))[0]
    if len(idx_phase) >= 3 and np.std(t_post[idx_phase]) > 1e-6:
        slope_p, _, _, _, _ = stats.linregress(t_post[idx_phase], p_post[idx_phase])
        slope_y, _, _, _, _ = stats.linregress(t_post[idx_phase], y_sig[idx_phase])
        v_active = max(0.0, float(slope_p if slope_p >= 0.005 else max(slope_p, 0.35 * slope_y)))
    elif len(idx_phase) >= 2:
        dt = t_post[idx_phase[-1]] - t_post[idx_phase[0]]
        dp = p_post[idx_phase[-1]] - p_post[idx_phase[0]]
        dy = y_sig[idx_phase[-1]] - y_sig[idx_phase[0]]
        s_p = dp / dt if dt > 0 else 0.0
        s_y = dy / dt if dt > 0 else 0.0
        v_active = max(0.0, float(s_p if s_p >= 0.005 else max(s_p, 0.35 * s_y)))
    else:
        v_active = 0.0

    if is_control_cell and max_p < 4.5 and (np.max(y_sig) - np.min(y_sig)) < 4.5:
        v_active = min(v_active, 0.005)

    return t_lag, t_max, v_active, max_p


def compute_categorical_kinetics_metrics(
    t_arr, v2l_pct_arr, cont_perm_arr, start_time=0.0,
    events_arr=None, v2l_events_arr=None, sarle_bc=0.33
):
    """
    Calcula parámetros bioestadísticos integrales basados en las proporciones categóricas
    de eventos Afectados (V2L %) vs No Afectados (V2R %).
    """
    base_mask = (t_arr >= -60.0) & (t_arr < 0.0)
    post_mask = t_arr >= 0.0
    t_p = np.asarray(t_arr[post_mask], dtype=float)
    v2l_p = np.clip(np.asarray(v2l_pct_arr[post_mask], dtype=float), 0.0, 100.0)
    cont_p = np.clip(np.asarray(cont_perm_arr[post_mask], dtype=float), 0.0, 100.0)

    if len(t_p) < 2:
        return {
            'Basal_V2L_Pct': 0.0, 'Min_V2L_Pct': 0.0, 'Max_V2L_Pct': 0.0, 'End_V2L_Pct': 0.0,
            'Min_Unaffected_V2R_Pct': 100.0, 'Max_Unaffected_V2R_Pct': 100.0, 'End_Unaffected_V2R_Pct': 100.0,
            'T_Effect_Min_s': 0.0, 'T_Effect_Max_s': 600.0,
            'T10_Effect_s': 600.0, 'T50_Effect_s': 600.0, 'T90_Effect_s': 600.0, 'Span_T10_T90_s': 0.0,
            'LT50_s': 600.0, 'T_Peak_Conversion_s': 600.0, 'V_Cat_Max_Inst': 0.0,
            'T_Stat_Onset_s': 600.0, 'AUC_V2L_Pct': 0.0, 'V_Cat_Active': 0.0,
            'Max_Entropy_Bits': 0.0, 'Hazard_Rate_s': 0.0, 'Sync_Index': 0.0,
            'Bimodality_BC': float(sarle_bc), 'Cohen_h_Max': 0.0, 'OR_vs_Basal': 1.0,
            'RR_vs_Basal': 1.0, 'Chi2_Max': 0.0, 'p_Fisher_Max': 1.0,
            'Affectation_State': 'No Afectado / Basal (NS)',
            'Mode_Delta': 0.0, 'Mechanism_Class': 'Inerte / Protegido'
        }

    basal_v2l = float(np.mean(v2l_pct_arr[base_mask])) if np.any(base_mask) else float(v2l_p[0])
    min_idx = int(np.argmin(v2l_p))
    max_idx = int(np.argmax(v2l_p))
    min_v2l = float(v2l_p[min_idx])
    max_v2l = float(v2l_p[max_idx])
    end_v2l = float(v2l_p[-1])

    min_v2r = float(100.0 - max_v2l)
    max_v2r = float(100.0 - min_v2l)
    end_v2r = float(100.0 - end_v2l)

    t_eff_min = float(t_p[min_idx])
    net_amp = max_v2l - min_v2l
    if net_amp >= 3.0:
        near_max_thresh = min_v2l + 0.98 * net_amp
        hits_max = t_p[v2l_p >= near_max_thresh]
        t_eff_max = float(hits_max[0]) if len(hits_max) > 0 else float(t_p[max_idx])
    else:
        t_eff_max = float(t_p[max_idx]) if max_v2l >= 5.0 else 600.0

    if net_amp >= 3.0 and max_v2l >= 5.0:
        y_10 = min_v2l + 0.10 * net_amp
        y_50 = min_v2l + 0.50 * net_amp
        y_90 = min_v2l + 0.90 * net_amp
        t10_eff = interpolate_crossing_time(t_p, v2l_p, y_10, start_time=start_time)
        t50_eff = interpolate_crossing_time(t_p, v2l_p, y_50, start_time=start_time)
        t90_eff = interpolate_crossing_time(t_p, v2l_p, y_90, start_time=start_time)
        span_10_90 = max(0.0, t90_eff - t10_eff)
    else:
        t10_eff, t50_eff, t90_eff, span_10_90 = 600.0, 600.0, 600.0, 0.0

    lt50 = interpolate_crossing_time(t_p, v2l_p, 50.0, start_time=0.0)

    dt_arr = np.gradient(t_p)
    dt_arr[dt_arr == 0] = 15.0
    inst_rate = np.gradient(v2l_p) / dt_arr
    valid_rate_mask = t_p >= start_time
    if np.any(valid_rate_mask):
        idx_sub = np.where(valid_rate_mask)[0]
        peak_rel = int(np.argmax(inst_rate[valid_rate_mask]))
        peak_rate_idx = int(idx_sub[peak_rel])
        v_cat_max_inst = max(0.0, float(inst_rate[peak_rate_idx]))
        t_peak_conv = float(t_p[peak_rate_idx]) if v_cat_max_inst >= 0.02 else 600.0
    else:
        v_cat_max_inst = 0.0
        t_peak_conv = 600.0

    auc_v2l = float(trapezoid(y=v2l_p, x=t_p) / (600.0 * 100.0) * 100.0)
    auc_cont = float(trapezoid(y=cont_p, x=t_p) / (600.0 * 100.0) * 100.0)

    _, _, v_cat, _ = extract_slope_based_parameters(
        t_arr, v2l_pct_arr, is_control_cell=(max_v2l < 5.0), start_time=start_time
    )

    p_prob = np.clip(v2l_p / 100.0, 1e-6, 1.0 - 1e-6)
    shannon_h = -(p_prob * np.log2(p_prob) + (1.0 - p_prob) * np.log2(1.0 - p_prob))
    shannon_h[(v2l_p < 1.0) | (v2l_p > 99.0)] = 0.0
    max_entropy = float(np.max(shannon_h))

    surv_prob = np.clip(1.0 - (v2l_p / 100.0), 0.02, 1.0)
    cum_hazard = -np.log(surv_prob)
    if len(t_p) >= 2 and np.std(t_p) > 1e-6:
        hazard_slope, _, _, _, _ = stats.linregress(t_p, cum_hazard)
        hazard_rate = max(0.0, float(hazard_slope))
    else:
        hazard_rate = 0.0

    sync_index = float(t50_eff / max(15.0, span_10_90)) if (t50_eff < 600.0 and span_10_90 > 0) else 0.0

    if events_arr is not None and v2l_events_arr is not None and len(events_arr) == len(t_arr):
        ev_all = np.asarray(events_arr, dtype=int)
        aff_all = np.asarray(v2l_events_arr, dtype=int)
        if np.any(base_mask):
            b_tot = int(np.sum(ev_all[base_mask]))
            b_aff = int(np.sum(aff_all[base_mask]))
        else:
            b_tot = int(ev_all[post_mask][0]) if len(ev_all[post_mask]) > 0 else 100
            b_aff = int(aff_all[post_mask][0]) if len(aff_all[post_mask]) > 0 else int(round(basal_v2l))
        b_unaff = max(0, b_tot - b_aff)

        ev_p = ev_all[post_mask]
        aff_p = aff_all[post_mask]
        m_tot = int(ev_p[max_idx])
        m_aff = int(aff_p[max_idx])
        m_unaff = max(0, m_tot - m_aff)

        chi2_max, p_fisher_max, or_basal, _, _, rr_basal = contingency_2x2_p_and_or(m_aff, m_unaff, b_aff, b_unaff)
        cohen_h_max = cohen_h_proportions(max_v2l, basal_v2l)

        t_stat_onset = 600.0
        for w_i in range(len(t_p)):
            if t_p[w_i] < start_time:
                continue
            w_aff = int(aff_p[w_i])
            w_unaff = max(0, int(ev_p[w_i]) - w_aff)
            _, pv_w, _, _, _, _ = contingency_2x2_p_and_or(w_aff, w_unaff, b_aff, b_unaff)
            h_w = cohen_h_proportions(v2l_p[w_i], basal_v2l)
            if pv_w < 0.001 and h_w >= 0.20:
                t_stat_onset = float(t_p[w_i])
                break
    else:
        cohen_h_max = cohen_h_proportions(max_v2l, basal_v2l)
        chi2_max, p_fisher_max, or_basal, _, _, rr_basal = contingency_2x2_p_and_or(
            int(round(max_v2l * 2.0)), int(round((100.0 - max_v2l) * 2.0)),
            int(round(basal_v2l * 2.0)), int(round((100.0 - basal_v2l) * 2.0))
        )
        onset_hits = t_p[(t_p >= start_time) & (v2l_p >= max(5.0, basal_v2l + 5.0))]
        t_stat_onset = float(onset_hits[0]) if len(onset_hits) > 0 else 600.0

    aff_state = classify_affectation_status(p_fisher_max, cohen_h_max, max_v2l)

    mode_delta = auc_cont - auc_v2l
    if max_v2l < 5.0 and np.max(cont_p) < 5.0:
        mech_class = "Inerte / Protegido"
    elif abs(mode_delta) <= 8.0 or sarle_bc >= 0.555:
        mech_class = "Todo-o-Nada (Poro Bimodal)"
    elif mode_delta > 8.0:
        mech_class = "Gradual (Fuga Parcial)"
    else:
        mech_class = "Umbral Temprano V2L"

    return {
        'Basal_V2L_Pct': basal_v2l,
        'Min_V2L_Pct': min_v2l,
        'Max_V2L_Pct': max_v2l,
        'End_V2L_Pct': end_v2l,
        'Min_Unaffected_V2R_Pct': min_v2r,
        'Max_Unaffected_V2R_Pct': max_v2r,
        'End_Unaffected_V2R_Pct': end_v2r,
        'T_Effect_Min_s': t_eff_min,
        'T_Effect_Max_s': t_eff_max,
        'T10_Effect_s': t10_eff,
        'T50_Effect_s': t50_eff,
        'T90_Effect_s': t90_eff,
        'Span_T10_T90_s': span_10_90,
        'LT50_s': lt50,
        'T_Peak_Conversion_s': t_peak_conv,
        'V_Cat_Max_Inst': v_cat_max_inst,
        'T_Stat_Onset_s': t_stat_onset,
        'AUC_V2L_Pct': auc_v2l,
        'V_Cat_Active': v_cat,
        'Max_Entropy_Bits': max_entropy,
        'Hazard_Rate_s': hazard_rate,
        'Sync_Index': sync_index,
        'Bimodality_BC': float(sarle_bc),
        'Cohen_h_Max': cohen_h_max,
        'OR_vs_Basal': or_basal,
        'RR_vs_Basal': rr_basal,
        'Chi2_Max': chi2_max,
        'p_Fisher_Max': p_fisher_max,
        'Affectation_State': aff_state,
        'Mode_Delta': mode_delta,
        'Mechanism_Class': mech_class
    }
