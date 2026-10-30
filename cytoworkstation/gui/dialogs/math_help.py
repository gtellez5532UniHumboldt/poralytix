"""
Modal de Referencia Matemática, Biofísica y Bioestadística de la Estación de Trabajo.
Soporte bilingüe Español e Inglés.
"""
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from pathlib import Path

from ...i18n import get_current_language


def open_mathematical_models_help_dialog(parent, initial_tab=0, lang=None):
    """Abre la ventana interactiva explicando modelos matemáticos, ecuaciones y estadística."""
    active_lang = lang or get_current_language()
    dlg = tk.Toplevel(parent)
    
    if active_lang == 'en':
        dlg.title("📐 Mathematical Models, Biophysical Equations & Statistical Calculations Reference")
        hdr_txt = "📐 MATHEMATICAL AND BIOSTATISTICAL REFERENCE OF THE CYTOMETRY WORKSTATION"
        tabs_info = [
            (" 📐 1. Pharmacodynamic & Kinetic Models ", [
                ("1.1 Regularized Allosteric Hill Model (Default)",
                 "Equation:\n    V_active(C) = ( V_max · C^n_H ) / ( K_d^n_H + C^n_H )\n\n"
                 "Parameters and Biophysical Meaning:\n"
                 "  • V_max (%/s): Maximum lysis rate at saturation.\n"
                 "  • K_d (µM): Semi-saturation concentration (kinetic EC50).\n"
                 "  • n_H: Hill coefficient in [0.5, 6.0]. n_H > 1 indicates positive cooperativity / oligomerization.\n"
                 "  • Lytic Efficiency: V_max / K_d (% · s⁻¹ · µM⁻¹)."),
                ("1.2 Hyperbolic Michaelis-Menten Model",
                 "Equation:\n    V_active(C) = ( V_max · C ) / ( K_d + C )\n\n"
                 "  • First-order hyperbolic saturation kinetics (fixed n_H = 1.0)."),
                ("1.3 Asymmetric Gompertz Sigmoidal Model",
                 "Equation:\n    V_active(C) = V_max · exp( -exp( -k_G · (C - K_d) ) )\n\n"
                 "  • Captures sharp threshold transitions once critical concentration K_d is exceeded."),
                ("1.4 Weibull / Stretched Exponential Model",
                 "Equation:\n    V_active(C) = V_max · [ 1 - exp( -( C / K_d )^m ) ]\n\n"
                 "  • Based on percolation theory and membrane resistance distributions."),
                ("1.5 Sequential Bi-Exponential Time Model",
                 "Equation:\n    P(t) = P_ads · (1 - exp(-k_ads · t)) + P_lys · (1 - exp(-k_lys · max(0, t - T_lag)))\n\n"
                 "  • Decouples early electrostatic adsorption from secondary pore-forming lytic phases.")
            ]),
            (" 🧬 2. Categorical Panel & Temporal Milestones ", [
                ("2.1 Temporal Milestones of Categorical Effect",
                 "  • Statistical T_onset (s): First time window showing effect with p < 0.001 and h >= 0.20.\n"
                 "  • Effect T50 (s): Interpolated time to reach 50% of the maximum effect in the curve.\n"
                 "  • Population LT50 (s): Time where absolute percentage of affected events exceeds 50%.\n"
                 "  • T_effect_max (s): Time window where maximal plateau effect is achieved (>= 98% amplitude)."),
                ("2.2 Mechanistic Classification (All-or-None vs Gradual)",
                 "  • Delta Mode = Continuous_AUC - Categorical_AUC.\n"
                 "  • All-or-None: Sudden bimodal pore depolarization (|Δ| <= 8% or Sarle BC >= 0.555).\n"
                 "  • Gradual: Progressive fluorophore leakage without clear bimodality (Δ > 8%).")
            ]),
            (" 📊 3. Biostatistics & Hypothesis Testing ", [
                ("3.1 Exact 95% Wilson Binomial Confidence Interval",
                 "Wilson Formula:\n    Center = (p + z² / (2N)) / (1 + z² / N)\n"
                 "    Half = (z / (1 + z²/N)) · sqrt( p(1-p)/N + z²/(4N²) )\n"
                 "  • Avoids degraded zero-variance bounds near 0% or 100% suffered by the Wald approximation."),
                ("3.2 Cohen's h Effect Size for Proportions",
                 "    h = 2 · arcsin(sqrt(P_exp)) - 2 · arcsin(sqrt(P_ref))\n"
                 "  • |h| < 0.20: Basal noise / No effect\n"
                 "  • 0.20 <= h < 0.50: Small effect\n"
                 "  • 0.50 <= h < 0.80: Moderate effect\n"
                 "  • h >= 0.80: Large effect / Massive lysis"),
                ("3.3 2x2 Tables: Odds Ratio, Chi-Square and Fisher's Exact Test",
                 "  • Odds Ratio with Haldane-Anscombe (+0.5) correction and Woolf 95% CI.\n"
                 "  • Pearson's Chi-square with Yates' continuity correction or Fisher's Exact Test.")
            ])
        ]
    else:
        dlg.title("📐 Manual de Modelos Matemáticos, Ecuaciones Biofísicas y Cálculos Estadísticos")
        hdr_txt = "📐 REFERENCIA MATEMÁTICA Y BIOESTADÍSTICA DE LA ESTACIÓN DE TRABAJO DE CITOMETRÍA"
        tabs_info = [
            (" 📐 1. Modelos Farmacodinámicos y Cinéticos ", [
                ("1.1 Modelo Alostérico de Hill Regularizado (Por defecto)",
                 "Ecuación:\n    V_active(C) = ( V_max · C^n_H ) / ( K_d^n_H + C^n_H )\n\n"
                 "Parámetros y Significado Biofísico:\n"
                 "  • V_max (%/s): Velocidad máxima de lisis a saturación.\n"
                 "  • K_d (µM): Concentración de semi-saturación (EC50 cinética).\n"
                 "  • n_H: Coeficiente de Hill en [0.5, 6.0]. n_H > 1 indica cooperatividad positiva / oligomerización.\n"
                 "  • Eficiencia Lítica: V_max / K_d (% · s⁻¹ · µM⁻¹)."),
                ("1.2 Modelo Hiperbólico de Michaelis-Menten",
                 "Ecuación:\n    V_active(C) = ( V_max · C ) / ( K_d + C )\n\n"
                 "  • Cinética de saturación hiperbólica de primer orden (n_H = 1.0 fijo)."),
                ("1.3 Modelo Sigmoidal Asimétrico de Gompertz",
                 "Ecuación:\n    V_active(C) = V_max · exp( -exp( -k_G · (C - K_d) ) )\n\n"
                 "  • Captura transiciones abruptas una vez superada la concentración umbral crítica K_d."),
                ("1.4 Modelo de Weibull / Exponencial Estirado",
                 "Ecuación:\n    V_active(C) = V_max · [ 1 - exp( -( C / K_d )^m ) ]\n\n"
                 "  • Basado en teoría de percolación y distribución de resistencia de membrana."),
                ("1.5 Modelo Bi-Exponencial Secuencial en Tiempo",
                 "Ecuación:\n    P(t) = P_ads · (1 - exp(-k_ads · t)) + P_lys · (1 - exp(-k_lys · max(0, t - T_lag)))\n\n"
                 "  • Desacopla la fase rápida de adsorción electrostática de la fase lítica secundaria.")
            ]),
            (" 🧬 2. Panel Categórico y Hitos Temporales ", [
                ("2.1 Hitos Temporales del Efecto Categórico",
                 "  • T_onset estadístico (s): Primera ventana temporal con afectación p < 0.001 y h >= 0.20.\n"
                 "  • T50 del Efecto (s): Tiempo interpolado para alcanzar el 50% del efecto propio de la curva.\n"
                 "  • LT50 poblacional (s): Instante donde el % absoluto de eventos afectados supera el 50%.\n"
                 "  • T_efecto_máx (s): Ventana de tiempo donde ocurre el efecto máximo (>= 98% de la amplitud)."),
                ("2.2 Clasificación Mecanística (Todo-o-Nada vs Gradual)",
                 "  • Modo Delta = AUC_Continua - AUC_Categórica.\n"
                 "  • Todo-o-Nada: Despolarización bimodal súbita del poro (|Δ| <= 8% o Sarle BC >= 0.555).\n"
                 "  • Gradual: Fuga gradual de fluoróforo sin bimodalidad neta (Δ > 8%).")
            ]),
            (" 📊 3. Bioestadística y Pruebas de Hipótesis ", [
                ("3.1 Intervalo Binomial Exacto de Wilson al 95%",
                 "Fórmula de Wilson:\n    Center = (p + z² / (2N)) / (1 + z² / N)\n"
                 "    Half = (z / (1 + z²/N)) · sqrt( p(1-p)/N + z²/(4N²) )\n"
                 "  • Evita sesgos degenerados cerca de 0% o 100% que sufre la aproximación clásica de Wald."),
                ("3.2 Tamaño del Efecto h de Cohen para Proporciones",
                 "    h = 2 · arcsin(sqrt(P_exp)) - 2 · arcsin(sqrt(P_ref))\n"
                 "  • |h| < 0.20: Ruido basal / Sin efecto\n"
                 "  • 0.20 <= h < 0.50: Afectación Leve\n"
                 "  • 0.50 <= h < 0.80: Afectación Moderada\n"
                 "  • h >= 0.80: Afectación Severa / Lisis Masiva"),
                ("3.3 Tablas 2x2: Odds Ratio, Chi2 y Fisher",
                 "  • Odds Ratio con corrección Haldane-Anscombe (+0.5) e IC 95% de Woolf.\n"
                 "  • Chi-cuadrado de Pearson con corrección de continuidad de Yates o Test Exacto de Fisher.")
            ])
        ]

    dlg.geometry("980x680")
    dlg.transient(parent)

    hdr = tk.Frame(dlg, bg="#002B49", pady=6, padx=12)
    hdr.pack(fill="x")
    tk.Label(
        hdr,
        text=hdr_txt,
        bg="#002B49", fg="white", font=("DejaVu Sans", 10, "bold")
    ).pack(side="left")

    nb = ttk.Notebook(dlg)
    nb.pack(fill="both", expand=True, padx=10, pady=8)

    for title, sections in tabs_info:
        frm = ttk.Frame(nb, padding=10)
        nb.add(frm, text=title)

        canvas = tk.Canvas(frm, bg="#FDFEFE", highlightthickness=0)
        scroll = ttk.Scrollbar(frm, orient="vertical", command=canvas.yview)
        sub = tk.Frame(canvas, bg="#FDFEFE")
        sub.bind("<Configure>", lambda e, c=canvas: c.configure(scrollregion=c.bbox("all")))
        canvas.create_window((0, 0), window=sub, anchor="nw")
        canvas.configure(yscrollcommand=scroll.set)
        canvas.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")

        for sec_t, sec_body in sections:
            lf = ttk.LabelFrame(sub, text=f" {sec_t} ", padding=8)
            lf.pack(fill="x", padx=10, pady=6)
            ttk.Label(lf, text=sec_body, justify="left", font=("Consolas", 9)).pack(anchor="w")

    nb.select(initial_tab)
    return dlg
