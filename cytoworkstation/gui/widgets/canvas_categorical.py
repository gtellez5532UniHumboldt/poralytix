"""
Lienzo interactivo Tkinter para el Panel de Análisis Categórico (Afectados V2L vs No Afectados V2R).
Soporta:
  1) Balance Temporal Complementario (% V2L vs % V2R + hitos temporales)
  2) Comparativa de Hitos Temporales (T_onset, T50, T_max)
  3) Mapa de Calor Espaciotemporal
  4) Velocidad Instantánea de Conversión d(%V2L)/dt
"""
import tkinter as tk
import numpy as np
import pandas as pd
from ...config import CURVE_PALETTE


class NativeCategoricalProportionsCanvas(tk.Frame):
    """
    Lienzo interactivo en Tkinter para análisis categórico y ventanas temporales.
    """
    def __init__(self, parent):
        super().__init__(parent, bg="white")
        self.canvas = tk.Canvas(self, bg="white", highlightthickness=0, height=235)
        self.canvas.pack(fill="both", expand=True)
        self.mode = "Balance Temporal: % Afectados (V2L) vs % No Afectados (V2R) y Hitos"
        self.win_df = pd.DataFrame()
        self.sum_df = pd.DataFrame()
        self.active_label = ""
        self.canvas.bind("<Configure>", lambda e: self.redraw())

    def set_panel_data(self, mode, win_df, sum_df, active_label=""):
        self.mode = mode
        self.win_df = win_df if win_df is not None else pd.DataFrame()
        self.sum_df = sum_df if sum_df is not None else pd.DataFrame()
        self.active_label = active_label
        self.redraw()

    def redraw(self):
        self.canvas.delete("all")
        W = max(520, self.canvas.winfo_width())
        H = max(195, self.canvas.winfo_height())

        if self.win_df.empty and self.sum_df.empty:
            self.canvas.create_text(
                W // 2, H // 2,
                text="Ejecuta el Pre-Análisis o carga el Set Demo para visualizar el Panel de Proporciones Categóricas.",
                font=("DejaVu Sans", 9, "bold"), fill="#566573"
            )
            return

        if "Mapa de Calor" in self.mode:
            self._draw_heatmap(W, H)
        elif "Comparativa de Hitos" in self.mode:
            self._draw_milestones_bars(W, H)
        elif "Velocidad" in self.mode:
            self._draw_conversion_rate_and_entropy(W, H)
        else:
            self._draw_complementary_proportions_curve(W, H)

    def _draw_complementary_proportions_curve(self, W, H):
        ml, mr, mt, mb = 68, 250, 28, 34
        x0, y0, x1, y1 = ml, mt, max(ml + 180, W - mr), max(mt + 90, H - mb)
        pw, ph = x1 - x0, y1 - y0
        xmin, xmax = -60.0, 600.0

        def to_px(tv, yv):
            px = x0 + ((tv - xmin) / (xmax - xmin)) * pw
            py = y1 - (np.clip(yv, 0.0, 100.0) / 100.0) * ph
            return px, py

        self.canvas.create_rectangle(x0, y0, x1, y1, fill="#FCFDFE", outline="#333333")
        title_txt = f"Proporción de Eventos Afectados (V2L %) vs. No Afectados (V2R %) — {self.active_label}"
        self.canvas.create_text((x0 + x1) // 2, y0 - 12, text=title_txt, font=("DejaVu Sans", 9, "bold"), fill="#002B49")

        if self.win_df.empty:
            return

        sub = self.win_df.copy()
        sub['CID'] = sub.apply(lambda r: f"{r['Peptide']} | {r['Cell_Line']} | {r['Concentration_uM']:.2f} µM", axis=1)
        unique_cids = sorted(sub['CID'].unique())
        chosen_cid = self.active_label if self.active_label in unique_cids else (unique_cids[-1] if unique_cids else "")
        if not chosen_cid:
            return
        c_df = sub[sub['CID'] == chosen_cid].sort_values('Time_Window')

        t_vals = c_df['Time_Window'].values.astype(float)
        v2l_vals = c_df['Mean_Affected_V2L_Pct'].values.astype(float)
        v2r_vals = c_df['Mean_Unaffected_V2R_Pct'].values.astype(float)
        sig_flags = c_df['Window_Status'].values if 'Window_Status' in c_df else []

        if len(t_vals) >= 2:
            poly_v2l = [to_px(t_vals[0], 0.0)[0], y1]
            for tv, yv in zip(t_vals, v2l_vals):
                poly_v2l.extend(to_px(tv, yv))
            poly_v2l.extend([to_px(t_vals[-1], 0.0)[0], y1])
            self.canvas.create_polygon(poly_v2l, fill="#FADBD8", outline="")

            poly_v2r = [to_px(t_vals[0], 100.0)[0], y0]
            for tv, yv in zip(t_vals, v2l_vals):
                poly_v2r.extend(to_px(tv, yv))
            poly_v2r.extend([to_px(t_vals[-1], 100.0)[0], y0])
            self.canvas.create_polygon(poly_v2r, fill="#D5F5E3", outline="")

            pts_v2l, pts_v2r = [], []
            for idx_pt, (tv, v_aff, v_unaff) in enumerate(zip(t_vals, v2l_vals, v2r_vals)):
                px, py_a = to_px(tv, v_aff)
                _, py_u = to_px(tv, v_unaff)
                pts_v2l.extend([px, py_a])
                pts_v2r.extend([px, py_u])
                self.canvas.create_oval(px - 2.5, py_a - 2.5, px + 2.5, py_a + 2.5, fill="#C0392B", outline="white")
                self.canvas.create_rectangle(px - 2.0, py_u - 2.0, px + 2.0, py_u + 2.0, fill="#1E8449", outline="")
                if idx_pt < len(sig_flags) and ("*" in str(sig_flags[idx_pt])):
                    self.canvas.create_text(px, max(y0 + 6, py_a - 7), text="*", font=("DejaVu Sans", 8, "bold"), fill="#922B21")

            self.canvas.create_line(pts_v2l, fill="#C0392B", width=2.5)
            self.canvas.create_line(pts_v2r, fill="#1E8449", width=2.0, dash=(4, 2))

        for tx in [-60, 0, 60, 120, 180, 240, 300, 360, 420, 480, 540, 600]:
            px, _ = to_px(tx, 0)
            self.canvas.create_line(px, y0, px, y1, fill="#D5D8DC" if tx != 0 else "#333333", dash=(2, 3) if tx != 0 else (4, 2))
            self.canvas.create_text(px, y1 + 10, text=f"{tx}s", font=("DejaVu Sans", 8))

        for gy in [0, 25, 50, 75, 100]:
            _, py = to_px(0, gy)
            self.canvas.create_line(x0, py, x1, py, fill="#BFC9CA" if gy == 50 else "#E5E8E8", dash=(4, 2) if gy == 50 else (2, 3))
            self.canvas.create_text(x0 - 6, py, text=f"{gy}%", anchor="e", font=("DejaVu Sans", 8))

        self.canvas.create_text(x0 - 44, (y0 + y1) // 2, text="Proporción (%)", font=("DejaVu Sans", 8, "bold"), fill="#1C2833")
        self.canvas.create_text((x0 + x1) // 2, y1 + 23, text="Ventana de Tiempo (segundos)", font=("DejaVu Sans", 8, "bold"), fill="#002B49")

        lx, ly = x1 + 10, y0 + 2
        self.canvas.create_text(lx, ly, text=f"Curva: {chosen_cid[:28]}", anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#002B49")
        ly += 16
        self.canvas.create_line(lx, ly + 6, lx + 16, ly + 6, fill="#C0392B", width=3)
        self.canvas.create_text(lx + 22, ly, text="% Afectados (V2L < Cut-Off)", anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#C0392B")
        ly += 15
        self.canvas.create_line(lx, ly + 6, lx + 16, ly + 6, fill="#1E8449", width=2, dash=(4, 2))
        self.canvas.create_text(lx + 22, ly, text="% No Afectados (V2R >= Cut-Off)", anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#1E8449")
        ly += 18

        if not self.sum_df.empty:
            s_sub = self.sum_df.copy()
            s_sub['CID'] = s_sub.apply(lambda r: f"{r['Peptide']} | {r['Cell_Line']} | {r['Concentration_uM']:.2f} µM", axis=1)
            row_m = s_sub[s_sub['CID'] == chosen_cid]
            if not row_m.empty:
                r0 = row_m.iloc[0]
                milestones = [
                    ("T_onset stat", float(r0['T_Stat_Onset_s']), "#2874A6", (3, 2)),
                    ("T50 del Efecto", float(r0['T50_Effect_s']), "#D35400", (5, 2)),
                    ("LT50 (50% Pob)", float(r0['LT50_s']), "#6C3483", (2, 2)),
                    ("T_efecto_máx", float(r0['T_Effect_Max_s']), "#922B21", (6, 2))
                ]
                for m_lbl, m_val, m_col, m_dash in milestones:
                    if 0.0 <= m_val < 600.0:
                        px_m, _ = to_px(m_val, 0)
                        self.canvas.create_line(px_m, y0, px_m, y1, fill=m_col, width=2, dash=m_dash)
                    self.canvas.create_text(
                        lx, ly,
                        text=f"• {m_lbl}: {m_val:.1f} s" if m_val < 600.0 else f"• {m_lbl}: >600 s",
                        anchor="nw", font=("DejaVu Sans", 8, "bold"), fill=m_col
                    )
                    ly += 14
                self.canvas.create_text(
                    lx, ly + 2,
                    text=f"Máx V2L: {r0['Max_V2L_Pct']:.1f}% | Mín: {r0['Min_V2L_Pct']:.1f}%\n"
                         f"h Cohen: {r0['Cohen_h_Max']:.2f} | OR: {r0['OR_vs_Basal']:.1f}x\n"
                         f"{r0['Affectation_State']}",
                    anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#1C2833"
                )

    def _draw_milestones_bars(self, W, H):
        if self.sum_df.empty:
            return
        ml, mr, mt, mb = 58, 20, 28, 42
        gap = 55
        pw = max(140, (W - ml - mr - gap) // 2)
        ph = max(90, H - mt - mb)
        b1 = (ml, mt, ml + pw, mt + ph)
        b2 = (ml + pw + gap, mt, ml + 2 * pw + gap, mt + ph)

        self.canvas.create_rectangle(*b1, fill="#FCFDFE", outline="#333333")
        self.canvas.create_rectangle(*b2, fill="#FCFDFE", outline="#333333")
        self.canvas.create_text((b1[0] + b1[2]) // 2, b1[1] - 12, text="A. Hitos Temporales: T_onset (Azul), T50 Efecto (Naranja), T_máx (Rojo)", font=("DejaVu Sans", 8, "bold"), fill="#002B49")
        self.canvas.create_text((b2[0] + b2[2]) // 2, b2[1] - 12, text="B. Proporción en Efecto Máximo: % Afectados V2L (Rojo) vs % No Afectados V2R (Verde)", font=("DejaVu Sans", 8, "bold"), fill="#002B49")

        df_p = self.sum_df.sort_values(['Peptide', 'Concentration_uM']).reset_index(drop=True)
        n_bars = len(df_p)
        if n_bars == 0:
            return

        for gy in [0, 150, 300, 450, 600]:
            py = b1[3] - (gy / 600.0) * ph
            self.canvas.create_line(b1[0], py, b1[2], py, fill="#E5E8E8", dash=(2, 3))
            self.canvas.create_text(b1[0] - 5, py, text=f"{gy}s", anchor="e", font=("DejaVu Sans", 8))

        for gy in [0, 25, 50, 75, 100]:
            py = b2[3] - (gy / 100.0) * ph
            self.canvas.create_line(b2[0], py, b2[2], py, fill="#E5E8E8", dash=(2, 3))
            self.canvas.create_text(b2[0] - 5, py, text=f"{gy}%", anchor="e", font=("DejaVu Sans", 8))

        slot_w = pw / max(1, n_bars)
        bw = max(3.0, min(14.0, slot_w * 0.24))

        for i, r in df_p.iterrows():
            cx1 = b1[0] + (i + 0.5) * slot_w
            cx2 = b2[0] + (i + 0.5) * slot_w
            short_lbl = f"{str(r['Peptide'])[:5]}\n{r['Concentration_uM']:.2f}"

            for offset, val_t, col_b in [(-bw, r['T_Stat_Onset_s'], "#2874A6"), (0.0, r['T50_Effect_s'], "#E67E22"), (bw, r['T_Effect_Max_s'], "#922B21")]:
                vy = b1[3] - (np.clip(float(val_t), 0.0, 600.0) / 600.0) * ph
                self.canvas.create_rectangle(cx1 + offset - bw * 0.45, vy, cx1 + offset + bw * 0.45, b1[3], fill=col_b, outline="")
            self.canvas.create_text(cx1, b1[3] + 15, text=short_lbl, font=("DejaVu Sans", 7), justify="center")

            mv2l = np.clip(float(r['Max_V2L_Pct']), 0.0, 100.0)
            py_split = b2[3] - (mv2l / 100.0) * ph
            self.canvas.create_rectangle(cx2 - bw * 1.1, py_split, cx2 + bw * 1.1, b2[3], fill="#C0392B", outline="white")
            self.canvas.create_rectangle(cx2 - bw * 1.1, b2[1], cx2 + bw * 1.1, py_split, fill="#27AE60", outline="white")
            self.canvas.create_text(cx2, max(b2[1] + 8, py_split - 6), text=f"{mv2l:.0f}%", font=("DejaVu Sans", 7, "bold"), fill="#1C2833")
            self.canvas.create_text(cx2, b2[3] + 15, text=short_lbl, font=("DejaVu Sans", 7), justify="center")

    def _draw_heatmap(self, W, H):
        if self.win_df.empty:
            return
        sub = self.win_df.copy()
        sub['Row_ID'] = sub.apply(lambda r: f"{r['Peptide']} | {r['Cell_Line']} | {r['Concentration_uM']:.2f}µM", axis=1)
        piv = sub.pivot_table(index='Row_ID', columns='Time_Window', values='Mean_Affected_V2L_Pct', aggfunc='mean').fillna(0.0)
        if piv.empty:
            return

        ml, mr, mt, mb = 175, 30, 26, 32
        x0, y0, x1, y1 = ml, mt, max(ml + 180, W - mr), max(mt + 80, H - mb)
        n_rows, n_cols = piv.shape
        cw = (x1 - x0) / max(1, n_cols)
        ch = (y1 - y0) / max(1, n_rows)

        self.canvas.create_text((x0 + x1) // 2, y0 - 12, text="Mapa de Calor: % de Eventos Afectados (V2L) por Ventana de Tiempo", font=("DejaVu Sans", 9, "bold"), fill="#002B49")

        def val_to_hex(pct):
            f = float(np.clip(pct / 100.0, 0.0, 1.0))
            if f < 0.5:
                u = f * 2.0
                r = int(213 * (1 - u) + 249 * u)
                g = int(245 * (1 - u) + 231 * u)
                b = int(227 * (1 - u) + 159 * u)
            else:
                u = (f - 0.5) * 2.0
                r = int(249 * (1 - u) + 192 * u)
                g = int(231 * (1 - u) + 57 * u)
                b = int(159 * (1 - u) + 43 * u)
            return f"#{r:02x}{g:02x}{b:02x}"

        for r_i, r_name in enumerate(piv.index):
            ry0 = y0 + r_i * ch
            self.canvas.create_text(x0 - 6, ry0 + ch * 0.5, text=str(r_name)[:26], anchor="e", font=("DejaVu Sans", 8, "bold"), fill="#1C2833")
            for c_i, t_col in enumerate(piv.columns):
                cx0 = x0 + c_i * cw
                val = float(piv.iloc[r_i, c_i])
                self.canvas.create_rectangle(cx0, ry0, cx0 + cw, ry0 + ch, fill=val_to_hex(val), outline="#F2F4F4")
                if cw >= 22 and ch >= 13:
                    txt_col = "white" if val > 70.0 else "#1C2833"
                    self.canvas.create_text(cx0 + cw * 0.5, ry0 + ch * 0.5, text=f"{val:.0f}", font=("DejaVu Sans", 7), fill=txt_col)

        for c_i, t_col in enumerate(piv.columns):
            if c_i % max(1, n_cols // 12) == 0:
                cx = x0 + (c_i + 0.5) * cw
                self.canvas.create_text(cx, y1 + 11, text=f"{int(t_col)}s", font=("DejaVu Sans", 8))

    def _draw_conversion_rate_and_entropy(self, W, H):
        if self.win_df.empty:
            return
        ml, mr, mt, mb = 68, 220, 28, 34
        x0, y0, x1, y1 = ml, mt, max(ml + 180, W - mr), max(mt + 90, H - mb)
        pw, ph = x1 - x0, y1 - y0
        xmin, xmax = -60.0, 600.0

        self.canvas.create_rectangle(x0, y0, x1, y1, fill="#FCFDFE", outline="#333333")
        self.canvas.create_text((x0 + x1) // 2, y0 - 12, text="Velocidad Instantánea de Conversión d(%V2L)/dt (%/s)", font=("DejaVu Sans", 9, "bold"), fill="#002B49")

        sub = self.win_df.copy()
        sub['CID'] = sub.apply(lambda r: f"{r['Peptide']} | {r['Cell_Line']} | {r['Concentration_uM']:.2f} µM", axis=1)
        unique_cids = sorted(sub['CID'].unique())
        max_rate = max(0.25, float(sub['Inst_Conversion_Rate_pct_s'].max()) * 1.15) if 'Inst_Conversion_Rate_pct_s' in sub else 1.0

        def to_px(tv, rv):
            px = x0 + ((tv - xmin) / (xmax - xmin)) * pw
            py = y1 - (np.clip(rv, 0.0, max_rate) / max_rate) * ph
            return px, py

        for gy in np.linspace(0, max_rate, 5):
            _, py = to_px(0, gy)
            self.canvas.create_line(x0, py, x1, py, fill="#E5E8E8", dash=(2, 3))
            self.canvas.create_text(x0 - 6, py, text=f"{gy:.2f}", anchor="e", font=("DejaVu Sans", 8))

        for tx in [0, 120, 240, 360, 480, 600]:
            px, _ = to_px(tx, 0)
            self.canvas.create_line(px, y0, px, y1, fill="#E5E8E8", dash=(2, 3))
            self.canvas.create_text(px, y1 + 10, text=f"{tx}s", font=("DejaVu Sans", 8))

        lx, ly = x1 + 10, y0 + 4
        self.canvas.create_text(lx, ly, text="Velocidad d(%V2L)/dt:", anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#002B49")
        ly += 16

        for idx_c, cid in enumerate(unique_cids[:8]):
            c_df = sub[sub['CID'] == cid].sort_values('Time_Window')
            col = CURVE_PALETTE[idx_c % len(CURVE_PALETTE)]
            pts = []
            if 'Inst_Conversion_Rate_pct_s' in c_df:
                for tv, rv in zip(c_df['Time_Window'].values, c_df['Inst_Conversion_Rate_pct_s'].values):
                    pts.extend(to_px(float(tv), max(0.0, float(rv))))
                if len(pts) >= 4:
                    self.canvas.create_line(pts, fill=col, width=2)
            self.canvas.create_line(lx, ly + 6, lx + 14, ly + 6, fill=col, width=2)
            self.canvas.create_text(lx + 18, ly, text=cid[:26], anchor="nw", font=("DejaVu Sans", 8), fill="#1C2833")
            ly += 14
