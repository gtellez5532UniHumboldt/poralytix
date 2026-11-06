"""
Lienzo interactivo Tkinter para Visualización Categórica (Barras Apiladas, Hitos, Heatmap y Tasa de Conversión).
Soporte bilingüe dinámico.
"""
import tkinter as tk
import numpy as np
import pandas as pd

from ...i18n import t


class NativeCategoricalProportionsCanvas(tk.Frame):
    """
    Lienzo nativo Tkinter con 4 Modos de visualización categórica intercambiables:
      - Modo 1: Balance Temporal (% Afectados V2L vs % No Afectados V2R con hitos T_onset, T50 y T_max).
      - Modo 2: Comparativa de Hitos Temporales por Dosis (Gráfico de barras de T_onset, T50, T_max).
      - Modo 3: Mapa de Calor Espaciotemporal (% V2L por Compuesto × Dosis × Ventana de Tiempo).
      - Modo 4: Velocidad Instantánea de Conversión d(%V2L)/dt y Entropía H(t).
    """
    def __init__(self, parent):
        super().__init__(parent, bg="white")
        self.canvas = tk.Canvas(self, bg="white", highlightthickness=0)
        self.canvas.pack(fill="both", expand=True)

        self.mode_str = "Balance Temporal"
        self.win_df = pd.DataFrame()
        self.sum_df = pd.DataFrame()
        self.active_label = ""

        self.canvas.bind("<Configure>", lambda e: self.redraw())

    def set_panel_data(self, mode_str, win_df, sum_df, active_label=""):
        self.mode_str = str(mode_str)
        self.win_df = win_df.copy() if win_df is not None else pd.DataFrame()
        self.sum_df = sum_df.copy() if sum_df is not None else pd.DataFrame()
        self.active_label = active_label
        self.redraw()

    def redraw(self):
        self.canvas.delete("all")
        W = max(520, self.canvas.winfo_width())
        H = max(260, self.canvas.winfo_height())

        if "Hitos" in self.mode_str or "Milestone" in self.mode_str:
            self._draw_milestones_comparison(W, H)
        elif "Calor" in self.mode_str or "Heatmap" in self.mode_str:
            self._draw_heatmap(W, H)
        elif "Velocidad" in self.mode_str or "Instantaneous" in self.mode_str:
            self._draw_conversion_rate_and_entropy(W, H)
        else:
            self._draw_temporal_balance(W, H)

    def _draw_temporal_balance(self, W, H):
        if self.win_df.empty:
            self.canvas.create_text(
                W // 2, H // 2,
                text="Sin datos categóricos calculados. Ejecuta el Pre-Análisis.",
                font=("DejaVu Sans", 9, "italic"), fill="#7F8C8D"
            )
            return

        ml, mr, mt, mb = 65, 230, 26, 32
        x0, y0, x1, y1 = ml, mt, max(ml + 180, W - mr), max(mt + 80, H - mb)
        pw, ph = x1 - x0, y1 - y0

        self.canvas.create_rectangle(x0, y0, x1, y1, fill="#FCFDFE", outline="#333333")

        t_min, t_max = -60.0, 600.0

        def to_px(tv, pct):
            px = x0 + ((tv - t_min) / (t_max - t_min)) * pw
            py = y1 - (np.clip(pct, 0.0, 100.0) / 100.0) * ph
            return px, py

        cids = sorted(self.win_df['Peptide'].astype(str) + " | " + self.win_df['Cell_Line'].astype(str) + " | " + self.win_df['Concentration_uM'].apply(lambda v: f"{v:.2f} µM"))
        chosen_cid = self.active_label if (self.active_label and self.active_label in cids) else (cids[0] if cids else "")

        title_txt = t("canvas_cat_curve").format(chosen_cid)
        self.canvas.create_text((x0 + x1) // 2, y0 - 12, text=title_txt, font=("DejaVu Sans", 9, "bold"), fill="#002B49")

        sub = self.win_df.copy()
        sub['CID'] = sub['Peptide'].astype(str) + " | " + sub['Cell_Line'].astype(str) + " | " + sub['Concentration_uM'].apply(lambda v: f"{v:.2f} µM")
        curve_data = sub[sub['CID'] == chosen_cid].sort_values('Time_Window').reset_index(drop=True)

        if not curve_data.empty:
            v2l_pts, v2r_pts = [], []
            for _, r in curve_data.iterrows():
                tw = float(r['Time_Window'])
                v2l_pct = float(r['Mean_Affected_V2L_Pct'])
                v2r_pct = float(r['Mean_Intact_V2R_Pct'])
                v2l_pts.append(to_px(tw, v2l_pct))
                v2r_pts.append(to_px(tw, v2r_pct))

            poly_pts = [to_px(curve_data.iloc[0]['Time_Window'], 0)[0], y1]
            for px, py in v2l_pts:
                poly_pts.extend([px, py])
            poly_pts.extend([to_px(curve_data.iloc[-1]['Time_Window'], 0)[0], y1])
            if len(poly_pts) >= 6:
                self.canvas.create_polygon(poly_pts, fill="#FDEDEC", outline="")

            line_v2l = []
            for px, py in v2l_pts:
                line_v2l.extend([px, py])
            if len(line_v2l) >= 4:
                self.canvas.create_line(line_v2l, fill="#C0392B", width=2.5)

            line_v2r = []
            for px, py in v2r_pts:
                line_v2r.extend([px, py])
            if len(line_v2r) >= 4:
                self.canvas.create_line(line_v2r, fill="#1E8449", width=2, dash=(4, 2))

            for _, r in curve_data.iterrows():
                tw = float(r['Time_Window'])
                v2l_pct = float(r['Mean_Affected_V2L_Pct'])
                v2r_pct = float(r['Mean_Intact_V2R_Pct'])
                px, py_a = to_px(tw, v2l_pct)
                _, py_r = to_px(tw, v2r_pct)
                self.canvas.create_oval(px - 3, py_a - 3, px + 3, py_a + 3, fill="#C0392B", outline="white")
                self.canvas.create_oval(px - 3, py_r - 3, px + 3, py_r + 3, fill="#1E8449", outline="white")
                if float(r.get('p_value_Fisher', 1.0)) < 0.05:
                    self.canvas.create_text(px, max(y0 + 6, py_a - 7), text="*", font=("DejaVu Sans", 8, "bold"), fill="#922B21")

        for tx in [-60, 0, 60, 120, 180, 240, 300, 360, 420, 480, 540, 600]:
            px, _ = to_px(tx, 0)
            self.canvas.create_line(px, y0, px, y1, fill="#D5D8DC" if tx != 0 else "#333333", dash=(2, 3) if tx != 0 else (4, 2))
            self.canvas.create_text(px, y1 + 10, text=f"{tx}s", font=("DejaVu Sans", 8))

        for gy in [0, 25, 50, 75, 100]:
            _, py = to_px(0, gy)
            self.canvas.create_line(x0, py, x1, py, fill="#BFC9CA" if gy == 50 else "#E5E8E8", dash=(4, 2) if gy == 50 else (2, 3))
            self.canvas.create_text(x0 - 6, py, text=f"{gy}%", anchor="e", font=("DejaVu Sans", 8))

        self.canvas.create_text(x0 - 44, (y0 + y1) // 2, text=t("canvas_cat_ylabel"), font=("DejaVu Sans", 8, "bold"), fill="#1C2833")
        self.canvas.create_text((x0 + x1) // 2, y1 + 23, text=t("canvas_cat_xlabel"), font=("DejaVu Sans", 8, "bold"), fill="#002B49")

        lx, ly = x1 + 10, y0 + 2
        self.canvas.create_text(lx, ly, text=t("canvas_cat_curve").format(chosen_cid[:28]), anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#002B49")
        ly += 16
        self.canvas.create_line(lx, ly + 6, lx + 16, ly + 6, fill="#C0392B", width=3)
        self.canvas.create_text(lx + 22, ly, text=t("canvas_cat_v2l"), anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#C0392B")
        ly += 15
        self.canvas.create_line(lx, ly + 6, lx + 16, ly + 6, fill="#1E8449", width=2, dash=(4, 2))
        self.canvas.create_text(lx + 22, ly, text=t("canvas_cat_v2r"), anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#1E8449")
        ly += 18

        if not self.sum_df.empty:
            s_sub = self.sum_df.copy()
            s_sub['CID'] = s_sub['Peptide'].astype(str) + " | " + s_sub['Cell_Line'].astype(str) + " | " + s_sub['Concentration_uM'].apply(lambda v: f"{v:.2f} µM")
            hit = s_sub[s_sub['CID'] == chosen_cid]
            if not hit.empty:
                r = hit.iloc[0]
                t_onset = r.get('T_Stat_Onset_s', np.nan)
                t50 = r.get('T50_Effect_s', np.nan)
                lt50 = r.get('LT50_Population_s', np.nan)
                t_max_eff = r.get('T_Effect_Max_s', np.nan)

                markers = [
                    (t_onset, "#2874A6", "T_onset", 1),
                    (t50, "#E67E22", "T50 Efecto", 2),
                    (lt50, "#884EA0", "LT50 Abs", 3),
                    (t_max_eff, "#922B21", "T_máx", 4)
                ]
                for tv, col, tag, rank in markers:
                    if tv is not None and np.isfinite(tv) and t_min <= tv <= t_max:
                        px, _ = to_px(tv, 0)
                        self.canvas.create_line(px, y0, px, y1, fill=col, width=1.5, dash=(4, 2))
                        self.canvas.create_polygon([px, y0, px - 4, y0 - 8, px + 4, y0 - 8], fill=col, outline="")
                        self.canvas.create_text(
                            lx, ly, text=f"• {tag}: {tv:.0f}s",
                            anchor="nw", font=("DejaVu Sans", 8, "bold"), fill=col
                        )
                        ly += 14

                cat_txt = str(r.get('Affectation_Status', 'Indeterminado'))
                ly += 4
                self.canvas.create_text(
                    lx, ly, text=f"Estado: {cat_txt}",
                    anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#1C2833"
                )

    def _draw_milestones_comparison(self, W, H):
        if self.sum_df.empty:
            return
        ml, mr, mt, mb = 65, 20, 26, 44
        gap = 25
        avail_w = max(240, W - ml - mr - gap)
        pw = avail_w / 2.0
        ph = max(90, H - mt - mb)
        b1 = (ml, mt, ml + pw, mt + ph)
        b2 = (ml + pw + gap, mt, ml + 2 * pw + gap, mt + ph)

        self.canvas.create_rectangle(*b1, fill="#FCFDFE", outline="#333333")
        self.canvas.create_rectangle(*b2, fill="#FCFDFE", outline="#333333")
        self.canvas.create_text((b1[0] + b1[2]) // 2, b1[1] - 12, text=t("canvas_cat_milestones_a"), font=("DejaVu Sans", 8, "bold"), fill="#002B49")
        self.canvas.create_text((b2[0] + b2[2]) // 2, b2[1] - 12, text=t("canvas_cat_milestones_b"), font=("DejaVu Sans", 8, "bold"), fill="#002B49")

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

        self.canvas.create_text((x0 + x1) // 2, y0 - 12, text=t("canvas_cat_heatmap_title"), font=("DejaVu Sans", 9, "bold"), fill="#002B49")

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
        self.canvas.create_text((x0 + x1) // 2, y0 - 12, text=t("canvas_cat_rate_title"), font=("DejaVu Sans", 9, "bold"), fill="#002B49")

        sub = self.win_df.copy()
        sub['CID'] = sub.apply(lambda r: f"{r['Peptide']} | {r['Cell_Line']} | {r['Concentration_uM']:.2f} µM", axis=1)
        unique_cids = sorted(sub['CID'].unique())
        max_rate = max(0.25, float(sub['Inst_Conversion_Rate_pct_s'].max()) * 1.15) if 'Inst_Conversion_Rate_pct_s' in sub else 1.0

        def to_px(tv, rv):
            px = x0 + ((tv - xmin) / (xmax - xmin)) * pw
            py = y1 - (np.clip(rv, 0.0, max_rate) / max_rate) * ph
            return px, py

        for tx in [0, 100, 200, 300, 400, 500, 600]:
            px, _ = to_px(tx, 0)
            self.canvas.create_line(px, y0, px, y1, fill="#E5E8E8", dash=(2, 3))
            self.canvas.create_text(px, y1 + 10, text=f"{tx}s", font=("DejaVu Sans", 8))

        for gy in np.linspace(0.0, max_rate, 5):
            _, py = to_px(0, gy)
            self.canvas.create_line(x0, py, x1, py, fill="#E5E8E8", dash=(2, 3))
            self.canvas.create_text(x0 - 6, py, text=f"{gy:.2f}", anchor="e", font=("DejaVu Sans", 8))

        palette = ["#1B4F72", "#2874A6", "#117864", "#D35400", "#922B21", "#6C3483", "#2E4053"]
        lx, ly = x1 + 10, y0 + 4
        self.canvas.create_text(lx, ly, text=t("canvas_cat_rate_legend"), anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#002B49")
        ly += 16

        for i, cid in enumerate(unique_cids[:7]):
            col = palette[i % len(palette)]
            c_data = sub[sub['CID'] == cid].sort_values('Time_Window')
            if len(c_data) < 2 or 'Inst_Conversion_Rate_pct_s' not in c_data:
                continue

            pts = []
            for _, r in c_data.iterrows():
                tw = float(r['Time_Window'])
                rate_v = float(r['Inst_Conversion_Rate_pct_s'])
                pts.extend(to_px(tw, rate_v))

            if len(pts) >= 4:
                self.canvas.create_line(pts, fill=col, width=2)
            for px, py in zip(pts[0::2], pts[1::2]):
                self.canvas.create_oval(px - 2.5, py - 2.5, px + 2.5, py + 2.5, fill=col, outline="white")

            self.canvas.create_line(lx, ly + 6, lx + 14, ly + 6, fill=col, width=2)
            self.canvas.create_text(lx + 18, ly, text=cid[:26], anchor="nw", font=("DejaVu Sans", 8), fill="#1C2833")
            ly += 15
