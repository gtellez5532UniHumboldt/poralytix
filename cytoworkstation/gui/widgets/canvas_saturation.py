"""
Lienzo interactivo Tkinter para la Previsualización de Saturación Farmacodinámica y Ajuste No Lineal.
Soporte bilingüe dinámico.
"""
import tkinter as tk
import numpy as np

from ...i18n import t


class NativeSaturationPreview(tk.Frame):
    """
    Lienzo nativo Tkinter para visualizar en tiempo real:
      - Curvas ajustadas de velocidad V_active vs Concentración (Hill, MM, Gompertz, Weibull).
      - Puntos experimentales por réplica y media ± SEM.
      - Línea vertical K_d / EC50 y asíntota V_max.
    """
    def __init__(self, parent):
        super().__init__(parent, bg="white")
        self.canvas = tk.Canvas(self, bg="white", highlightthickness=0)
        self.canvas.pack(fill="both", expand=True)

        self.models_data = []  # [{'name', 'color', 'x_fit', 'y_fit', 'vmax', 'kd', 'nh', 'r2', 'is_active'}, ...]
        self.points_data = []  # [{'conc', 'mean', 'sem', 'reps': [v1, v2, ...]}, ...]
        self.title = ""
        self.equation_str = ""
        self.global_r2_str = ""
        self.x_label = "Concentración (µM)"
        self.y_label = "V_active (%/s)"

        self.canvas.bind("<Configure>", lambda e: self.redraw())

    def set_preview_data(self, models_data, points_data, title="", equation_str="", global_r2_str="", x_label=None, y_label=None):
        self.models_data = models_data
        self.points_data = points_data
        self.title = title
        self.equation_str = equation_str
        self.global_r2_str = global_r2_str
        if x_label:
            self.x_label = x_label
        if y_label:
            self.y_label = y_label
        self.redraw()

    def redraw(self):
        self.canvas.delete("all")
        W = max(380, self.canvas.winfo_width())
        H = max(260, self.canvas.winfo_height())

        ml, mr, mt, mb = 64, 25, 54, 42
        pw = max(180, W - ml - mr)
        ph = max(110, H - mt - mb)
        x0, y0, x1, y1 = ml, mt, ml + pw, mt + ph

        self.canvas.create_rectangle(x0, 5, x1, y0 - 6, fill="#EBF5FB", outline="#85C1E9")
        self.canvas.create_text(
            (x0 + x1) // 2, 16, text=self.title,
            font=("DejaVu Sans", 9, "bold"), fill="#002B49"
        )
        self.canvas.create_text(
            (x0 + x1) // 2, 32, text=t("canvas_sat_fn").format(self.equation_str),
            font=("Consolas", 8, "bold"), fill="#1B4F72"
        )
        if self.global_r2_str:
            self.canvas.create_text(
                (x0 + x1) // 2, 47, text=self.global_r2_str,
                font=("DejaVu Sans", 8, "bold"), fill="#117864"
            )

        self.canvas.create_rectangle(x0, y0, x1, y1, fill="#FCFDFE", outline="#333333")

        all_x, all_y = [0.0], [0.0]
        for p in self.points_data:
            all_x.append(float(p['conc']))
            all_y.append(float(p['mean'] + p.get('sem', 0.0)))
            for r_v in p.get('reps', []):
                all_y.append(float(r_v))
        for m in self.models_data:
            if len(m['x_fit']) > 0:
                all_x.extend([float(np.min(m['x_fit'])), float(np.max(m['x_fit']))])
            if len(m['y_fit']) > 0:
                all_y.extend([float(np.min(m['y_fit'])), float(np.max(m['y_fit']))])

        xmax = max(1.0, float(np.max(all_x)) * 1.15)
        ymax = max(0.1, float(np.max(all_y)) * 1.20)

        def to_px(xv, yv):
            px = x0 + (np.clip(float(xv), 0.0, xmax) / xmax) * pw
            py = y1 - (np.clip(float(yv), 0.0, ymax) / ymax) * ph
            return px, py

        for gx in np.linspace(0.0, xmax, 6):
            px, _ = to_px(gx, 0)
            self.canvas.create_line(px, y0, px, y1, fill="#E5E8E8", dash=(2, 3))
            lbl_x = f"{gx:.1f}" if xmax < 15 else f"{gx:.0f}"
            self.canvas.create_text(px, y1 + 12, text=lbl_x, font=("DejaVu Sans", 8))

        for gy in np.linspace(0.0, ymax, 5):
            _, py = to_px(0, gy)
            self.canvas.create_line(x0, py, x1, py, fill="#E5E8E8", dash=(2, 3))
            lbl_y = f"{gy:.2f}" if ymax < 2 else f"{gy:.1f}"
            self.canvas.create_text(x0 - 6, py, text=lbl_y, anchor="e", font=("DejaVu Sans", 8))

        self.canvas.create_text((x0 + x1) // 2, y1 + 29, text=self.x_label, font=("DejaVu Sans", 8, "bold"), fill="#002B49")
        self.canvas.create_text(x0 - 48, (y0 + y1) // 2, text=self.y_label, font=("DejaVu Sans", 8, "bold"), fill="#002B49")

        for m in self.models_data:
            xf, yf = m['x_fit'], m['y_fit']
            if len(xf) < 2:
                continue
            pts = []
            for xv, yv in zip(xf, yf):
                pts.extend(to_px(xv, yv))
            col = m['color']
            w_px = 3 if m.get('is_active', False) else 1
            dash_p = None if m.get('is_active', False) else (4, 3)
            self.canvas.create_line(pts, fill=col, width=w_px, dash=dash_p)

            kd_v = m.get('kd')
            if m.get('is_active', False) and kd_v is not None and np.isfinite(kd_v) and 0 < kd_v < xmax:
                px_kd, _ = to_px(kd_v, 0)
                self.canvas.create_line(px_kd, y0, px_kd, y1, fill="#E67E22", dash=(3, 2), width=1)
                self.canvas.create_text(px_kd + 4, y0 + 12, text=f"Kd = {kd_v:.2f} µM", font=("DejaVu Sans", 7, "bold"), fill="#E67E22", anchor="w")

        for p in self.points_data:
            cx = float(p['conc'])
            my = float(p['mean'])
            sem = float(p.get('sem', 0.0))

            px, py = to_px(cx, my)
            _, py_top = to_px(cx, my + sem)
            _, py_bot = to_px(cx, max(0.0, my - sem))

            self.canvas.create_line(px, py_bot, px, py_top, fill="#1C2833", width=1)
            self.canvas.create_line(px - 3, py_top, px + 3, py_top, fill="#1C2833", width=1)
            self.canvas.create_line(px - 3, py_bot, px + 3, py_bot, fill="#1C2833", width=1)

            for r_v in p.get('reps', []):
                _, r_py = to_px(cx, r_v)
                self.canvas.create_oval(px - 2.5, r_py - 2.5, px + 2.5, r_py + 2.5, fill="#AED6F1", outline="#2874A6")

            self.canvas.create_oval(px - 4.5, py - 4.5, px + 4.5, py + 4.5, fill="#1B4F72", outline="white", width=1.5)

        ly = y0 + 6
        for m in self.models_data:
            col = m['color']
            name = m['name']
            r2_val = m.get('r2')
            r2_txt = f" (R² = {r2_val:.3f})" if r2_val is not None and np.isfinite(r2_val) else ""
            line1 = f"{name}{r2_txt}"
            eq_eval = m.get('eq_eval', '')

            self.canvas.create_line(x0 + 8, ly + 6, x0 + 20, ly + 6, fill=col, width=3 if m.get('is_active') else 1)
            self.canvas.create_text(x0 + 12, ly, text=line1, anchor="nw", font=("DejaVu Sans", 8, "bold"), fill=col)
            if eq_eval and m.get('is_active'):
                self.canvas.create_text(x0 + 22, ly + 13, text=eq_eval, anchor="nw", font=("Consolas", 8), fill="#1C2833")
                ly += 13
            ly += 15
