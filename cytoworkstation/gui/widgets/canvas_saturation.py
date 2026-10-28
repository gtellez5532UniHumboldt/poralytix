"""
Lienzo interactivo Tkinter para la vista previa de la curva de saturación y ajuste del modelo matemático.
Muestra la función matemática, parámetros ajustados (Vmax, Kd, nH) y métricas R².
"""
import tkinter as tk
import numpy as np


class NativeSaturationPreview(tk.Frame):
    """
    Lienzo nativo en Tkinter para previsualizar la curva ajustada del modelo farmacodinámico.
    """
    def __init__(self, parent):
        super().__init__(parent, bg="white")
        self.canvas = tk.Canvas(self, bg="white", highlightthickness=0)
        self.canvas.pack(fill="both", expand=True)
        self.curves = []
        self.title = "Previsualización de la Función Matemática y Bondad de Ajuste (R²)"
        self.equation_str = "V_active(C) = (V_max · C^nH) / (K_d^nH + C^nH)"
        self.global_r2_str = ""
        self.x_label = "Concentración (µM)"
        self.y_label = "V_active (%/s)"
        self.canvas.bind("<Configure>", lambda e: self.redraw())

    def set_curves(
        self, curves, title, equation_str="",
        x_label="Concentración (µM)", y_label="V_active (%/s)",
        global_r2_str=""
    ):
        self.curves = curves
        self.title = title
        if equation_str:
            self.equation_str = equation_str
        self.x_label = x_label
        self.y_label = y_label
        self.global_r2_str = global_r2_str
        self.redraw()

    def redraw(self):
        self.canvas.delete("all")
        W = max(360, self.canvas.winfo_width())
        H = max(260, self.canvas.winfo_height())
        ml, mr, mt, mb = 68, 22, 62, 46
        pw = max(140, W - ml - mr)
        ph = max(110, H - mt - mb)
        x0, y0, x1, y1 = ml, mt, ml + pw, mt + ph

        self.canvas.create_rectangle(x0, 5, x1, y0 - 6, fill="#EBF5FB", outline="#85C1E9")
        self.canvas.create_text(
            (x0 + x1) // 2, 16, text=self.title,
            font=("DejaVu Sans", 9, "bold"), fill="#002B49"
        )
        self.canvas.create_text(
            (x0 + x1) // 2, 32, text=f"Función del Modelo: {self.equation_str}",
            font=("Consolas", 8, "bold"), fill="#1B4F72"
        )
        if self.global_r2_str:
            self.canvas.create_text(
                (x0 + x1) // 2, 47, text=self.global_r2_str,
                font=("DejaVu Sans", 8, "bold"), fill="#117864"
            )

        self.canvas.create_rectangle(x0, y0, x1, y1, fill="#FCFDFE", outline="#333333")

        all_x, all_y = [0.0], [0.0]
        for c in self.curves:
            all_x.extend(c.get('x_pts', []))
            all_x.extend(c.get('raw_x', []))
            y_pts = np.asarray(c.get('y_pts', []), dtype=float)
            y_sem = np.asarray(c.get('y_sem', np.zeros_like(y_pts)), dtype=float)
            if len(y_pts) > 0:
                all_y.extend(y_pts + y_sem)
            all_y.extend(c.get('raw_y', []))
            all_x.extend(c.get('x_fit', []))
            all_y.extend(c.get('y_fit', []))

        xmax = max(all_x) * 1.12 if max(all_x) > 0 else 5.0
        ymax = max(all_y) * 1.25 if max(all_y) > 0 else 0.5

        def to_px(xv, yv):
            px = x0 + (xv / xmax) * pw
            py = y1 - (yv / ymax) * ph
            return px, py

        for gx in np.linspace(0, xmax, 6):
            px, _ = to_px(gx, 0)
            self.canvas.create_line(px, y0, px, y1, fill="#E5E8E8", dash=(2, 3))
            lbl_x = f"{gx:.0f}" if xmax >= 50 else f"{gx:.2f}"
            self.canvas.create_text(px, y1 + 12, text=lbl_x, font=("DejaVu Sans", 8))

        for gy in np.linspace(0, ymax, 5):
            _, py = to_px(0, gy)
            self.canvas.create_line(x0, py, x1, py, fill="#E5E8E8", dash=(2, 3))
            lbl_y = f"{gy:.1f}" if ymax >= 5.0 else f"{gy:.3f}"
            self.canvas.create_text(x0 - 6, py, text=lbl_y, anchor="e", font=("DejaVu Sans", 8))

        self.canvas.create_text((x0 + x1) // 2, y1 + 29, text=self.x_label, font=("DejaVu Sans", 8, "bold"), fill="#002B49")
        self.canvas.create_text(x0 - 48, (y0 + y1) // 2, text=self.y_label, font=("DejaVu Sans", 8, "bold"), fill="#002B49")

        for c in self.curves:
            col = c['color']
            dash_pat = c.get('dash', None)

            raw_x = c.get('raw_x', [])
            raw_y = c.get('raw_y', [])
            for rx, ry in zip(raw_x, raw_y):
                prx, pry = to_px(rx, ry)
                self.canvas.create_oval(prx - 2, pry - 2, prx + 2, pry + 2, outline=col, fill="")

            if len(c.get('x_fit', [])) >= 2:
                fit_pts = []
                for xv, yv in zip(c['x_fit'], c['y_fit']):
                    fit_pts.extend(to_px(np.clip(xv, 0.0, xmax), np.clip(yv, 0.0, ymax)))
                if dash_pat:
                    self.canvas.create_line(fit_pts, fill=col, width=2, dash=dash_pat)
                else:
                    self.canvas.create_line(fit_pts, fill=col, width=2)

            kd_val = c.get('kd_marker', None)
            uhalf_val = c.get('y_half', None)
            if kd_val is not None and 0 < kd_val <= xmax and uhalf_val is not None and 0 <= uhalf_val <= ymax:
                pkx, pky = to_px(kd_val, uhalf_val)
                self.canvas.create_line(pkx, pky, pkx, y1, fill=col, dash=(3, 2), width=1)
                self.canvas.create_line(x0, pky, pkx, pky, fill=col, dash=(3, 2), width=1)

            x_pts = c.get('x_pts', [])
            y_pts = c.get('y_pts', [])
            y_sem = c.get('y_sem', [0.0] * len(x_pts))
            for xv, yv, sv in zip(x_pts, y_pts, y_sem):
                px, py = to_px(xv, yv)
                if sv > 0:
                    _, py_top = to_px(xv, min(ymax, yv + sv))
                    _, py_bot = to_px(xv, max(0.0, yv - sv))
                    self.canvas.create_line(px, py_top, px, py_bot, fill=col, width=1)
                    self.canvas.create_line(px - 3, py_top, px + 3, py_top, fill=col, width=1)
                    self.canvas.create_line(px - 3, py_bot, px + 3, py_bot, fill=col, width=1)
                self.canvas.create_oval(px - 4, py - 4, px + 4, py + 4, fill=col, outline="white")

        if self.curves:
            box_h = min(ph - 12, 10 + 28 * len(self.curves))
            self.canvas.create_rectangle(
                x0 + 6, y0 + 6, min(x1 - 6, x0 + max(280, int(pw * 0.78))), y0 + 6 + box_h,
                fill="#FDFEFE", outline="#D6DBDF"
            )
            ly = y0 + 10
            for c in self.curves:
                col = c['color']
                r2_badge = c.get('r2_badge', '')
                eq_eval = c.get('eq_eval', c.get('param_str', ''))
                line1 = f"● {c['label']}   {r2_badge}".strip()
                self.canvas.create_text(x0 + 12, ly, text=line1, anchor="nw", font=("DejaVu Sans", 8, "bold"), fill=col)
                if eq_eval:
                    self.canvas.create_text(x0 + 22, ly + 13, text=eq_eval, anchor="nw", font=("Consolas", 8), fill="#1C2833")
                ly += 27
