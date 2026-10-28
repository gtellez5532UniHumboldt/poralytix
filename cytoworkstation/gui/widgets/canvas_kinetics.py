"""
Lienzo científico interactivo Tkinter para cinéticas de permeabilización, FSC y SSC.
Permite arrastrar con el ratón las líneas verticales T_lag y T_max directamente sobre la gráfica.
"""
import tkinter as tk
import numpy as np


class NativeInteractiveMultiPlot(tk.Frame):
    """
    Lienzo científico interactivo en Tkinter.
    Dibuja simultáneamente:
      1. Permeabilización (%) o Afectados V2L (%) vs Tiempo con T_lag y T_max arrastrables.
      2. FSC vs Tiempo (tamaño celular / tumefacción).
      3. SSC vs Tiempo (granularidad interna).
    """
    def __init__(self, parent, on_threshold_drag_callback=None):
        super().__init__(parent, bg="white")
        self.on_threshold_drag_callback = on_threshold_drag_callback

        self.canvas = tk.Canvas(self, bg="white", highlightthickness=0, cursor="crosshair")
        self.canvas.pack(fill="both", expand=True)

        self.series_data = []
        self.active_tlag = None
        self.active_tmax = None
        self.active_title = ""
        self.scatter_norm = False
        self.display_mode = "Continua (Mediana FITC %)"
        self.dragging = None

        self.xmin, self.xmax = -60.0, 600.0
        self.plot_boxes = []

        self.canvas.bind("<Configure>", lambda e: self.redraw())
        self.canvas.bind("<ButtonPress-1>", self._on_press)
        self.canvas.bind("<B1-Motion>", self._on_motion)
        self.canvas.bind("<ButtonRelease-1>", self._on_release)

    def set_data(self, series_data, active_tlag, active_tmax, active_title, scatter_norm=False, display_mode="Continua (Mediana FITC %)"):
        self.series_data = series_data
        self.active_tlag = active_tlag
        self.active_tmax = active_tmax
        self.active_title = active_title
        self.scatter_norm = scatter_norm
        self.display_mode = display_mode
        self.redraw()

    def _x_to_px(self, x_val, box):
        x0, _, x1, _ = box
        return x0 + ((x_val - self.xmin) / (self.xmax - self.xmin)) * max(1.0, (x1 - x0))

    def _px_to_x(self, px, box):
        x0, _, x1, _ = box
        frac = (px - x0) / max(1.0, (x1 - x0))
        return np.clip(self.xmin + frac * (self.xmax - self.xmin), 0.0, 600.0)

    def _y_to_px(self, y_val, ymin, ymax, box):
        _, y0, _, y1 = box
        frac = (y_val - ymin) / max(1e-9, (ymax - ymin))
        return y1 - np.clip(frac, 0.0, 1.0) * max(1.0, (y1 - y0))

    def redraw(self):
        self.canvas.delete("all")
        W = max(480, self.canvas.winfo_width())
        H = max(360, self.canvas.winfo_height())

        ml = 76
        mr = 235
        mt = 28
        mb = 38
        gap = 38
        pw = max(200, W - ml - mr)
        ph = max(90, (H - mt - mb - 2 * gap) // 3)

        box1 = (ml, mt, ml + pw, mt + ph)
        box2 = (ml, mt + ph + gap, ml + pw, mt + 2 * ph + gap)
        box3 = (ml, mt + 2 * (ph + gap), ml + pw, mt + 3 * ph + 2 * gap)
        self.plot_boxes = [box1, box2, box3]

        fsc_vals, ssc_vals = [], []
        for s in self.series_data:
            fsc_vals.extend(s['fsc_m'])
            ssc_vals.extend(s['ssc_m'])

        fsc_min = min(fsc_vals) * 0.92 if fsc_vals else 0.0
        fsc_max = max(fsc_vals) * 1.08 if fsc_vals else 100.0
        if np.isclose(fsc_min, fsc_max):
            fsc_max = fsc_min + 10.0

        ssc_min = min(ssc_vals) * 0.92 if ssc_vals else 0.0
        ssc_max = max(ssc_vals) * 1.08 if ssc_vals else 100.0
        if np.isclose(ssc_min, ssc_max):
            ssc_max = ssc_min + 10.0

        if "Categórica" in self.display_mode:
            p1_title = f"1. Proporción de Eventos Afectados (V2L % < Cut-Off) vs. Tiempo — {self.active_title}"
            p1_ylabel = "Afectados V2L (%)"
        elif "Dual" in self.display_mode:
            p1_title = f"1. Superposición Dual: Continua (Sólida) vs. Afectados V2L % (Punteada) — {self.active_title}"
            p1_ylabel = "Permeab / V2L (%)"
        else:
            p1_title = f"1. Cinética de Permeabilización (Mediana Continua) vs. Tiempo — {self.active_title}"
            p1_ylabel = "Permeab. (%)"

        configs = [
            (box1, -5.0, 105.0, p1_title, p1_ylabel),
            (box2, fsc_min, fsc_max, "2. Dispersión Frontal (FSC vs. Tiempo — Tamaño Celular / Tumefacción)", "FSC (%)" if self.scatter_norm else "Mean FSC-A"),
            (box3, ssc_min, ssc_max, "3. Dispersión Lateral (SSC vs. Tiempo — Granularidad / Complejidad Interna)", "SSC (%)" if self.scatter_norm else "Mean SSC-A")
        ]

        for idx_p, (box, ymin, ymax, title, ylabel) in enumerate(configs):
            x0, y0, x1, y1 = box
            self.canvas.create_rectangle(x0, y0, x1, y1, fill="#FCFDFE", outline="#333333", width=1)
            self.canvas.create_text((x0 + x1) // 2, y0 - 11, text=title, font=("DejaVu Sans", 9, "bold"), fill="#002B49")
            self.canvas.create_text(x0 - 44, (y0 + y1) // 2, text=ylabel, font=("DejaVu Sans", 8, "bold"), fill="#1C2833")

            for tx in [0, 100, 200, 300, 400, 500, 600]:
                px = self._x_to_px(tx, box)
                self.canvas.create_line(px, y0, px, y1, fill="#E5E8E8", dash=(2, 3))
                if idx_p == 2:
                    self.canvas.create_text(px, y1 + 11, text=str(tx), font=("DejaVu Sans", 8), fill="#333333")

            px0 = self._x_to_px(0.0, box)
            self.canvas.create_line(px0, y0, px0, y1, fill="#333333", dash=(4, 3), width=1)

            for gy in np.linspace(ymin, ymax, 5):
                py = self._y_to_px(gy, ymin, ymax, box)
                self.canvas.create_line(x0, py, x1, py, fill="#E5E8E8", dash=(2, 3))
                lbl_y = f"{gy:.0f}" if abs(ymax - ymin) > 25 else f"{gy:.1f}"
                self.canvas.create_text(x0 - 6, py, text=lbl_y, anchor="e", font=("DejaVu Sans", 8), fill="#333333")

        self.canvas.create_text((box3[0] + box3[2]) // 2, box3[3] + 25, text="Tiempo (segundos)", font=("DejaVu Sans", 9, "bold"), fill="#002B49")

        if self.active_tlag is not None and self.active_tmax is not None:
            px_lag = self._x_to_px(self.active_tlag, box1)
            px_max = self._x_to_px(self.active_tmax, box1)
            self.canvas.create_rectangle(
                px_lag, box1[1] + 1, px_max, box1[3] - 1,
                fill="#EBF5FB", outline=""
            )

        def _soft_fill(hex_c):
            h = hex_c.lstrip('#')
            if len(h) != 6:
                return "#EBF5FB"
            r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
            return f"#{int(r*0.20 + 255*0.80):02x}{int(g*0.20 + 255*0.80):02x}{int(b*0.20 + 255*0.80):02x}"

        for s in self.series_data:
            t_arr = s['t']
            col = s['color']
            fill_c = _soft_fill(col)
            is_ctrl = s['is_ctrl']
            dash_pat = (3, 3) if is_ctrl else None

            primary_y_key = 'v2l_m' if "Categórica" in self.display_mode else 'perm_m'
            primary_se_key = 'v2l_se' if "Categórica" in self.display_mode else 'perm_se'

            for box, y_key, se_key, ymin, ymax in [
                (box1, primary_y_key, primary_se_key, -5.0, 105.0),
                (box2, 'fsc_m', 'fsc_se', fsc_min, fsc_max),
                (box3, 'ssc_m', 'ssc_se', ssc_min, ssc_max)
            ]:
                y_arr = s[y_key]
                se_arr = s[se_key]
                if len(t_arr) < 2:
                    continue

                poly_pts = []
                for t_v, y_v, se_v in zip(t_arr, y_arr, se_arr):
                    poly_pts.extend([self._x_to_px(t_v, box), self._y_to_px(y_v + se_v, ymin, ymax, box)])
                for t_v, y_v, se_v in zip(reversed(t_arr), reversed(y_arr), reversed(se_arr)):
                    poly_pts.extend([self._x_to_px(t_v, box), self._y_to_px(y_v - se_v, ymin, ymax, box)])
                if len(poly_pts) >= 6:
                    self.canvas.create_polygon(poly_pts, fill=fill_c, outline="")

                line_pts = []
                for t_v, y_v in zip(t_arr, y_arr):
                    px = self._x_to_px(t_v, box)
                    py = self._y_to_px(y_v, ymin, ymax, box)
                    line_pts.extend([px, py])
                    if is_ctrl:
                        self.canvas.create_oval(px - 2.5, py - 2.5, px + 2.5, py + 2.5, fill=col, outline="white")
                    else:
                        self.canvas.create_rectangle(px - 2.5, py - 2.5, px + 2.5, py + 2.5, fill=col, outline=col)

                if len(line_pts) >= 4:
                    self.canvas.create_line(line_pts, fill=col, width=2, dash=dash_pat)

            if "Dual" in self.display_mode and len(t_arr) >= 2 and 'v2l_m' in s:
                v2l_pts = []
                for t_v, y_v in zip(t_arr, s['v2l_m']):
                    v2l_pts.extend([self._x_to_px(t_v, box1), self._y_to_px(y_v, -5.0, 105.0, box1)])
                if len(v2l_pts) >= 4:
                    self.canvas.create_line(v2l_pts, fill=col, width=2, dash=(6, 3))

        if self.active_tlag is not None and self.active_tmax is not None:
            px_lag = self._x_to_px(self.active_tlag, box1)
            px_max = self._x_to_px(self.active_tmax, box1)
            self.canvas.create_line(px_lag, box1[1], px_lag, box1[3], fill="#0055CC", width=2, dash=(4, 2))
            self.canvas.create_rectangle(px_lag - 5, box1[1] + 4, px_lag + 5, box1[1] + 18, fill="#0055CC", outline="white")
            self.canvas.create_text(px_lag + 6, box1[1] + 11, text=f"T_lag ({self.active_tlag:.0f}s)", anchor="w", font=("DejaVu Sans", 8, "bold"), fill="#0055CC")

            self.canvas.create_line(px_max, box1[1], px_max, box1[3], fill="#CC2200", width=2, dash=(6, 2, 2, 2))
            self.canvas.create_rectangle(px_max - 5, box1[1] + 22, px_max + 5, box1[1] + 36, fill="#CC2200", outline="white")
            self.canvas.create_text(px_max + 6, box1[1] + 29, text=f"T_max ({self.active_tmax:.0f}s)", anchor="w", font=("DejaVu Sans", 8, "bold"), fill="#CC2200")

        lx = box1[2] + 12
        ly = box1[1] + 4
        self.canvas.create_text(lx, ly, text="Leyenda de Curvas:", anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#002B49")
        ly += 16
        for s in self.series_data[:12]:
            col = s['color']
            dash_pat = (3, 2) if s['is_ctrl'] else None
            self.canvas.create_line(lx, ly + 6, lx + 18, ly + 6, fill=col, width=2, dash=dash_pat)
            self.canvas.create_text(lx + 24, ly, text=s['label'], anchor="nw", font=("DejaVu Sans", 8), fill="#1C2833")
            ly += 16

        if self.active_tlag is not None:
            ly += 6
            self.canvas.create_line(lx, ly + 6, lx + 18, ly + 6, fill="#0055CC", width=2, dash=(3, 2))
            self.canvas.create_text(lx + 24, ly, text=f"T_lag ({self.active_tlag:.0f}s) [Arrastrable]", anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#0055CC")
            ly += 16
            self.canvas.create_line(lx, ly + 6, lx + 18, ly + 6, fill="#CC2200", width=2, dash=(5, 2))
            self.canvas.create_text(lx + 24, ly, text=f"T_max ({self.active_tmax:.0f}s) [Arrastrable]", anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#CC2200")

    def _on_press(self, event):
        if not self.plot_boxes or self.active_tlag is None or self.active_tmax is None:
            return
        box1 = self.plot_boxes[0]
        if not (box1[0] <= event.x <= box1[2] and box1[1] <= event.y <= box1[3]):
            return
        px_lag = self._x_to_px(self.active_tlag, box1)
        px_max = self._x_to_px(self.active_tmax, box1)
        d_lag = abs(event.x - px_lag)
        d_max = abs(event.x - px_max)
        if min(d_lag, d_max) <= 18:
            self.dragging = 'tlag' if d_lag <= d_max else 'tmax'

    def _on_motion(self, event):
        if not self.dragging or not self.plot_boxes:
            return
        box1 = self.plot_boxes[0]
        new_t = round(self._px_to_x(event.x, box1) / 15.0) * 15.0
        if self.dragging == 'tlag':
            self.active_tlag = min(new_t, self.active_tmax - 15.0)
        else:
            self.active_tmax = max(new_t, self.active_tlag + 15.0)
        self.redraw()
        if self.on_threshold_drag_callback:
            self.on_threshold_drag_callback(self.active_tlag, self.active_tmax, final=False)

    def _on_release(self, event):
        if self.dragging:
            self.dragging = None
            if self.on_threshold_drag_callback:
                self.on_threshold_drag_callback(self.active_tlag, self.active_tmax, final=True)
