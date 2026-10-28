"""
Lienzo interactivo Tkinter para Histograma de Fluorescencia y Cut-Offs arrastrables con el ratón.
"""
import tkinter as tk
import numpy as np


class NativeHistogramCutoffCanvas(tk.Frame):
    """
    Lienzo interactivo en Tkinter para visualizar el Histograma de Fluorescencia
    (X = log10 Fluorescencia, Y = Número de eventos o %) y arrastrar los 2 Cut-Offs:
      - Cut-Off 1 (Exclusión V1L, línea morada): Omite eventos < excl_log.
      - Cut-Off 2 (Discriminación V2L/V2R, línea roja): Clasifica en Afectados vs No Afectados.
    """
    def __init__(self, parent, on_cutoff_change_callback=None):
        super().__init__(parent, bg="white")
        self.on_cutoff_change_callback = on_cutoff_change_callback
        self.canvas = tk.Canvas(self, bg="white", highlightthickness=0, cursor="sb_h_double_arrow", height=175)
        self.canvas.pack(fill="both", expand=True)

        self.bins_edges = np.linspace(1.5, 6.4, 85)
        self.hist_series = []
        self.excl_log = 3.00
        self.cutoff_log = 4.85
        self.y_mode = "pct"
        self.v2l_summary_text = ""
        self.scope_text = "Global"
        self.dragging_which = None
        self.box = (60, 24, 500, 135)

        self.canvas.bind("<Configure>", lambda e: self.redraw())
        self.canvas.bind("<ButtonPress-1>", self._on_press)
        self.canvas.bind("<B1-Motion>", self._on_motion)
        self.canvas.bind("<ButtonRelease-1>", self._on_release)

    def set_histogram_data(self, bins_edges, hist_series, cutoff_log, excl_log=3.00, y_mode="pct", v2l_summary_text="", scope_text="Global"):
        self.bins_edges = bins_edges
        self.hist_series = hist_series
        self.cutoff_log = float(cutoff_log)
        self.excl_log = float(excl_log)
        self.y_mode = y_mode
        self.v2l_summary_text = v2l_summary_text
        self.scope_text = scope_text
        self.redraw()

    def _x_to_px(self, x_val):
        x0, _, x1, _ = self.box
        xmin, xmax = float(self.bins_edges[0]), float(self.bins_edges[-1])
        return x0 + ((x_val - xmin) / max(1e-9, xmax - xmin)) * max(1.0, x1 - x0)

    def _px_to_x(self, px):
        x0, _, x1, _ = self.box
        xmin, xmax = float(self.bins_edges[0]), float(self.bins_edges[-1])
        frac = (px - x0) / max(1.0, x1 - x0)
        return float(np.clip(xmin + frac * (xmax - xmin), xmin, xmax - 0.05))

    @staticmethod
    def _blend_hex(hex_col, alpha=0.20):
        h = hex_col.lstrip('#')
        if len(h) != 6:
            return "#D6EAF8"
        r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
        r_b = int(r * alpha + 255 * (1.0 - alpha))
        g_b = int(g * alpha + 255 * (1.0 - alpha))
        b_b = int(b * alpha + 255 * (1.0 - alpha))
        return f"#{r_b:02x}{g_b:02x}{b_b:02x}"

    def redraw(self):
        self.canvas.delete("all")
        W = max(480, self.canvas.winfo_width())
        H = max(150, self.canvas.winfo_height())
        ml, mr, mt, mb = 68, 245, 24, 30
        x0, y0, x1, y1 = ml, mt, max(ml + 160, W - mr), max(mt + 60, H - mb)
        self.box = (x0, y0, x1, y1)

        px_excl = np.clip(self._x_to_px(self.excl_log), x0, x1)
        px_cut = np.clip(self._x_to_px(self.cutoff_log), px_excl, x1)

        if px_excl > x0:
            self.canvas.create_rectangle(x0, y0, px_excl, y1, fill="#EAECEE", outline="")
        self.canvas.create_rectangle(px_excl, y0, px_cut, y1, fill="#FDEDEC", outline="")
        self.canvas.create_rectangle(px_cut, y0, x1, y1, fill="#EAFAF1", outline="")

        max_y = max([float(np.max(s['y'])) for s in self.hist_series if len(s['y']) > 0] + [5.0]) * 1.16
        centers = 0.5 * (self.bins_edges[:-1] + self.bins_edges[1:])

        for s in self.hist_series:
            col = s['color']
            fill_col = self._blend_hex(col, 0.22)
            y_vals = s['y']
            if len(y_vals) != len(centers):
                continue
            poly_pts = [self._x_to_px(centers[0]), y1]
            for cx, cy in zip(centers, y_vals):
                px = self._x_to_px(cx)
                py = y1 - np.clip(cy / max_y, 0.0, 1.0) * (y1 - y0)
                poly_pts.extend([px, py])
            poly_pts.extend([self._x_to_px(centers[-1]), y1])
            if len(poly_pts) >= 6:
                self.canvas.create_polygon(poly_pts, fill=fill_col, outline="")

        xmin, xmax = float(self.bins_edges[0]), float(self.bins_edges[-1])
        for gx in np.linspace(xmin, xmax, 8):
            px = self._x_to_px(gx)
            self.canvas.create_line(px, y0, px, y1, fill="#D5D8DC", dash=(2, 3))
            self.canvas.create_text(px, y1 + 9, text=f"{gx:.1f}", font=("DejaVu Sans", 8), fill="#333333")

        for gy in np.linspace(0.0, max_y, 4):
            py = y1 - (gy / max_y) * (y1 - y0)
            self.canvas.create_line(x0, py, x1, py, fill="#E5E8E8", dash=(2, 3))
            lbl = f"{gy:.0f}%" if self.y_mode == "pct" else f"{gy:.0f}"
            self.canvas.create_text(x0 - 5, py, text=lbl, anchor="e", font=("DejaVu Sans", 8), fill="#333333")

        for s in self.hist_series:
            col = s['color']
            y_vals = s['y']
            if len(y_vals) != len(centers):
                continue
            line_pts = []
            for cx, cy in zip(centers, y_vals):
                px = self._x_to_px(cx)
                py = y1 - np.clip(cy / max_y, 0.0, 1.0) * (y1 - y0)
                line_pts.extend([px, py])
            dash_pat = (4, 2) if s.get('dashed', False) else None
            if len(line_pts) >= 4:
                self.canvas.create_line(line_pts, fill=col, width=2, dash=dash_pat)

        self.canvas.create_rectangle(x0, y0, x1, y1, fill="", outline="#333333", width=1)

        self.canvas.create_text(
            (x0 + x1) // 2, y0 - 11,
            text=f"Histograma Dinámico (Muestras Seleccionadas) — Exclusión V1L < {self.excl_log:.2f} | Cut-Off V2L/V2R = {self.cutoff_log:.2f} [{self.scope_text}]",
            font=("DejaVu Sans", 9, "bold"), fill="#002B49"
        )
        if px_excl - x0 > 55:
            self.canvas.create_text(x0 + 5, y0 + 8, text="V1L: Excluidos", anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#4A235A")
        self.canvas.create_text(max(x0 + 6, px_excl + 6), y0 + 8, text="◄ V2L: Afectados", anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#B03A2E")
        self.canvas.create_text(x1 - 6, y0 + 8, text="V2R: Sanos (Intactos) ►", anchor="ne", font=("DejaVu Sans", 8, "bold"), fill="#1E8449")

        ylabel = "Eventos (%)" if self.y_mode == "pct" else "N° Eventos"
        self.canvas.create_text(x0 - 40, (y0 + y1) // 2, text=ylabel, font=("DejaVu Sans", 8, "bold"), fill="#1C2833")
        self.canvas.create_text((x0 + x1) // 2, y1 + 20, text="log10(Fluorescencia Canal Seleccionado)", font=("DejaVu Sans", 8, "bold"), fill="#002B49")

        self.canvas.create_line(px_excl, y0, px_excl, y1, fill="#4A235A", width=2, dash=(3, 2))
        self.canvas.create_rectangle(px_excl - 5, y0 + 22, px_excl + 5, y0 + 36, fill="#4A235A", outline="white")
        self.canvas.create_text(px_excl + 7, y0 + 29, text=f"Excl V1L ({self.excl_log:.2f})", anchor="w", font=("DejaVu Sans", 8, "bold"), fill="#4A235A")

        self.canvas.create_line(px_cut, y0, px_cut, y1, fill="#C0392B", width=3, dash=(5, 2))
        self.canvas.create_rectangle(px_cut - 6, y0 + 4, px_cut + 6, y0 + 18, fill="#C0392B", outline="white")
        self.canvas.create_text(px_cut + 8, y0 + 11, text=f"Cut-Off V2L ({self.cutoff_log:.2f})", anchor="w", font=("DejaVu Sans", 8, "bold"), fill="#C0392B")

        lx = x1 + 10
        ly = y0 + 2
        self.canvas.create_text(lx, ly, text=f"Muestras en Histograma ({len(self.hist_series)}):", anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#002B49")
        ly += 15
        for s in self.hist_series[:6]:
            dash_p = (3, 2) if s.get('dashed', False) else None
            self.canvas.create_line(lx, ly + 6, lx + 16, ly + 6, fill=s['color'], width=3, dash=dash_p)
            self.canvas.create_text(lx + 20, ly, text=s['label'][:28], anchor="nw", font=("DejaVu Sans", 8), fill="#1C2833")
            ly += 14
        if len(self.hist_series) > 6:
            self.canvas.create_text(lx, ly, text=f"+ {len(self.hist_series) - 6} curvas adicionales...", anchor="nw", font=("DejaVu Sans", 8), fill="#566573")
            ly += 14
        if self.v2l_summary_text:
            ly += 2
            self.canvas.create_text(lx, ly, text=self.v2l_summary_text, anchor="nw", font=("DejaVu Sans", 8, "bold"), fill="#922B21")

    def _on_press(self, event):
        x0, y0, x1, y1 = self.box
        if x0 <= event.x <= x1 and y0 <= event.y <= y1:
            px_excl = self._x_to_px(self.excl_log)
            px_cut = self._x_to_px(self.cutoff_log)
            if abs(event.x - px_excl) <= 14 and abs(event.x - px_excl) < abs(event.x - px_cut):
                self.dragging_which = 'excl'
            else:
                self.dragging_which = 'v2l'
            self._update_from_mouse(event.x, final=False)

    def _on_motion(self, event):
        if self.dragging_which:
            self._update_from_mouse(event.x, final=False)

    def _on_release(self, event):
        if self.dragging_which:
            self._update_from_mouse(event.x, final=True)
            self.dragging_which = None

    def _update_from_mouse(self, px, final=False):
        val = round(self._px_to_x(px), 2)
        if self.dragging_which == 'excl':
            self.excl_log = min(val, self.cutoff_log - 0.15)
        else:
            self.cutoff_log = max(val, self.excl_log + 0.15)
        self.redraw()
        if self.on_cutoff_change_callback:
            self.on_cutoff_change_callback(self.cutoff_log, self.excl_log, final=final)
