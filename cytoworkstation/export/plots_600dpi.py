"""
Motor de Renderizado Gráfico a 600 DPI (TIFF / PNG) para Publicaciones Científicas.
Genera figuras con antialiasing, bandas de error SEM, marcadores diferenciados y tipografía nítida.
"""
from pathlib import Path
from datetime import datetime
import struct
import zlib
import numpy as np
import pandas as pd

from ..config import HAS_MPL

if HAS_MPL:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg


FONT_5X7 = {
    ' ': [0, 0, 0, 0, 0, 0, 0],
    '0': [14, 17, 19, 21, 25, 17, 14], '1': [4, 12, 4, 4, 4, 4, 14],
    '2': [14, 17, 1, 6, 8, 16, 31], '3': [31, 2, 4, 2, 1, 17, 14],
    '4': [2, 6, 10, 18, 31, 2, 2], '5': [31, 16, 30, 1, 1, 17, 14],
    '6': [6, 8, 16, 30, 17, 17, 14], '7': [31, 1, 2, 4, 8, 8, 8],
    '8': [14, 17, 17, 14, 17, 17, 14], '9': [14, 17, 17, 15, 1, 2, 12],
    'A': [14, 17, 17, 31, 17, 17, 17], 'B': [30, 17, 17, 30, 17, 17, 30],
    'C': [14, 17, 16, 16, 16, 17, 14], 'D': [28, 18, 17, 17, 17, 18, 28],
    'E': [31, 16, 16, 30, 16, 16, 31], 'F': [31, 16, 16, 30, 16, 16, 16],
    'G': [14, 17, 16, 23, 17, 17, 15], 'H': [17, 17, 17, 31, 17, 17, 17],
    'I': [14, 4, 4, 4, 4, 4, 14], 'J': [7, 2, 2, 2, 2, 18, 12],
    'K': [17, 18, 20, 24, 20, 18, 17], 'L': [16, 16, 16, 16, 16, 16, 31],
    'M': [17, 27, 21, 21, 17, 17, 17], 'N': [17, 25, 21, 19, 17, 17, 17],
    'O': [14, 17, 17, 17, 17, 17, 14], 'P': [30, 17, 17, 30, 16, 16, 16],
    'Q': [14, 17, 17, 17, 21, 18, 13], 'R': [30, 17, 17, 30, 20, 18, 17],
    'S': [15, 16, 16, 14, 1, 1, 30], 'T': [31, 4, 4, 4, 4, 4, 4],
    'U': [17, 17, 17, 17, 17, 17, 14], 'V': [17, 17, 17, 17, 17, 10, 4],
    'W': [17, 17, 17, 21, 21, 21, 10], 'X': [17, 17, 10, 4, 10, 17, 17],
    'Y': [17, 17, 10, 4, 4, 4, 4], 'Z': [31, 1, 2, 4, 8, 16, 31],
    '.': [0, 0, 0, 0, 0, 12, 12], '-': [0, 0, 0, 31, 0, 0, 0],
    '_': [0, 0, 0, 0, 0, 0, 31], '/': [1, 2, 2, 4, 8, 8, 16],
    '%': [24, 25, 2, 4, 8, 19, 3], '(': [2, 4, 8, 8, 8, 4, 2],
    ')': [8, 4, 2, 2, 2, 4, 8], ':': [0, 12, 12, 0, 12, 12, 0],
    '|': [4, 4, 4, 4, 4, 4, 4], '+': [0, 4, 4, 31, 4, 4, 0],
    '=': [0, 0, 31, 0, 31, 0, 0], '<': [2, 4, 8, 16, 8, 4, 2],
    '>': [8, 4, 2, 1, 2, 4, 8], '[': [14, 8, 8, 8, 8, 8, 14],
    ']': [14, 2, 2, 2, 2, 2, 14], ',': [0, 0, 0, 0, 6, 4, 8]
}


def draw_bitmap_text(buf, x, y, text, rgb=(20, 20, 20), scale=3):
    """Dibuja texto nítido directamente en un buffer numpy RGB (H, W, 3)."""
    H, W, _ = buf.shape
    cx = int(x)
    cy = int(y)
    clean_txt = str(text).replace('µ', 'U').upper()
    for ch in clean_txt:
        glyph = FONT_5X7.get(ch, FONT_5X7[' '])
        for row_idx, row_bits in enumerate(glyph):
            for col_idx in range(5):
                if (row_bits >> (4 - col_idx)) & 1:
                    y0 = cy + row_idx * scale
                    x0 = cx + col_idx * scale
                    if 0 <= y0 < H - scale and 0 <= x0 < W - scale:
                        buf[y0:y0 + scale, x0:x0 + scale] = rgb
        cx += 6 * scale


def hex_to_rgb(hex_str):
    h = hex_str.lstrip('#')
    if len(h) == 6:
        return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return 0, 68, 136


def save_rgb_buffer_as_png(filepath, buf, width, height, dpi=600):
    """Codifica una matriz numpy RGB (H, W, 3) uint8 en un archivo PNG certificado a 600 DPI."""
    ppm = int(round(dpi / 0.0254))
    raw_rows = b"".join(b"\x00" + buf[y].tobytes() for y in range(height))
    compressed = zlib.compress(raw_rows, level=6)

    def make_chunk(chunk_type, data):
        crc = zlib.crc32(chunk_type + data) & 0xffffffff
        return struct.pack(">I", len(data)) + chunk_type + data + struct.pack(">I", crc)

    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    phys = struct.pack(">IIB", ppm, ppm, 1)
    png_bytes = (
        b"\x89PNG\r\n\x1a\n" +
        make_chunk(b"IHDR", ihdr) +
        make_chunk(b"pHYs", phys) +
        make_chunk(b"IDAT", compressed) +
        make_chunk(b"IEND", b"")
    )
    target_path = Path(filepath)
    try:
        with open(target_path, "wb") as f:
            f.write(png_bytes)
    except PermissionError:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        target_path = target_path.with_name(f"{target_path.stem}_{ts}{target_path.suffix}")
        with open(target_path, "wb") as f:
            f.write(png_bytes)
    return str(target_path)


def save_rgb_buffer_as_tiff(filepath, buf, width, height, dpi=600):
    """Codifica una matriz numpy RGB (H, W, 3) uint8 en TIFF 6.0 Baseline RGB a 600 DPI."""
    target_path = Path(filepath).with_suffix(".tiff")
    raw_rgb = np.ascontiguousarray(buf[:height, :width, :3], dtype=np.uint8).tobytes()

    num_tags = 12
    ifd_offset = 8
    after_ifd = ifd_offset + 2 + (num_tags * 12) + 4
    bits_per_sample_offset = after_ifd
    xres_offset = bits_per_sample_offset + 6
    yres_offset = xres_offset + 8
    pixel_data_offset = yres_offset + 8

    def ifd_entry(tag, ttype, count, val_or_offset):
        return struct.pack("<HHII", tag, ttype, count, val_or_offset)

    entries = [
        ifd_entry(256, 4, 1, int(width)),
        ifd_entry(257, 4, 1, int(height)),
        ifd_entry(258, 3, 3, bits_per_sample_offset),
        ifd_entry(259, 3, 1, 1),
        ifd_entry(262, 3, 1, 2),
        ifd_entry(273, 4, 1, pixel_data_offset),
        ifd_entry(277, 3, 1, 3),
        ifd_entry(278, 4, 1, int(height)),
        ifd_entry(279, 4, 1, len(raw_rgb)),
        ifd_entry(282, 5, 1, xres_offset),
        ifd_entry(283, 5, 1, yres_offset),
        ifd_entry(296, 3, 1, 2),
    ]

    header = struct.pack("<2sHI", b"II", 42, ifd_offset)
    ifd_block = struct.pack("<H", num_tags) + b"".join(entries) + struct.pack("<I", 0)
    extra_vals = (
        struct.pack("<HHH", 8, 8, 8) +
        struct.pack("<II", int(dpi), 1) +
        struct.pack("<II", int(dpi), 1)
    )
    tiff_bytes = header + ifd_block + extra_vals + raw_rgb
    try:
        with open(target_path, "wb") as f:
            f.write(tiff_bytes)
    except PermissionError:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        target_path = target_path.with_name(f"{target_path.stem}_{ts}.tiff")
        with open(target_path, "wb") as f:
            f.write(tiff_bytes)
    return str(target_path)


def save_mpl_figure_600dpi_tiff(fig, filepath):
    """Guarda una figura Matplotlib en TIFF a 600 DPI con compresión LZW o estándar."""
    target_path = Path(filepath).with_suffix(".tiff")
    save_kwargs = dict(format="tiff", dpi=600, bbox_inches="tight", facecolor="white")
    try:
        try:
            fig.savefig(target_path, pil_kwargs={"compression": "tiff_lzw"}, **save_kwargs)
        except TypeError:
            fig.savefig(target_path, **save_kwargs)
    except PermissionError:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        target_path = target_path.with_name(f"{target_path.stem}_{ts}.tiff")
        try:
            fig.savefig(target_path, pil_kwargs={"compression": "tiff_lzw"}, **save_kwargs)
        except TypeError:
            fig.savefig(target_path, **save_kwargs)
    return target_path


def render_highres_plot_png(filepath, series_list, vlines, title, xlabel, ylabel, xlim, ylim):
    """Renderiza una figura cinética en TIFF a 600 DPI y devuelve el buffer RGB para el PDF."""
    target_path = Path(filepath).with_suffix(".tiff")
    if HAS_MPL:
        try:
            fig = Figure(figsize=(8.6, 5.0), dpi=600, facecolor='white')
            canvas = FigureCanvasAgg(fig)
            ax = fig.add_subplot(111)

            for idx_s, s in enumerate(series_list):
                xs_arr = np.asarray(s['x'], dtype=float)
                ys_arr = np.asarray(s['y'], dtype=float)
                se_arr = np.asarray(s.get('se', np.zeros_like(ys_arr)), dtype=float)
                col = s.get('color', '#004488')
                is_dashed = bool(s.get('dashed', False))
                ls = ':' if is_dashed else '--'
                mk = 'o' if is_dashed else ('x' if 'THP' in str(s.get('label', '')) else 's')
                lbl = str(s.get('label', f'Serie {idx_s+1}'))

                ax.plot(
                    xs_arr, ys_arr, label=lbl, color=col, linestyle=ls,
                    marker=mk, markersize=4.5, linewidth=1.85, alpha=0.95
                )
                if len(se_arr) == len(ys_arr) and np.any(se_arr > 0):
                    ax.fill_between(xs_arr, ys_arr - se_arr, ys_arr + se_arr, color=col, alpha=0.18, linewidth=0)

            for v_item in vlines:
                vx, vcol, vdash = v_item[0], v_item[1], v_item[2]
                vlbl = v_item[3] if len(v_item) > 3 else None
                if xlim[0] <= vx <= xlim[1]:
                    v_ls = '--' if vx == 0.0 else (':' if vdash else '-.')
                    v_lw = 1.2 if vx == 0.0 else 1.45
                    ax.axvline(vx, color=vcol, linestyle=v_ls, linewidth=v_lw, alpha=0.9, label=vlbl)

            if ylim[0] < 0 <= ylim[1]:
                ax.axhline(0, color='gray', linestyle=':', linewidth=0.9, alpha=0.6)

            ax.set_title(title, weight='bold', fontsize=10.5, color='#002B49', pad=10)
            ax.set_xlabel(xlabel, weight='bold', fontsize=10, color='#1C2833')
            ax.set_ylabel(ylabel, weight='bold', fontsize=10, color='#1C2833')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.grid(True, linestyle=':', alpha=0.45, color='#BFC9CA')
            ax.spines['top'].set_visible(True)
            ax.spines['right'].set_visible(True)

            if series_list:
                ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, facecolor='white', edgecolor='#D5D8DC', fontsize=7.5)
            fig.tight_layout()

            save_mpl_figure_600dpi_tiff(fig, target_path)

            fig.set_dpi(220)
            canvas.draw()
            rgba = np.asarray(canvas.buffer_rgba())
            return rgba[..., :3].copy()
        except Exception as e:
            print(f"[Aviso Matplotlib Plot]: {e}")

    # Fallback nativo
    W, H = 2400, 1450
    buf = np.full((H, W, 3), 255, dtype=np.uint8)
    ml, mr, mt, mb = 250, 660, 140, 190
    pw, ph = W - ml - mr, H - mt - mb

    def to_px(x, y):
        px = int(ml + (x - xlim[0]) / max(1e-9, xlim[1] - xlim[0]) * pw)
        py = int(mt + ph - (y - ylim[0]) / max(1e-9, ylim[1] - ylim[0]) * ph)
        return np.clip(px, ml, ml + pw - 1), np.clip(py, mt, mt + ph - 1)

    draw_bitmap_text(buf, ml, 55, title[:62], rgb=(0, 43, 73), scale=5)
    save_rgb_buffer_as_tiff(filepath, buf, W, H, dpi=600)
    return buf[::3, ::3].copy()


def render_highres_histogram_png(filepath, bins_edges, hist_series, cutoff_log, title, ylabel, excl_log=3.00):
    """Genera el histograma de fluorescencia con zonas sombreadas V1L, V2L y V2R a 600 DPI."""
    target_path = Path(filepath).with_suffix(".tiff")
    if HAS_MPL:
        try:
            fig = Figure(figsize=(8.8, 4.8), dpi=600, facecolor='white')
            canvas = FigureCanvasAgg(fig)
            ax = fig.add_subplot(111)

            xlim = (float(bins_edges[0]), float(bins_edges[-1]))
            max_y = max([float(np.max(s['y'])) for s in hist_series if len(s['y']) > 0] + [10.0]) * 1.18
            centers = 0.5 * (bins_edges[:-1] + bins_edges[1:])

            if excl_log > xlim[0]:
                ax.axvspan(xlim[0], excl_log, color='#8E44AD', alpha=0.10, label='Zona Excluida V1L')
            ax.axvspan(excl_log, cutoff_log, color='#E74C3C', alpha=0.12, label='Zona Afectados V2L')
            ax.axvspan(cutoff_log, xlim[1], color='#27AE60', alpha=0.10, label='Zona Intactas V2R')

            for s in hist_series:
                col = s.get('color', '#004488')
                lbl = s.get('label', '')
                ax.plot(centers, s['y'], color=col, linewidth=1.8, label=lbl)

            if excl_log > xlim[0]:
                ax.axvline(excl_log, color='#4A235A', linestyle=':', linewidth=2.0, label=f'Cut-Off 1 V1L ({excl_log:.2f})')
            ax.axvline(cutoff_log, color='#C0392B', linestyle='--', linewidth=2.2, label=f'Cut-Off 2 V2L ({cutoff_log:.2f})')

            ax.set_title(title, weight='bold', fontsize=10.0, color='#002B49', pad=10)
            ax.set_xlabel("log10(Fluorescencia Canal FITC)", weight='bold', fontsize=10)
            ax.set_ylabel(ylabel, weight='bold', fontsize=10)
            ax.set_xlim(xlim)
            ax.set_ylim(0, max_y)
            ax.grid(True, linestyle=':', alpha=0.45, color='#BFC9CA')
            ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, facecolor='white', edgecolor='#D5D8DC', fontsize=7.2)
            fig.tight_layout()

            save_mpl_figure_600dpi_tiff(fig, target_path)
            fig.set_dpi(220)
            canvas.draw()
            rgba = np.asarray(canvas.buffer_rgba())
            return rgba[..., :3].copy()
        except Exception as e:
            print(f"[Aviso Matplotlib Histogram]: {e}")

    # Fallback
    W, H = 2400, 1350
    buf = np.full((H, W, 3), 255, dtype=np.uint8)
    save_rgb_buffer_as_tiff(filepath, buf, W, H, dpi=600)
    return buf[::3, ::3].copy()


def render_highres_categorical_milestones_png(filepath, df_cat_summary, target_cell="THP-1"):
    """
    Renderiza en formato TIFF a 600 DPI un panel doble para el Reporte Final integrando
    los Cut-Offs específicos configurados por el usuario para cada muestra:
      - Subpanel A: Hitos Temporales del Efecto Categórico (T_onset estadístico, T50 del efecto y T_efecto_máx) por Compuesto y Concentración.
      - Subpanel B: Proporciones Complementarias en el Efecto Máximo (% Afectados V2L vs % No Afectados V2R) y % Mínimo.
    """
    target_path = Path(filepath).with_suffix(".tiff")
    if df_cat_summary is None or df_cat_summary.empty:
        return None

    sub_df = df_cat_summary[df_cat_summary['Cell_Line'] == target_cell].copy()
    if sub_df.empty:
        sub_df = df_cat_summary.copy()
    sub_df = sub_df.sort_values(['Peptide', 'Concentration_uM']).reset_index(drop=True)

    labels = []
    for _, r in sub_df.iterrows():
        cut_str = f"\n[V2L<{float(r['V2L_Cutoff_Log10']):.2f}]" if 'V2L_Cutoff_Log10' in r and pd.notna(r['V2L_Cutoff_Log10']) else ""
        labels.append(f"{r['Peptide']}\n{r['Concentration_uM']:.2f} µM{cut_str}")

    t_onset = sub_df['T_Stat_Onset_s'].values.astype(float)
    t50_eff = sub_df['T50_Effect_s'].values.astype(float)
    t_max_eff = sub_df['T_Effect_Max_s'].values.astype(float)
    max_v2l = sub_df['Max_V2L_Pct'].values.astype(float)
    min_v2l = sub_df['Min_V2L_Pct'].values.astype(float)
    min_v2r = sub_df['Min_Unaffected_V2R_Pct'].values.astype(float)

    if HAS_MPL and len(labels) > 0:
        try:
            fig = Figure(figsize=(9.5, 4.9), dpi=600, facecolor='white')
            canvas = FigureCanvasAgg(fig)
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)

            x_idx = np.arange(len(labels))
            bw = 0.25

            ax1.bar(x_idx - bw, t_onset, width=bw, color='#2874A6', label='T_onset estadístico (s)', alpha=0.9)
            ax1.bar(x_idx, t50_eff, width=bw, color='#E67E22', label='T50 del Efecto (s)', alpha=0.9)
            ax1.bar(x_idx + bw, t_max_eff, width=bw, color='#922B21', label='T_efecto_máx (s)', alpha=0.9)
            ax1.set_xticks(x_idx)
            ax1.set_xticklabels(labels, fontsize=7.0, rotation=25, ha='right')
            ax1.set_ylabel("Tiempo (segundos)", weight='bold', fontsize=9.5)
            ax1.set_ylim(0, 650)
            ax1.set_title(f"A. Hitos Temporales del Efecto ({target_cell})", weight='bold', fontsize=10, color='#002B49')
            ax1.grid(True, axis='y', linestyle=':', alpha=0.5)
            ax1.legend(loc='upper right', fontsize=7.5, frameon=True)

            ax2.bar(x_idx, max_v2l, width=0.52, color='#C0392B', label='% Afectados Máx (V2L)', alpha=0.88)
            ax2.bar(x_idx, min_v2r, width=0.52, bottom=max_v2l, color='#27AE60', label='% No Afectados Remanente (V2R)', alpha=0.78)
            ax2.plot(x_idx, min_v2l, 'o--', color='#1B4F72', markersize=4.5, linewidth=1.3, label='% Afectados Mín (V2L)')

            for i, (mv, st) in enumerate(zip(max_v2l, sub_df['Affectation_State'].values)):
                sig_tok = "***" if "***" in str(st) else ("**" if "**" in str(st) else ("*" if "*" in str(st) else "NS"))
                ax2.text(i, min(96.0, mv + 2.5), f"{mv:.1f}%\n({sig_tok})", ha='center', va='bottom', fontsize=6.8, weight='bold', color='#1C2833')

            ax2.set_xticks(x_idx)
            ax2.set_xticklabels(labels, fontsize=7.0, rotation=25, ha='right')
            ax2.set_ylabel("Proporción de Eventos (%)", weight='bold', fontsize=9.5)
            ax2.set_ylim(0, 112)
            ax2.set_title(f"B. Balance Afectados (V2L) vs No Afectados (V2R)", weight='bold', fontsize=10, color='#002B49')
            ax2.grid(True, axis='y', linestyle=':', alpha=0.5)
            ax2.legend(loc='upper right', fontsize=7.2, frameon=True)

            fig.tight_layout()
            save_mpl_figure_600dpi_tiff(fig, target_path)

            fig.set_dpi(220)
            canvas.draw()
            rgba = np.asarray(canvas.buffer_rgba())
            return rgba[..., :3].copy()
        except Exception as e:
            print(f"[Aviso Matplotlib Milestones]: {e}")

    # Fallback nativo
    series_fallback = [
        {'x': np.arange(len(labels)), 'y': t50_eff, 'se': np.zeros_like(t50_eff), 'color': '#E67E22', 'dashed': False, 'label': 'T50 Efecto (s)'},
        {'x': np.arange(len(labels)), 'y': t_max_eff, 'se': np.zeros_like(t_max_eff), 'color': '#922B21', 'dashed': False, 'label': 'T_efecto_max (s)'}
    ]
    return render_highres_plot_png(
        target_path, series_fallback, [],
        f"Hitos Temporales T50 y T_max ({target_cell})", "CONDICION (COMPUESTO / DOSIS)", "TIEMPO (S)",
        (0.0, max(1.0, float(len(labels) - 1))), (0.0, 620.0)
    )


def render_highres_categorical_heatmap_png(filepath, df_cat_windows, target_cell="THP-1"):
    """
    Renderiza en formato TIFF a 600 DPI un Mapa de Calor (Heatmap) de la Proporción de Eventos Afectados (% V2L)
    por Compuesto, Concentración (con su Cut-Off específico V2L configurado por el usuario) y Ventana de Tiempo (-60 a 600 s).
    """
    target_path = Path(filepath).with_suffix(".tiff")
    if df_cat_windows is None or df_cat_windows.empty:
        return None

    sub = df_cat_windows[df_cat_windows['Cell_Line'] == target_cell].copy()
    if sub.empty:
        sub = df_cat_windows.copy()

    def _fmt_row_id(r):
        cut_part = f" [V2L<{float(r['V2L_Cutoff_Log10']):.2f}]" if 'V2L_Cutoff_Log10' in r and pd.notna(r['V2L_Cutoff_Log10']) else ""
        return f"{r['Peptide']} | {r['Concentration_uM']:.2f} µM ({r['Cell_Line']}){cut_part}"

    sub['Row_ID'] = sub.apply(_fmt_row_id, axis=1)
    piv = sub.pivot_table(index='Row_ID', columns='Time_Window', values='Mean_Affected_V2L_Pct', aggfunc='mean').fillna(0.0)
    if piv.empty:
        return None

    if HAS_MPL:
        try:
            fig = Figure(figsize=(9.4, max(3.9, 0.44 * len(piv.index) + 1.8)), dpi=600, facecolor='white')
            canvas = FigureCanvasAgg(fig)
            ax = fig.add_subplot(111)

            mat = piv.values
            im = ax.imshow(mat, aspect='auto', cmap='YlOrRd', vmin=0.0, vmax=100.0)
            ax.set_yticks(np.arange(len(piv.index)))
            ax.set_yticklabels(list(piv.index), fontsize=7.8, weight='bold')

            cols_t = list(piv.columns)
            step_x = max(1, len(cols_t) // 12)
            xt_idx = np.arange(0, len(cols_t), step_x)
            ax.set_xticks(xt_idx)
            ax.set_xticklabels([f"{int(cols_t[i])}s" for i in xt_idx], fontsize=8)

            ax.set_title("Mapa de Calor Espaciotemporal: % de Eventos Afectados (V2L) por Ventana de Tiempo y Cut-Off Configurado", weight='bold', fontsize=10.2, color='#002B49', pad=10)
            ax.set_xlabel("Ventana de Tiempo (segundos)", weight='bold', fontsize=9.5)
            ax.set_ylabel("Compuesto | Concentración [Cut-Off V2L]", weight='bold', fontsize=9.5)

            cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
            cbar.set_label("% Eventos Afectados (V2L)", fontsize=8.5, weight='bold')

            fig.tight_layout()
            save_mpl_figure_600dpi_tiff(fig, target_path)

            fig.set_dpi(220)
            canvas.draw()
            rgba = np.asarray(canvas.buffer_rgba())
            return rgba[..., :3].copy()
        except Exception as e:
            print(f"[Aviso Matplotlib Heatmap]: {e}")

    return None


def render_highres_target_vs_control_overlap_png(filepath, df_full, target_cell="THP-1", control_cell="PBMCs", df_kinetics=None):
    """
    Exporta una figura comparativa de alta resolución en formato TIFF a 600 DPI con las curvas superpuestas
    de Células Diana (target_cell, trazo continuo y marcador sólido) y Células Control
    (control_cell, trazo discontinuo y marcador hueco) a través de los bins de tiempo, integrando
    los Cut-Offs (V2L, T_lag y T_max) configurados por el usuario para cada muestra:
      - Panel A: Cinética Continua de Permeabilización Superpuesta (% Mediana FITC ± SEM) y líneas T_lag/T_max.
      - Panel B: Cinética Categórica de Eventos Afectados Superpuesta (% V2L ± SEM) con su Cut-Off específico.
    Retorna además el buffer RGB (H, W, 3) uint8 para incrustarlo en el reporte PDF.
    """
    if df_full is None or df_full.empty:
        return None
    target_path = Path(filepath).with_suffix(".tiff")

    if HAS_MPL:
        try:
            fig = Figure(figsize=(12.6, 5.6), dpi=600, facecolor='white')
            canvas = FigureCanvasAgg(fig)
            ax1 = fig.add_subplot(1, 2, 1)
            ax2 = fig.add_subplot(1, 2, 2)

            for ax in (ax1, ax2):
                ax.set_facecolor('#FAFCFC')
                ax.grid(True, linestyle='--', linewidth=0.5, color='#D5DBDB', alpha=0.75)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.axvline(0.0, color='#333333', linestyle='--', linewidth=1.1, alpha=0.8)

            palette = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#8E44AD', '#C0392B']
            combos = sorted(
                df_full[['Peptide', 'Concentration_uM']].drop_duplicates().itertuples(index=False),
                key=lambda x: (str(x.Peptide), float(x.Concentration_uM))
            )

            any_plotted = False
            for idx_c, combo in enumerate(combos):
                pep = combo.Peptide
                conc = float(combo.Concentration_uM)
                col_hex = palette[idx_c % len(palette)]
                sub_c = df_full[(df_full['Peptide'] == pep) & (np.isclose(df_full['Concentration_uM'], conc))]

                # Célula Diana (Target Cell) -> línea sólida con marcador lleno + Cut-Offs configurados
                sub_t = sub_c[sub_c['Cell_Line'] == target_cell]
                agg_t = None
                if not sub_t.empty:
                    agg_t = sub_t.groupby('Time_Window').agg(
                        m_perm=('Permeabilization_Pct', 'mean'),
                        se_perm=('Permeabilization_Pct', 'sem'),
                        m_v2l=('Affected_V2L_Pct', 'mean'),
                        se_v2l=('Affected_V2L_Pct', 'sem')
                    ).reset_index().sort_values('Time_Window')
                    tx = agg_t['Time_Window'].values
                    ty_p = agg_t['m_perm'].values
                    tse_p = agg_t['se_perm'].fillna(0.0).values
                    ty_v = agg_t['m_v2l'].values
                    tse_v = agg_t['se_v2l'].fillna(0.0).values

                    cut_t_v2l = float(sub_t['Applied_V2L_Cutoff_Log10'].iloc[0]) if 'Applied_V2L_Cutoff_Log10' in sub_t else 4.85
                    tlag_t, tmax_t = None, None
                    if df_kinetics is not None and not df_kinetics.empty:
                        kt_sub = df_kinetics[(df_kinetics['Peptide'] == pep) & (df_kinetics['Cell_Line'] == target_cell) & (np.isclose(df_kinetics['Concentration_uM'], conc))]
                        if not kt_sub.empty:
                            tlag_t = float(kt_sub['T_lag'].mean())
                            tmax_t = float(kt_sub['T_max'].mean())

                    cut_info_t = f"Cut<{cut_t_v2l:.2f}"
                    if tlag_t is not None and tmax_t is not None:
                        cut_info_t += f", Tlag={tlag_t:.0f}s, Tmax={tmax_t:.0f}s"
                        if 0.0 < tlag_t < 600.0:
                            ax1.axvline(tlag_t, color=col_hex, linestyle=':', linewidth=1.1, alpha=0.65)
                            ax2.axvline(tlag_t, color=col_hex, linestyle=':', linewidth=1.1, alpha=0.65)
                        if 0.0 < tmax_t < 600.0:
                            ax1.axvline(tmax_t, color=col_hex, linestyle='-.', linewidth=1.1, alpha=0.65)
                            ax2.axvline(tmax_t, color=col_hex, linestyle='-.', linewidth=1.1, alpha=0.65)

                    lbl_t = f"{pep} {conc:.2f}µM [{target_cell} | {cut_info_t}]"
                    ax1.plot(tx, ty_p, linestyle='-', marker='o', markersize=3.2, linewidth=1.9, color=col_hex, label=lbl_t)
                    ax1.fill_between(tx, ty_p - tse_p, ty_p + tse_p, color=col_hex, alpha=0.14)
                    ax2.plot(tx, ty_v, linestyle='-', marker='o', markersize=3.2, linewidth=1.9, color=col_hex, label=lbl_t)
                    ax2.fill_between(tx, ty_v - tse_v, ty_v + tse_v, color=col_hex, alpha=0.14)
                    any_plotted = True

                # Célula Control (Control Cell) -> línea discontinua con marcador hueco + Cut-Offs configurados
                sub_ctrl = sub_c[sub_c['Cell_Line'] == control_cell]
                agg_ctrl = None
                if not sub_ctrl.empty:
                    agg_ctrl = sub_ctrl.groupby('Time_Window').agg(
                        m_perm=('Permeabilization_Pct', 'mean'),
                        se_perm=('Permeabilization_Pct', 'sem'),
                        m_v2l=('Affected_V2L_Pct', 'mean'),
                        se_v2l=('Affected_V2L_Pct', 'sem')
                    ).reset_index().sort_values('Time_Window')
                    cx = agg_ctrl['Time_Window'].values
                    cy_p = agg_ctrl['m_perm'].values
                    cse_p = agg_ctrl['se_perm'].fillna(0.0).values
                    cy_v = agg_ctrl['m_v2l'].values
                    cse_v = agg_ctrl['se_v2l'].fillna(0.0).values

                    cut_c_v2l = float(sub_ctrl['Applied_V2L_Cutoff_Log10'].iloc[0]) if 'Applied_V2L_Cutoff_Log10' in sub_ctrl else 4.85
                    lbl_ctrl = f"{pep} {conc:.2f}µM [{control_cell} | Cut<{cut_c_v2l:.2f}]"
                    ax1.plot(
                        cx, cy_p, linestyle='--', marker='s', markersize=3.0,
                        markerfacecolor='white', markeredgecolor=col_hex,
                        linewidth=1.5, color=col_hex, alpha=0.88, label=lbl_ctrl
                    )
                    ax1.fill_between(cx, cy_p - cse_p, cy_p + cse_p, color=col_hex, alpha=0.07)
                    ax2.plot(
                        cx, cy_v, linestyle='--', marker='s', markersize=3.0,
                        markerfacecolor='white', markeredgecolor=col_hex,
                        linewidth=1.5, color=col_hex, alpha=0.88, label=lbl_ctrl
                    )
                    ax2.fill_between(cx, cy_v - cse_v, cy_v + cse_v, color=col_hex, alpha=0.07)
                    any_plotted = True

                # Si existen ambas curvas en este compuesto/dosis, sombrear sutilmente la brecha de selectividad
                if agg_t is not None and agg_ctrl is not None:
                    merged_tc = pd.merge(
                        agg_t[['Time_Window', 'm_perm', 'm_v2l']],
                        agg_ctrl[['Time_Window', 'm_perm', 'm_v2l']],
                        on='Time_Window', suffixes=('_t', '_c')
                    )
                    if not merged_tc.empty:
                        ax1.fill_between(
                            merged_tc['Time_Window'].values,
                            merged_tc['m_perm_c'].values,
                            merged_tc['m_perm_t'].values,
                            color=col_hex, alpha=0.05
                        )
                        ax2.fill_between(
                            merged_tc['Time_Window'].values,
                            merged_tc['m_v2l_c'].values,
                            merged_tc['m_v2l_t'].values,
                            color=col_hex, alpha=0.05
                        )

            # Si las etiquetas de célula no coincidieron con target_cell/control_cell exactos, graficar todas las células presentes
            if not any_plotted:
                for idx_cell, ((pep, cell, conc), grp) in enumerate(df_full.groupby(['Peptide', 'Cell_Line', 'Concentration_uM'])):
                    col_hex = palette[idx_cell % len(palette)]
                    ls = '--' if cell == control_cell else '-'
                    agg_any = grp.groupby('Time_Window').agg(
                        m_perm=('Permeabilization_Pct', 'mean'),
                        m_v2l=('Affected_V2L_Pct', 'mean')
                    ).reset_index().sort_values('Time_Window')
                    ax1.plot(agg_any['Time_Window'].values, agg_any['m_perm'].values, linestyle=ls, marker='o', markersize=3.0, linewidth=1.6, color=col_hex, label=f"{pep} {conc:.2f}µM ({cell})")
                    ax2.plot(agg_any['Time_Window'].values, agg_any['m_v2l'].values, linestyle=ls, marker='o', markersize=3.0, linewidth=1.6, color=col_hex, label=f"{pep} {conc:.2f}µM ({cell})")

            ax1.set_title(f"A. Superposición Diana ({target_cell}) vs Control ({control_cell}) — Permeabilización Continua", fontsize=9.5, weight='bold', color='#002B49')
            ax1.set_xlabel("Bin de Tiempo (segundos)", fontsize=8.5, weight='bold')
            ax1.set_ylabel("Permeabilización Continua (%)", fontsize=8.5, weight='bold')
            ax1.set_xlim(-60.0, 600.0)
            ax1.set_ylim(-5.0, 105.0)
            ax1.legend(loc='best', fontsize=6.2, frameon=True, facecolor='white', framealpha=0.88, ncol=1)

            ax2.set_title(f"B. Superposición Diana ({target_cell}) vs Control ({control_cell}) — % Afectados (V2L) por Bin", fontsize=9.5, weight='bold', color='#002B49')
            ax2.set_xlabel("Bin de Tiempo (segundos)", fontsize=8.5, weight='bold')
            ax2.set_ylabel("Eventos Afectados V2L (%)", fontsize=8.5, weight='bold')
            ax2.set_xlim(-60.0, 600.0)
            ax2.set_ylim(-5.0, 105.0)
            ax2.legend(loc='best', fontsize=6.2, frameon=True, facecolor='white', framealpha=0.88, ncol=1)

            fig.suptitle(
                f"Análisis Comparativo en Gráfica Superpuesta por Bin de Tiempo (con Cut-Offs Configurados): Diana ({target_cell}) vs. Control ({control_cell})",
                fontsize=10.2, weight='bold', color='#002B49', y=0.99
            )
            fig.tight_layout(rect=[0, 0, 1, 0.95])
            save_mpl_figure_600dpi_tiff(fig, target_path)

            fig.set_dpi(220)
            canvas.draw()
            rgba = np.asarray(canvas.buffer_rgba())
            return rgba[..., :3].copy()
        except Exception as e:
            print(f"[Aviso Matplotlib Target vs Control Overlap]: {e}")

    return None
