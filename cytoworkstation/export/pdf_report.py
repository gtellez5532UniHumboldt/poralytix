"""
Generador de Reportes Científicos Auditables en formato PDF vectorial (PDF 1.4).
Incrusta tablas bioestadísticas completas y figuras de alta resolución (XObjects FlateDecode).
"""
import zlib
from pathlib import Path
from datetime import datetime
import numpy as np


def generate_native_pdf_report(pdf_path, title_str, subtitle_str, categories_list):
    """
    Genera un archivo PDF 1.4 multipágina estructurado por CATEGORÍAS DE ANÁLISIS,
    incrustando todas las tablas bioestadísticas y todas las figuras RGB.
    
    Parameters
    ----------
    pdf_path : str or Path
        Ruta de destino del reporte PDF
    title_str : str
        Título principal del reporte institucional
    subtitle_str : str
        Subtítulo, fecha y autor/laboratorio
    categories_list : list of dict
        Cada dict contiene:
          - category_title: str
          - description: str
          - tables: [(sec_title, headers, rows), ...]
          - images: [(img_caption, rgb_array_uint8), ...]
          
    Returns
    -------
    str
        Ruta final del archivo PDF generado
    """
    def esc_pdf(txt):
        s = str(txt).replace('\\', '\\\\').replace('(', '\\(').replace(')', '\\)')
        return s.encode('latin1', errors='replace').decode('latin1')

    pages_data = []
    cur_ops = []
    cur_imgs = []
    y_cursor = 752
    img_counter = 0

    def start_new_page():
        nonlocal cur_ops, cur_imgs, y_cursor
        if cur_ops or cur_imgs:
            pages_data.append(("\n".join(cur_ops), list(cur_imgs)))
        cur_ops = []
        cur_imgs = []
        y_cursor = 752

    # Encabezado institucional en la primera página
    cur_ops.append("0.0 0.17 0.29 rg")
    cur_ops.append("36 734 540 32 re f")
    cur_ops.append("1.0 1.0 1.0 rg")
    cur_ops.append(f"BT /F2 11.5 Tf 44 751 Td ({esc_pdf(title_str)}) Tj ET")
    cur_ops.append(f"BT /F1 7.5 Tf 44 739 Td ({esc_pdf(subtitle_str[:105])}) Tj ET")
    y_cursor = 716

    for cat in categories_list:
        if y_cursor < 160:
            start_new_page()

        # Banner de Categoría
        cur_ops.append("0.07 0.31 0.45 rg")
        cur_ops.append(f"36 {y_cursor - 16} 540 20 re f")
        cur_ops.append("1.0 1.0 1.0 rg")
        cur_ops.append(f"BT /F2 9.5 Tf 42 {y_cursor - 10} Td ({esc_pdf(cat['category_title'])}) Tj ET")
        y_cursor -= 24

        if cat.get('description'):
            cur_ops.append("0.2 0.2 0.2 rg")
            cur_ops.append(f"BT /F1 7.5 Tf 38 {y_cursor} Td ({esc_pdf(cat['description'][:115])}) Tj ET")
            y_cursor -= 14

        # Tablas
        for sec_title, headers, rows in cat.get('tables', []):
            if not rows:
                continue
            needed_h = 28 + min(len(rows), 15) * 12
            if y_cursor - needed_h < 48:
                start_new_page()

            cur_ops.append("0.0 0.17 0.29 rg")
            cur_ops.append(f"BT /F2 8.5 Tf 36 {y_cursor} Td ({esc_pdf(sec_title)}) Tj ET")
            y_cursor -= 13

            n_cols = max(1, len(headers))
            col_w = 540.0 / n_cols
            cur_ops.append("0.11 0.22 0.36 rg")
            cur_ops.append(f"36 {y_cursor - 9} 540 12 re f")
            cur_ops.append("1.0 1.0 1.0 rg")
            for c_i, h_txt in enumerate(headers):
                cx = 38 + c_i * col_w
                cur_ops.append(f"BT /F2 6.5 Tf {cx:.1f} {y_cursor - 6} Td ({esc_pdf(str(h_txt)[:20])}) Tj ET")
            y_cursor -= 12

            for r_i, row in enumerate(rows):
                if y_cursor < 46:
                    start_new_page()
                if r_i % 2 == 1:
                    cur_ops.append("0.95 0.96 0.98 rg")
                    cur_ops.append(f"36 {y_cursor - 9} 540 11 re f")
                cur_ops.append("0.1 0.1 0.1 rg")
                for c_i, val in enumerate(row):
                    cx = 38 + c_i * col_w
                    cur_ops.append(f"BT /F1 6.3 Tf {cx:.1f} {y_cursor - 6} Td ({esc_pdf(str(val)[:22])}) Tj ET")
                y_cursor -= 11
            y_cursor -= 10

        # Imágenes
        for img_caption, rgb_buf in cat.get('images', []):
            if rgb_buf is None:
                continue
            img_h_px, img_w_px, _ = rgb_buf.shape
            disp_w = 460.0
            disp_h = disp_w * (img_h_px / max(1.0, float(img_w_px)))
            if disp_h > 250.0:
                disp_h = 250.0
                disp_w = disp_h * (img_w_px / max(1.0, float(img_h_px)))

            if y_cursor - (disp_h + 26) < 42:
                start_new_page()

            cur_ops.append("0.0 0.17 0.29 rg")
            cur_ops.append(f"BT /F2 8 Tf 42 {y_cursor} Td ({esc_pdf(img_caption)}) Tj ET")
            y_cursor -= 8

            img_counter += 1
            im_name = f"Im{img_counter}"
            x_pos = 36.0 + (540.0 - disp_w) / 2.0
            y_pos = y_cursor - disp_h
            cur_ops.append(f"q {disp_w:.2f} 0 0 {disp_h:.2f} {x_pos:.2f} {y_pos:.2f} cm /{im_name} Do Q")
            cur_imgs.append((im_name, rgb_buf))
            y_cursor = y_pos - 14

        y_cursor -= 8

    if cur_ops or cur_imgs:
        pages_data.append(("\n".join(cur_ops), list(cur_imgs)))

    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"",
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>",
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>"
    ]

    page_obj_ids = []
    for p_content, p_imgs in pages_data:
        xobj_entries = []
        for im_name, rgb_arr in p_imgs:
            ih, iw, _ = rgb_arr.shape
            comp_bytes = zlib.compress(np.ascontiguousarray(rgb_arr, dtype=np.uint8).tobytes(), level=6)
            im_obj_id = len(objects) + 1
            im_header = (
                f"<< /Type /XObject /Subtype /Image /Width {iw} /Height {ih} "
                f"/ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode /Length {len(comp_bytes)} >>\nstream\n"
            ).encode('ascii')
            objects.append(im_header + comp_bytes + b"\nendstream")
            xobj_entries.append(f"/{im_name} {im_obj_id} 0 R")

        stream_bytes = p_content.encode('latin1', errors='replace')
        content_id = len(objects) + 1
        objects.append(f"<< /Length {len(stream_bytes)} >>\nstream\n".encode('ascii') + stream_bytes + b"\nendstream")

        page_id = len(objects) + 1
        page_obj_ids.append(page_id)
        xobj_dict_str = f" /XObject << {' '.join(xobj_entries)} >>" if xobj_entries else ""
        page_dict = (
            f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] "
            f"/Resources << /Font << /F1 3 0 R /F2 4 0 R >>{xobj_dict_str} >> /Contents {content_id} 0 R >>"
        ).encode('ascii')
        objects.append(page_dict)

    kids_str = " ".join(f"{pid} 0 R" for pid in page_obj_ids)
    objects[1] = f"<< /Type /Pages /Kids [{kids_str}] /Count {len(page_obj_ids)} >>".encode('ascii')

    out = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offsets = []
    for idx, obj_bytes in enumerate(objects, start=1):
        offsets.append(len(out))
        out.extend(f"{idx} 0 obj\n".encode('ascii'))
        out.extend(obj_bytes)
        out.extend(b"\nendobj\n")

    xref_pos = len(out)
    out.extend(f"xref\n0 {len(objects) + 1}\n0000000000 65535 f \n".encode('ascii'))
    for off in offsets:
        out.extend(f"{off:010d} 00000 n \n".encode('ascii'))
    out.extend(
        f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref_pos}\n%%EOF\n".encode('ascii')
    )

    target_pdf = Path(pdf_path)
    try:
        with open(target_pdf, "wb") as f:
            f.write(out)
    except PermissionError:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        target_pdf = target_pdf.with_name(f"{target_pdf.stem}_{ts}{target_pdf.suffix}")
        with open(target_pdf, "wb") as f:
            f.write(out)
    return str(target_pdf)
