"""
Modal para la edición individual de umbrales morfológicos y Cut-Offs por muestra.
Soporte bilingüe Español e Inglés.
"""
import tkinter as tk
from tkinter import ttk, messagebox

from ...i18n import get_current_language, t


def open_sample_cutoff_editor_dialog(parent, row_id, sample_name, init_values, on_save_callback, lang=None):
    """
    Abre diálogo modal para ajustar Cut-Offs y ventanas específicas de una muestra experimental.
    
    Parameters
    ----------
    parent : tk.Widget
        Ventana padre
    row_id : str
        Identificador de la fila en la tabla
    sample_name : str
        Nombre identificador de la muestra
    init_values : dict or list
        Valores iniciales [v1l, v2l, fsc, ssc, tlag, tmax]
    on_save_callback : callable
        Función a invocar al guardar con los 6 valores numéricos
    lang : str, optional
        Idioma de la interfaz ('es' o 'en')
    """
    active_lang = lang or get_current_language()
    dlg = tk.Toplevel(parent)
    dlg.title(t('editor_title', lang=active_lang).format(sample_name))
    dlg.geometry("450x330")
    dlg.transient(parent)
    dlg.grab_set()

    frm = ttk.LabelFrame(dlg, text=t('editor_group', lang=active_lang).format(sample_name), padding=10)
    frm.pack(fill="both", expand=True, padx=10, pady=8)

    fields = [
        (t('editor_v1l', lang=active_lang), init_values[0]),
        (t('editor_v2l', lang=active_lang), init_values[1]),
        (t('editor_fsc', lang=active_lang), init_values[2]),
        (t('editor_ssc', lang=active_lang), init_values[3]),
        (t('editor_tlag', lang=active_lang), init_values[4]),
        (t('editor_tmax', lang=active_lang), init_values[5]),
    ]
    entries = []
    for idx_f, (lbl_txt, init_v) in enumerate(fields):
        ttk.Label(frm, text=lbl_txt, font=("DejaVu Sans", 8, "bold")).grid(row=idx_f, column=0, sticky="w", pady=4, padx=4)
        ent = ttk.Entry(frm, width=14)
        ent.insert(0, str(init_v))
        ent.grid(row=idx_f, column=1, sticky="ew", pady=4, padx=4)
        entries.append(ent)

    def save_modal():
        try:
            v1l_v = float(entries[0].get().strip())
            v2l_v = float(entries[1].get().strip())
            fsc_v = float(entries[2].get().strip())
            ssc_v = float(entries[3].get().strip())
            tlag_v = float(entries[4].get().strip())
            tmax_v = float(entries[5].get().strip())
        except ValueError as err:
            err_title = "Valor Inválido" if active_lang == 'es' else "Invalid Value"
            err_msg = f"Todos los campos deben ser numéricos:\n{err}" if active_lang == 'es' else f"All fields must be numerical:\n{err}"
            messagebox.showerror(err_title, err_msg, parent=dlg)
            return
        dlg.destroy()
        if on_save_callback:
            on_save_callback(row_id, v1l_v, v2l_v, fsc_v, ssc_v, tlag_v, tmax_v)

    btn_f = ttk.Frame(frm)
    btn_f.grid(row=len(fields), column=0, columnspan=2, pady=10)
    tk.Button(
        btn_f,
        text=t('editor_btn_save', lang=active_lang),
        bg="#117864", fg="white",
        font=("DejaVu Sans", 8, "bold"), padx=10,
        command=save_modal
    ).pack(side="left", padx=4)
    ttk.Button(
        btn_f,
        text=t('editor_btn_cancel', lang=active_lang),
        command=dlg.destroy
    ).pack(side="left", padx=4)
    return dlg
