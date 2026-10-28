"""
Modal para la edición individual de umbrales morfológicos y Cut-Offs por muestra.
"""
import tkinter as tk
from tkinter import ttk, messagebox


def open_sample_cutoff_editor_dialog(parent, row_id, sample_name, init_values, on_save_callback):
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
    """
    dlg = tk.Toplevel(parent)
    dlg.title(f"Editar Cut-Offs de Muestra: {sample_name}")
    dlg.geometry("420x310")
    dlg.transient(parent)
    dlg.grab_set()

    frm = ttk.LabelFrame(dlg, text=f" Parámetros Individuales: {sample_name} ", padding=10)
    frm.pack(fill="both", expand=True, padx=10, pady=8)

    fields = [
        ("Exclusión Detritos V1L (< log10):", init_values[0]),
        ("Cut-Off Afectados V2L/V2R (log10):", init_values[1]),
        ("Umbral Morfológico FSC-A Mín (>):", init_values[2]),
        ("Umbral Morfológico SSC-A Mín (>):", init_values[3]),
        ("Inicio Fase Activa T_lag (s):", init_values[4]),
        ("Fin Fase Activa / Meseta T_max (s):", init_values[5]),
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
            messagebox.showerror("Valor Inválido", f"Todos los campos deben ser numéricos:\n{err}", parent=dlg)
            return
        dlg.destroy()
        if on_save_callback:
            on_save_callback(row_id, v1l_v, v2l_v, fsc_v, ssc_v, tlag_v, tmax_v)

    btn_f = ttk.Frame(frm)
    btn_f.grid(row=len(fields), column=0, columnspan=2, pady=10)
    tk.Button(btn_f, text="✅ Aplicar a Esta Muestra", bg="#117864", fg="white", font=("DejaVu Sans", 8, "bold"), padx=10, command=save_modal).pack(side="left", padx=4)
    ttk.Button(btn_f, text="Cancelar", command=dlg.destroy).pack(side="left", padx=4)
    return dlg
