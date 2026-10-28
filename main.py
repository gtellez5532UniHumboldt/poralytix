#!/usr/bin/env python3
"""
Punto de Entrada Oficial (GUI) - Estación de Trabajo de Cinética de Citometría de Flujo (.fcs).
Ejecutar con:
    py -3.11 main.py
"""
import sys
import os
import traceback
import tkinter as tk
from tkinter import messagebox

# Asegurar que el directorio de Refactorizacion esté en el PYTHONPATH
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)

from cytoworkstation.gui.app import FlowCytometryKineticWorkstation


def main():
    try:
        root = tk.Tk()
        app = FlowCytometryKineticWorkstation(root)
        root.mainloop()
    except Exception as exc:
        err_msg = traceback.format_exc()
        print(f"\n[ERROR CRÍTICO AL INICIAR LA APLICACIÓN]:\n{err_msg}", file=sys.stderr)
        try:
            r_err = tk.Tk()
            r_err.withdraw()
            messagebox.showerror(
                "Error al Iniciar la Estación de Trabajo",
                f"Ocurrió un error inesperado al iniciar la interfaz:\n\n{exc}\n\n{err_msg}"
            )
            r_err.destroy()
        except Exception:
            pass
        sys.exit(1)


if __name__ == "__main__":
    main()
