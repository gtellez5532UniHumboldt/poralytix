#!/usr/bin/env python3
"""
Punto de Entrada Oficial (GUI) - Estación de Trabajo de Cinética de Citometría de Flujo (.fcs).
Ejecutar con:
    py -3.11 main.py
    py -3.11 main.py --lang en
"""
import sys
import os
import argparse
import traceback
import tkinter as tk
from tkinter import messagebox

# Asegurar que el directorio de Refactorizacion esté en el PYTHONPATH
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)

from cytoworkstation.gui.app import FlowCytometryKineticWorkstation
from cytoworkstation.i18n import SUPPORTED_LANGUAGES


def main():
    parser = argparse.ArgumentParser(description="Estación de Trabajo de Cinética de Citometría de Flujo (.fcs)")
    parser.add_argument("--lang", choices=list(SUPPORTED_LANGUAGES.keys()), default="es",
                        help="Idioma de la interfaz gráfica (es: Español, en: English)")
    args, _ = parser.parse_known_args()

    try:
        root = tk.Tk()
        app = FlowCytometryKineticWorkstation(root)
        if args.lang and args.lang != "es":
            app.set_language(args.lang)
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
