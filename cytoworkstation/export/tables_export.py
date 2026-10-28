"""
Módulo de exportación de tablas científicas en formatos CSV y TSV (Tab-Separated Values).
Exporta las tablas estructuradas de las Categorías 0 a V del análisis de citometría.
"""
from pathlib import Path
from datetime import datetime
import pandas as pd


def export_analysis_tables(tables_dict, output_dir):
    """
    Exporta un diccionario de DataFrames a CSV y TSV tabulados con manejo de errores de permisos.
    
    Parameters
    ----------
    tables_dict : dict
        Diccionario { 'Nombre_Tabla': pd.DataFrame }
    output_dir : str or Path
        Directorio destino para los archivos generados
        
    Returns
    -------
    list
        Lista de rutas completas de los archivos exportados exitosamente
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    exported_files = []

    for name, df in tables_dict.items():
        if df is None:
            continue
        df_out = df.copy() if isinstance(df, pd.DataFrame) else pd.DataFrame(df)
        if df_out.empty and len(df_out.columns) == 0:
            continue

        csv_p = out_dir / f"{name}.csv"
        tsv_p = out_dir / f"{name}.tsv"

        try:
            df_out.to_csv(csv_p, index=False)
        except PermissionError:
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            csv_p = out_dir / f"{name}_{ts}.csv"
            df_out.to_csv(csv_p, index=False)

        try:
            df_out.to_csv(tsv_p, sep="\t", index=False)
        except PermissionError:
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            tsv_p = out_dir / f"{name}_{ts}.tsv"
            df_out.to_csv(tsv_p, sep="\t", index=False)

        exported_files.extend([str(csv_p), str(tsv_p)])

    return exported_files
