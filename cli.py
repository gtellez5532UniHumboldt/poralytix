#!/usr/bin/env python3
"""
Interfaz de Línea de Comandos (CLI / Modo Batch Headless)
Permite procesar directorios completos de archivos .fcs (como CZS5c-, 813, 820, sat1)
sin necesidad de abrir la interfaz gráfica Tkinter.

Uso:
    py -3.11 cli.py --input-dir ../CZS5c- --output-dir ../Resultados_CZS5 --target-cell THP-1 --control-cell PBMCs
"""
import sys
import os
import argparse
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)

from cytoworkstation.config import PEPTIDE_ALIASES
from cytoworkstation.io.fcs_parser import parse_fcs_native, resolve_perm_channel
from cytoworkstation.io.metadata import guess_metadata_from_filename
from cytoworkstation.core.kinetics import extract_tube_kinetics, extract_slope_based_parameters, compute_categorical_kinetics_metrics
from cytoworkstation.core.models import fit_dose_response_model
from cytoworkstation.export.tables_export import export_analysis_tables
from cytoworkstation.export.pdf_report import generate_native_pdf_report


def run_batch_analysis(
    input_dir,
    output_dir,
    target_cell="THP-1",
    control_cell="PBMCs",
    perm_channel="FITC (Por defecto)",
    model_name="Allosteric Hill",
    fsc_min=50000.0,
    ssc_min=10000.0,
    v1l_excl=3.00,
    v2l_cutoff=4.85,
    bin_size=15
):
    in_path = Path(input_dir)
    if not in_path.exists():
        raise FileNotFoundError(f"El directorio de entrada no existe: {in_path}")

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    fcs_files = sorted(list(in_path.glob("**/*.fcs")) + list(in_path.glob("**/*.FCS")))
    if not fcs_files:
        print(f"No se encontraron archivos .fcs en {in_path}")
        return

    print(f"\n=======================================================")
    print(f"  CITOMETRÍA CINÉTICA — PROCESAMIENTO BATCH HEADLESS")
    print(f"=======================================================")
    print(f"  • Directorio Origen:  {in_path.resolve()}")
    print(f"  • Archivos .fcs:      {len(fcs_files)}")
    print(f"  • Célula Diana:       {target_cell}")
    print(f"  • Célula Control:     {control_cell}")
    print(f"  • Salida:             {out_path.resolve()}")
    print(f"=======================================================\n")

    # 1. Extracción de Metadatos
    records = []
    for fp in fcs_files:
        meta = guess_metadata_from_filename(str(fp), target_cell, control_cell)
        records.append(meta)

    df_meta = pd.DataFrame(records)
    print("Muestras detectadas:")
    for _, r in df_meta.iterrows():
        print(f"  - {r['Filename']}: Compuesto={r['Compound']}, Célula={r['Cell_Line']}, Dosis={r['Concentration_uM']}µM, Rol={r['Sample_Role']}")

    # 2. Extracción Cinética por Tubo
    kinetics_rows = []
    full_event_rows = []

    print("\nProcesando cinéticas por tubo...")
    for meta in records:
        fp = meta['Filepath']
        is_cc = "CC" in meta['Sample_Role']
        try:
            binned_df, cache = extract_tube_kinetics(
                fp, is_cc, bin_size=bin_size,
                fsc_min=fsc_min, ssc_min=ssc_min,
                perm_channel=perm_channel,
                log_excl_v1l=v1l_excl,
                log_cutoff=v2l_cutoff
            )
            for _, brow in binned_df.iterrows():
                row_dict = dict(brow)
                row_dict.update({
                    'Filename': meta['Filename'],
                    'Compound': meta['Compound'],
                    'Cell_Line': meta['Cell_Line'],
                    'Sample_Role': meta['Sample_Role'],
                    'Concentration_uM': meta['Concentration_uM'],
                    'Replicate': meta['Replicate']
                })
                full_event_rows.append(row_dict)

            # Extraer pendiente y parámetros si es experimental
            if not is_cc and len(binned_df) >= 2:
                t_arr = binned_df['Time_Window'].values
                p_arr = binned_df['Affected_V2L_Pct'].values
                tlag, tmax, v_act, max_p = extract_slope_based_parameters(
                    t_arr, p_arr, is_control_cell=(meta['Cell_Line'] == control_cell),
                    start_time=meta['Exclusion_s']
                )
                cat_metrics = compute_categorical_kinetics_metrics(
                    t_arr, p_arr, binned_df['Median_FITC'].values,
                    start_time=meta['Exclusion_s']
                )
                kin_rec = dict(meta)
                kin_rec.update({
                    'T_lag': tlag,
                    'T_max': tmax,
                    'V_active': v_act,
                    'Max_Lysis_Pct': max_p,
                    'V2L_Cutoff_Log10': v2l_cutoff,
                    'V1L_Excl_Log10': v1l_excl
                })
                kin_rec.update(cat_metrics)
                kinetics_rows.append(kin_rec)
        except Exception as err:
            print(f"  [Error procesando {meta['Filename']}]: {err}")

    df_kinetics = pd.DataFrame(kinetics_rows)
    df_full = pd.DataFrame(full_event_rows)

    # 3. Ajuste de Modelos Matemáticos por Compuesto y Línea Celular
    saturation_records = []
    if not df_kinetics.empty:
        for (pep, cell), grp in df_kinetics.groupby(['Compound', 'Cell_Line']):
            grp_dose = grp.groupby('Concentration_uM').agg(
                v_mean=('V_active', 'mean'),
                v_sem=('V_active', 'sem'),
                max_p=('Max_Lysis_Pct', 'mean')
            ).reset_index().sort_values('Concentration_uM')

            if len(grp_dose) >= 2:
                fit_res = fit_dose_response_model(
                    grp_dose['Concentration_uM'].values,
                    grp_dose['v_mean'].values,
                    analyzer_name=model_name,
                    raw_x=grp['Concentration_uM'].values,
                    raw_y=grp['V_active'].values
                )
                if fit_res:
                    saturation_records.append({
                        'Compound': pep,
                        'Cell_Line': cell,
                        'Model_Used': fit_res['model_lbl'],
                        'V_max': fit_res['vmax'],
                        'K_d_uM': fit_res['kd'],
                        'n_H': fit_res['nh'],
                        'R2_Group': fit_res['r2_grp'],
                        'R2_Replicates': fit_res['r2_rep'],
                        'RMSE': fit_res['rmse'],
                        'Equation': fit_res['eq_eval']
                    })

    df_saturation = pd.DataFrame(saturation_records)

    # 4. Exportar Tablas CSV / TSV
    tables_to_export = {
        'Metadata_Samples_Summary': df_meta,
        'Kinetics_Summary_Parameters': df_kinetics,
        'Model_Saturation_Parameters': df_saturation,
        'Full_Binned_Kinetics': df_full
    }
    exp_files = export_analysis_tables(tables_to_export, out_path)
    print(f"\nTablas exportadas ({len(exp_files)} archivos en total).")

    # 5. Generar Reporte PDF
    pdf_out = out_path / f"Auditoria_Cinetica_{in_path.name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf"
    cat_tables = []
    if not df_saturation.empty:
        cat_tables.append((
            "Ajustes de Modelos Farmacodinámicos",
            list(df_saturation.columns),
            df_saturation.round(4).astype(str).values.tolist()
        ))
    if not df_kinetics.empty:
        cols_k = ['Compound', 'Cell_Line', 'Concentration_uM', 'T_lag', 'T_max', 'V_active', 'Affectation_State']
        avail_k = [c for c in cols_k if c in df_kinetics.columns]
        cat_tables.append((
            "Parámetros Cinéticos y Categóricos",
            avail_k,
            df_kinetics[avail_k].round(3).astype(str).values.tolist()
        ))

    categories_list = [
        {
            'category_title': f"Categoría 0: Muestras Procesadas ({in_path.name})",
            'description': f"Total de muestras FCS procesadas en modo batch: {len(fcs_files)}.",
            'tables': [("Listado de Muestras y Metadatos", list(df_meta.columns), df_meta.astype(str).values.tolist())],
            'images': []
        },
        {
            'category_title': "Categoría I: Resultados Farmacodinámicos y Cinéticos",
            'description': "Parámetros de ajuste no lineal de modelos y tasas de permeabilización calculadas.",
            'tables': cat_tables,
            'images': []
        }
    ]

    saved_pdf = generate_native_pdf_report(
        pdf_out,
        f"AUDITORÍA CINÉTICA AUTOMATIZADA: {in_path.name.upper()}",
        f"Generado en Modo Batch Headless el {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        categories_list
    )
    print(f"Reporte PDF generado exitosamente:\n  -> {saved_pdf}\n")
    print("[OK] Procesamiento batch finalizado con exito.")


def main():
    parser = argparse.ArgumentParser(
        description="Analizador de Cinética de Citometría de Flujo - Modo Batch Headless"
    )
    parser.add_argument("--input-dir", "-i", required=True, help="Directorio con archivos .fcs")
    parser.add_argument("--output-dir", "-o", default="Resultados_Exportados", help="Directorio destino")
    parser.add_argument("--target-cell", default="THP-1", help="Línea celular blanco diana (default: THP-1)")
    parser.add_argument("--control-cell", default="PBMCs", help="Línea celular sana control (default: PBMCs)")
    parser.add_argument("--model", default="Allosteric Hill", help="Modelo matemático a ajustar")
    parser.add_argument("--v2l-cutoff", type=float, default=4.85, help="Cut-Off V2L en log10")
    parser.add_argument("--v1l-excl", type=float, default=3.00, help="Cut-Off V1L de exclusión en log10")

    args = parser.parse_args()
    run_batch_analysis(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        target_cell=args.target_cell,
        control_cell=args.control_cell,
        model_name=args.model,
        v2l_cutoff=args.v2l_cutoff,
        v1l_excl=args.v1l_excl
    )


if __name__ == "__main__":
    main()
