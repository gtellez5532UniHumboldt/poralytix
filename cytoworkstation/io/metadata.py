"""
Módulo de deducción y edición de metadatos experimentales a partir de nombres de archivo y carpetas.
Resuelve de forma robusta Compuesto, Línea Celular, Concentración, Rol (CC vs EXP) y Réplica.
"""
import os
import re
import numpy as np
from ..config import PEPTIDE_ALIASES


KNOWN_CELL_NAMES = {'pbmc', 'pbmcs', 'thp1', 'thp-1', 'erythrocyte', 'erythrocytes', 'rbc', 'jurkat', 'k562'}
GENERIC_FOLDERS = {
    'data', 'fcs', 'files', 'desktop', 'downloads', 'demo_fcs_data',
    'muestras', '', 'pbmc', 'pbmcs', 'thp1', 'thp-1'
}


def guess_metadata_from_filename(filepath, default_target="THP-1", default_control="PBMCs"):
    """
    Deduce Compuesto, Línea Celular, Rol (CC vs Experimental), Concentración (µM),
    Réplica (R1..Rn) y Ventana de Exclusión Inicial desde el nombre del archivo y la estructura de carpetas.
    
    Parameters
    ----------
    filepath : str
        Ruta completa del archivo .fcs
    default_target : str
        Línea celular diana/tumoral por defecto (ej. THP-1)
    default_control : str
        Línea celular sana/control por defecto (ej. PBMCs)
        
    Returns
    -------
    dict
        Metadatos parseados con tipos normalizados
    """
    fname = os.path.basename(filepath)
    parent_dir = os.path.dirname(filepath)
    folder = os.path.basename(parent_dir)
    grandparent_dir = os.path.basename(os.path.dirname(parent_dir)) if parent_dir else ""
    base = os.path.splitext(fname)[0]
    parts = re.split(r'[-_ ]+', base)
    upper_tokens = [p.upper() for p in parts]

    # 1. Determinación de Línea Celular
    base_upper = base.upper()
    folder_upper = folder.upper()
    if 'PBMC' in base_upper or 'PBMC' in folder_upper or (default_control and default_control.upper() in base_upper):
        cell = default_control
    elif 'THP' in base_upper or 'THP' in folder_upper or (default_target and default_target.upper() in base_upper):
        cell = default_target
    else:
        cell = default_target

    # 2. Rol: Control Basal CC (-60 a 0 s) vs Experimental (0 a 600 s)
    is_cc = (
        any(k in upper_tokens for k in ('CC', 'BASE', 'BASELINE', 'CTRL', 'PRE')) or
        'BASE' in base_upper or
        '_CC_' in fname or
        '-CC-' in fname or
        fname.endswith('_CC.fcs')
    )
    role = 'CC (Baseline -60 to 0s)' if is_cc else 'Experimental (0 to 600s)'

    # 3. Determinación del Compuesto
    compound = None
    # Revisar tokens del nombre del archivo primero si coinciden con alias conocidos
    for p in parts:
        p_clean = p.lower()
        if p_clean in PEPTIDE_ALIASES:
            compound = PEPTIDE_ALIASES[p_clean]
            break
        # Manejo de casos como CZS5, CZS5c, CZS-5
        if re.match(r'^czs\d+.*$', p_clean):
            compound = PEPTIDE_ALIASES.get(p_clean, p.upper())
            break

    # Si no se encontró en tokens del archivo, revisar carpeta contenedora
    if not compound:
        folder_clean = folder.lower()
        if folder_clean in PEPTIDE_ALIASES:
            compound = PEPTIDE_ALIASES[folder_clean]
        elif folder_clean not in GENERIC_FOLDERS and folder_clean not in KNOWN_CELL_NAMES:
            compound = PEPTIDE_ALIASES.get(folder_clean, folder)
        elif grandparent_dir and grandparent_dir.lower() in PEPTIDE_ALIASES:
            compound = PEPTIDE_ALIASES[grandparent_dir.lower()]
        elif grandparent_dir and grandparent_dir.lower() not in GENERIC_FOLDERS and grandparent_dir.lower() not in KNOWN_CELL_NAMES:
            compound = PEPTIDE_ALIASES.get(grandparent_dir.lower(), grandparent_dir)

    # Fallback si aún no se deduce
    if not compound:
        for p in parts:
            if p.upper() not in ('PBMC', 'PBMCS', 'THP1', 'THP', 'CC', 'BASE', 'CTRL', 'EXP'):
                if not re.match(r'^\d+(\.\d+)?$', p):
                    compound = PEPTIDE_ALIASES.get(p.lower(), p)
                    break
        if not compound:
            compound = parts[0] if parts else "Compound_A"

    # 4. Concentración molar (µM)
    base_no_cell = re.sub(r'THP-?1|PBMC[s]?|Act8-13|Act8-20|Satanin[-_ ]?1|CZS5[c-]?', '', base, flags=re.IGNORECASE)
    parts_conc = [t for t in re.split(r'[-_ ]+', base_no_cell) if t]
    conc = 0.0

    # Priorizar números con decimales explícitos (ej. 0.5, 0.70, 1.08, 1.40, 2.16, 2.80, 4.33, 4.35, 5.70)
    for p in parts_conc:
        if re.match(r'^\d+\.\d+$', p):
            try:
                val = float(p)
                if val < 500:
                    conc = val
                    break
            except ValueError:
                pass

    # Si no hay decimal explícito, buscar enteros excluyendo réplicas y números identificadores
    if conc == 0.0:
        candidates = parts_conc[:-1] if len(parts_conc) > 1 else parts_conc
        for p in candidates:
            if re.match(r'^\d+$', p):
                try:
                    val = float(p)
                    if val < 500 and p not in ('813', '820'):
                        conc = val
                        break
                except ValueError:
                    pass

    # 5. Réplica (R1..Rn)
    rep = "R1"
    # Si el último token es un número o R1..R9
    if parts:
        last = parts[-1]
        m = re.match(r'^[rR]?(\d+)$', last)
        if m:
            rep = f"R{m.group(1)}"
        else:
            for p in reversed(parts):
                m = re.match(r'^[rR](\d+)$', p)
                if m:
                    rep = f"R{m.group(1)}"
                    break

    # 6. Ventana de Exclusión Inicial (s) para artefactos tempranos de inyección
    excl_s = 0.0
    if cell == default_target:
        if compound in ('Satanin 1', 'sat1') and 0.0 < conc <= 2.80:
            excl_s = 60.0
        elif compound in ('Act8-20', '820') and np.isclose(conc, 4.35, atol=0.1):
            excl_s = 45.0

    return {
        'Filepath': str(filepath),
        'Filename': fname,
        'Compound': str(compound),
        'Cell_Line': str(cell),
        'Sample_Role': str(role),
        'Concentration_uM': float(conc),
        'Replicate': str(rep),
        'Exclusion_s': float(excl_s)
    }
