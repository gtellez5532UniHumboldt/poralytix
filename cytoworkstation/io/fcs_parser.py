"""
Parser binario nativo de archivos FCS (Flow Cytometry Standard 2.0, 3.0 y 3.1).
No requiere librerías externas de citometría como FlowKit o fcsparser.
"""
import os
import numpy as np
import pandas as pd


def parse_fcs_native(filepath):
    """
    Lee archivos binarios .fcs sin librerías externas.
    
    Parameters
    ----------
    filepath : str or Path
        Ruta absoluta o relativa al archivo .fcs
        
    Returns
    -------
    meta : dict
        Diccionario con las claves/tokens del segmento TEXT del FCS
    df_events : pd.DataFrame
        DataFrame conteniendo los eventos medidos por canal
    """
    filepath = str(filepath)
    with open(filepath, 'rb') as f:
        header = f.read(58).decode('ascii', errors='ignore')
        if not header.startswith('FCS'):
            raise ValueError(f"Cabecera FCS inválida en {os.path.basename(filepath)}")

        text_start = int(header[10:18].strip())
        text_end = int(header[18:26].strip())
        data_start_hdr = int(header[26:34].strip()) if header[26:34].strip() else 0

        f.seek(text_start)
        text_seg = f.read(text_end - text_start + 1)
        delimiter = chr(text_seg[0])
        tokens = text_seg.decode('latin1', errors='ignore').split(delimiter)

        meta = {}
        for i in range(1, len(tokens) - 1, 2):
            meta[tokens[i].strip().upper()] = tokens[i + 1].strip()

        data_start = int(meta.get('$DATASTART', data_start_hdr))
        total_events = int(meta.get('$TOT', 0))
        num_channels = int(meta.get('$PAR', 0))
        datatype = meta.get('$DATATYPE', 'F').upper()
        byte_order = meta.get('$BYTEORD', '4,3,2,1')

        channel_names = []
        for p in range(1, num_channels + 1):
            p_short = meta.get(f'$P{p}S', '').strip()
            p_name = meta.get(f'$P{p}N', '').strip()
            if p_name and any(k in p_name.upper() for k in ('FITC', 'FL1', 'FSC', 'SSC', 'TIME', 'BL1')):
                channel_names.append(p_name)
            elif p_short:
                channel_names.append(p_short)
            elif p_name:
                channel_names.append(p_name)
            else:
                channel_names.append(f'Param_{p}')

        endian = '>' if byte_order in ('4,3,2,1', '2,1') else '<'
        if datatype == 'F':
            dtype_str = f'{endian}f4'
        elif datatype == 'D':
            dtype_str = f'{endian}f8'
        elif datatype == 'I':
            bits = int(meta.get('$P1B', 16))
            dtype_str = f'{endian}u4' if bits == 32 else f'{endian}u2'
        else:
            dtype_str = f'{endian}f4'

        f.seek(data_start)
        count_expected = total_events * num_channels
        raw_arr = np.fromfile(f, dtype=dtype_str, count=count_expected)
        if len(raw_arr) < count_expected and num_channels > 0:
            total_events = len(raw_arr) // num_channels
            raw_arr = raw_arr[:total_events * num_channels]

        raw_arr = np.ascontiguousarray(raw_arr, dtype=np.float64)
        df_events = pd.DataFrame(raw_arr.reshape((total_events, num_channels)), columns=channel_names)
        return meta, df_events


def get_fcs_channels(filepath):
    """Obtiene la lista de nombres de canales disponibles en un archivo .fcs."""
    try:
        _, df_ev = parse_fcs_native(filepath)
        return [str(c).strip() for c in df_ev.columns]
    except Exception:
        return []


def resolve_perm_channel(columns, preferred_channel="FITC (Por defecto)"):
    """
    Resuelve el canal de fluorescencia para permeabilidad (prioriza FITC / FL1).
    """
    cols = {str(c).upper().strip(): c for c in columns}
    pref = (preferred_channel or "").strip().upper()

    if pref and not pref.startswith("FITC (") and not pref.startswith("AUTO"):
        if pref in cols:
            return cols[pref]
        for k, orig in cols.items():
            if pref in k:
                return orig

    priority_list = ['FITC-A', 'FL1-A', 'FITC-H', 'FL1-H', 'FITC', 'FL1', 'BL1-A', 'BL1-H', 'GFP-A']
    for cand in priority_list:
        if cand in cols:
            return cols[cand]

    for k, orig in cols.items():
        if any(tok in k for tok in ('FITC', 'FL1', 'CALCEIN', 'BL1', 'GFP')):
            return orig

    return None


def resolve_channel(columns, candidates):
    """Busca el primer canal que coincida con la lista de candidatos prioritarios."""
    cols = {str(c).upper().strip(): c for c in columns}
    for cand in candidates:
        cu = cand.upper()
        if cu in cols:
            return cols[cu]
    for cand in candidates:
        cu = cand.upper()
        for k, orig in cols.items():
            if cu in k:
                return orig
    return None
