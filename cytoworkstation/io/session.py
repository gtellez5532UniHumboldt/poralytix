"""
Módulo de serialización y deserialización de estado de sesiones de trabajo (.json).
Permite persistir parámetros globales, archivos, metadatos, cut-offs y ventanas de análisis.
"""
import json
from pathlib import Path
from datetime import datetime


def save_session_dict(state_dict, fpath):
    """
    Guarda un diccionario de estado de sesión a archivo JSON en UTF-8.
    
    Parameters
    ----------
    state_dict : dict
        Diccionario con la configuración completa de la sesión
    fpath : str or Path
        Ruta de destino del archivo .json
    """
    fpath = Path(fpath)
    fpath.parent.mkdir(parents=True, exist_ok=True)
    if 'session_version' not in state_dict:
        state_dict['session_version'] = '2.0'
    if 'saved_at' not in state_dict:
        state_dict['saved_at'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    fpath.write_text(json.dumps(state_dict, indent=2, ensure_ascii=False), encoding="utf-8")
    return str(fpath)


def load_session_dict(fpath):
    """
    Lee y valida un archivo JSON de sesión de trabajo.
    
    Parameters
    ----------
    fpath : str or Path
        Ruta al archivo .json
        
    Returns
    -------
    dict
        Diccionario de estado deserializado
    """
    fpath = Path(fpath)
    if not fpath.exists():
        raise FileNotFoundError(f"No existe el archivo de sesión: {fpath}")
    raw_text = fpath.read_text(encoding="utf-8")
    state_dict = json.loads(raw_text)
    return state_dict
