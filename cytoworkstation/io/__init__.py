"""
Módulo I/O: Lectura binaria FCS, extracción de metadatos y gestión de sesiones.
"""
from .fcs_parser import parse_fcs_native, resolve_channel
from .metadata import guess_metadata_from_filename
from .session import save_session_dict, load_session_dict

__all__ = [
    'parse_fcs_native',
    'resolve_channel',
    'guess_metadata_from_filename',
    'save_session_dict',
    'load_session_dict'
]
