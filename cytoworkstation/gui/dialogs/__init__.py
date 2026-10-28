"""
Diálogos modales y ventanas auxiliares para la interfaz gráfica.
"""
from .math_help import open_mathematical_models_help_dialog
from .user_guide import open_user_guide_dialog
from .sample_editor import open_sample_cutoff_editor_dialog

__all__ = [
    'open_mathematical_models_help_dialog',
    'open_user_guide_dialog',
    'open_sample_cutoff_editor_dialog'
]
