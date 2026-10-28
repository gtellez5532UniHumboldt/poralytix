"""
Widgets interactivos para la interfaz gráfica científica en Tkinter.
"""
from .canvas_histogram import NativeHistogramCutoffCanvas
from .canvas_kinetics import NativeInteractiveMultiPlot
from .canvas_saturation import NativeSaturationPreview
from .canvas_categorical import NativeCategoricalProportionsCanvas

__all__ = [
    'NativeHistogramCutoffCanvas',
    'NativeInteractiveMultiPlot',
    'NativeSaturationPreview',
    'NativeCategoricalProportionsCanvas'
]
