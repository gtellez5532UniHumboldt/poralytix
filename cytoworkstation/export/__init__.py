"""
Módulo de exportación: Tablas científicas CSV/TSV, figuras 600 DPI y reporte PDF nativo.
"""
from .tables_export import export_analysis_tables
from .plots_600dpi import (
    save_rgb_buffer_as_png,
    save_rgb_buffer_as_tiff,
    save_mpl_figure_600dpi_tiff,
    render_highres_plot_png,
    render_highres_histogram_png
)
from .pdf_report import generate_native_pdf_report

__all__ = [
    'export_analysis_tables',
    'save_rgb_buffer_as_png',
    'save_rgb_buffer_as_tiff',
    'save_mpl_figure_600dpi_tiff',
    'render_highres_plot_png',
    'render_highres_histogram_png',
    'generate_native_pdf_report'
]
