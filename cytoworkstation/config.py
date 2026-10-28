"""
Módulo de Configuración, Paletas Cromáticas y Constantes Biofísicas.
"""
import sys

AVOGADRO = 6.02214076e23

# Mapeo canónico de nombres y alias de compuestos/péptidos
PEPTIDE_ALIASES = {
    'sat1': 'Satanin 1',
    'satanin 1': 'Satanin 1',
    '813': 'Act8-13',
    'act8-13': 'Act8-13',
    '820': 'Act8-20',
    'act8-20': 'Act8-20',
    'czs5': 'CZS5',
    'czs5c': 'CZS5',
    'czs5c-': 'CZS5'
}

# Rampas cromáticas para publicaciones científicas (600 DPI)
DOSE_RAMPS = {
    'Satanin 1': ['#85C1E9', '#3498DB', '#1B4F72', '#0B172A', '#050A14'],
    'sat1': ['#85C1E9', '#3498DB', '#1B4F72', '#0B172A', '#050A14'],
    'Act8-13': ['#F9E79F', '#F39C12', '#B9770E', '#5F3A04', '#3B2202'],
    '813': ['#F9E79F', '#F39C12', '#B9770E', '#5F3A04', '#3B2202'],
    'Act8-20': ['#F1948A', '#E74C3C', '#922B21', '#4A110A', '#2B0804'],
    '820': ['#F1948A', '#E74C3C', '#922B21', '#4A110A', '#2B0804'],
    'CZS5': ['#A3E4D7', '#48C9B0', '#1ABC9C', '#117864', '#0B4F41'],
    'czs5': ['#A3E4D7', '#48C9B0', '#1ABC9C', '#117864', '#0B4F41'],
    'czs5c-': ['#A3E4D7', '#48C9B0', '#1ABC9C', '#117864', '#0B4F41']
}

GENERIC_RAMPS = [
    ['#5DADE2', '#2E86C1', '#1B4F72', '#0B172A', '#050A14'],
    ['#F5B041', '#F39C12', '#B9770E', '#7E5109', '#4A2F04'],
    ['#EC7063', '#E74C3C', '#B03A2E', '#78281F', '#4A110A'],
    ['#48C9B0', '#1ABC9C', '#148F77', '#0E6251', '#07362C'],
    ['#AF7AC5', '#884EA0', '#633974', '#4A235A', '#2C1338']
]

CURVE_PALETTE = [
    "#004488", "#E66100", "#117733", "#BB5566", "#5D3A9B",
    "#2874A6", "#B9770E", "#922B21", "#148F77", "#633974"
]

MAIN_COLORS = ['#004488', '#DDAA33', '#BB5566', '#117733', '#882255', '#44AA99']

# Canales estándar buscados en archivos .fcs
DEFAULT_PERM_CHANNELS = ('FITC', 'FL1', 'BL1', 'GFP', 'CALCEIN')
DEFAULT_FSC_CHANNELS = ('FSC-A', 'FSC-H', 'FSC')
DEFAULT_SSC_CHANNELS = ('SSC-A', 'SSC-H', 'SSC')
DEFAULT_TIME_CHANNELS = ('TIME', 'TIME-A', 'TIEMPO')

# Verificación de librerías de exportación de alta calidad
HAS_MPL = False
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    plt.rcParams.update({
        'font.sans-serif': 'DejaVu Sans',
        'font.family': 'sans-serif',
        'font.size': 10,
        'axes.labelsize': 11,
        'axes.titlesize': 11.5,
        'xtick.labelsize': 9.5,
        'ytick.labelsize': 9.5,
        'legend.fontsize': 8.5,
        'figure.titlesize': 12.5,
        'savefig.dpi': 600,
        'figure.dpi': 600
    })
    HAS_MPL = True
except Exception:
    HAS_MPL = False

HAS_REPORTLAB = False
try:
    import reportlab
    from reportlab.lib.pagesizes import letter
    from reportlab.lib import colors as rl_colors
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table as RLTable, TableStyle, Image as RLImage, PageBreak
    )
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    HAS_REPORTLAB = True
except Exception:
    HAS_REPORTLAB = False
