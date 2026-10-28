"""
Generador de archivos FCS sintéticos (FCS 3.0 estándar) para validación y pruebas automatizadas.
"""
from pathlib import Path
import numpy as np


def write_synthetic_fcs_file(filepath, time_s, fsc_arr, ssc_arr, fitc_arr):
    """Escribe un archivo binario .fcs estándar FCS 3.0 en disco."""
    n_events = len(time_s)
    num_par = 4
    timestep = 0.00002
    raw_time = (time_s / timestep).astype(np.float32)

    mat = np.column_stack([
        fsc_arr.astype(np.float32),
        ssc_arr.astype(np.float32),
        fitc_arr.astype(np.float32),
        raw_time
    ])
    data_bytes = mat.astype('>f4').tobytes()

    text_dict = {
        '$BEGINANALYSIS': '0', '$ENDANALYSIS': '0',
        '$DATATYPE': 'F', '$BYTEORD': '4,3,2,1', '$MODE': 'L',
        '$PAR': str(num_par), '$TOT': str(n_events), '$TIMESTEP': str(timestep),
        '$P1N': 'FSC-A', '$P1S': 'FSC-A', '$P1B': '32', '$P1E': '0,0', '$P1R': '1048576',
        '$P2N': 'SSC-A', '$P2S': 'SSC-A', '$P2B': '32', '$P2E': '0,0', '$P2R': '1048576',
        '$P3N': 'FITC-A', '$P3S': 'FITC-A', '$P3B': '32', '$P3E': '0,0', '$P3R': '1048576',
        '$P4N': 'TIME', '$P4S': 'TIME', '$P4B': '32', '$P4E': '0,0', '$P4R': '1048576',
    }

    delim = '/'
    body = delim + delim.join(f"{k}{delim}{v}" for k, v in text_dict.items()) + delim
    text_start = 256
    data_start = text_start + len(body) + 64
    text_dict['$DATASTART'] = str(data_start)
    text_dict['$DATAEND'] = str(data_start + len(data_bytes) - 1)
    text_bytes = (delim + delim.join(f"{k}{delim}{v}" for k, v in text_dict.items()) + delim).encode('latin1')
    text_end = text_start + len(text_bytes) - 1

    header = f"FCS3.0    {text_start:8d}{text_end:8d}{data_start:8d}{data_start + len(data_bytes) - 1:8d}       0       0"
    header_bytes = header.encode('ascii')[:58].ljust(text_start, b' ')

    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, 'wb') as f:
        f.write(header_bytes)
        f.write(text_bytes)
        f.write(b' ' * (data_start - (text_start + len(text_bytes))))
        f.write(data_bytes)


def generate_demo_fcs_dataset(base_folder):
    """
    Genera un conjunto de archivos .fcs sintéticos realistas
    (Satanin 1 y Act8-20 en THP-1 y PBMCs con controles basales CC).
    """
    demo_dir = Path(base_folder) / "demo_fcs_data"
    demo_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(12345)
    created_files = []

    configs = [
        ('sat1', [0.70, 1.40, 2.80, 5.70]),
        ('820', [0.50, 1.09, 2.18, 4.35])
    ]

    for folder_code, doses in configs:
        sub_dir = demo_dir / folder_code
        sub_dir.mkdir(exist_ok=True)
        for cell in ['THP-1', 'PBMCs']:
            for d in doses:
                for rep in [1, 2, 3]:
                    # 1. Archivo Control Basal CC (-60 a 0 s)
                    n_cc = 1200
                    t_cc = np.sort(rng.uniform(0, 60, n_cc))
                    fsc_cc = rng.normal(380000, 35000, n_cc)
                    ssc_cc = rng.normal(140000, 18000, n_cc)
                    fitc_cc = rng.normal(250000, 18000, n_cc)
                    cc_name = f"{cell}_CC_{d:.2f}_{rep}.fcs"
                    cc_path = sub_dir / cc_name
                    write_synthetic_fcs_file(cc_path, t_cc, fsc_cc, ssc_cc, fitc_cc)
                    created_files.append(str(cc_path))

                    # 2. Archivo Experimental (0 a 600 s)
                    n_exp = 4000
                    t_exp = np.sort(rng.uniform(0, 600, n_exp))
                    if cell == 'PBMCs':
                        perm_curve = np.clip(rng.normal(0.01, 0.015, n_exp), 0.0, 0.04)
                        swell = np.ones(n_exp)
                    else:
                        if d == doses[0]:
                            perm_curve = 0.10 / (1.0 + np.exp(-(t_exp - 420.0) / 45.0))
                            swell = 1.0 + 0.08 * perm_curve
                        elif d == doses[1]:
                            perm_curve = 0.18 / (1.0 + np.exp(-(t_exp - 510.0) / 35.0))
                            swell = 1.0 + 0.12 * perm_curve
                        elif d == doses[2]:
                            perm_curve = 0.77 / (1.0 + np.exp(-(t_exp - 235.0) / 28.0))
                            swell = 1.0 + 0.28 * perm_curve
                        else:
                            shock = 0.05 * np.exp(-((t_exp - 25.0) / 12.0) ** 2)
                            perm_curve = shock + 0.79 / (1.0 + np.exp(-(t_exp - 135.0) / 18.0))
                            swell = 1.0 + 0.36 * perm_curve

                    perm_curve = np.clip(perm_curve + rng.normal(0, 0.015, n_exp), 0.0, 0.96)
                    is_lysed_event = rng.random(n_exp) < perm_curve
                    fitc_intact = 250000.0 * (1.0 - 0.25 * perm_curve) * rng.normal(1.0, 0.12, n_exp)
                    fitc_lysed = 18000.0 * rng.normal(1.0, 0.22, n_exp)
                    fitc_exp = np.where(is_lysed_event, np.maximum(800.0, fitc_lysed), np.maximum(5000.0, fitc_intact))
                    fsc_exp = 380000.0 * swell * rng.normal(1.0, 0.06, n_exp)
                    ssc_exp = 140000.0 * (1.0 - 0.18 * perm_curve) * rng.normal(1.0, 0.06, n_exp)

                    exp_name = f"{cell}_{d:.2f}_{rep}.fcs"
                    exp_path = sub_dir / exp_name
                    write_synthetic_fcs_file(exp_path, t_exp, fsc_exp, ssc_exp, fitc_exp)
                    created_files.append(str(exp_path))

    return created_files
