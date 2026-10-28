"""
Núcleo biofísico, bioestadístico y modelos matemáticos de citometría cinética.
"""
from .statistics import (
    wilson_binomial_ci,
    cohen_h_proportions,
    contingency_2x2_p_and_or,
    cochran_armitage_trend_test,
    interpolate_crossing_time,
    classify_affectation_status,
    adjust_pvalues
)
from .models import (
    hill_model,
    michaelis_menten_model,
    gompertz_dose_model,
    weibull_dose_model,
    biexponential_time_model,
    fit_dose_response_model
)
from .kinetics import (
    extract_tube_kinetics,
    apply_gate_and_cutoffs_to_cache,
    extract_slope_based_parameters,
    compute_categorical_kinetics_metrics
)
from .synthetic_demo import (
    write_synthetic_fcs_file,
    generate_demo_fcs_dataset
)

__all__ = [
    'wilson_binomial_ci',
    'cohen_h_proportions',
    'contingency_2x2_p_and_or',
    'cochran_armitage_trend_test',
    'interpolate_crossing_time',
    'classify_affectation_status',
    'adjust_pvalues',
    'hill_model',
    'michaelis_menten_model',
    'gompertz_dose_model',
    'weibull_dose_model',
    'biexponential_time_model',
    'fit_dose_response_model',
    'extract_tube_kinetics',
    'apply_gate_and_cutoffs_to_cache',
    'extract_slope_based_parameters',
    'compute_categorical_kinetics_metrics',
    'write_synthetic_fcs_file',
    'generate_demo_fcs_dataset'
]
