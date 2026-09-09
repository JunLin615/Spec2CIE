"""SPEC2CIE: FDTD reflectance spectra to CIE 1976 u'v' analysis."""

from .colorimetry import (
    AnalysisSettings,
    AnalyzedSpectrum,
    ColorimetricResult,
    ILLUMINANT_OPTIONS,
    OBSERVER_OPTIONS,
    analyze_fdtd_file,
    analyze_fdtd_files,
    analyze_reflectance_spectrum,
    export_results_csv,
    xyz_to_uv1976,
    xyz_to_xy,
)
from .plotting import (
    MARKER_SEQUENCE,
    PlotStyle,
    PointAppearance,
    ZoomRegion,
    plot_cie1976_main,
    plot_cie1976_zoom,
    plot_spectra,
)
from .txt_to_interpolated_csv_gui import (
    SpectrumData,
    convert_batch,
    convert_fdtd_txt_to_csv,
    load_spectrum_txt,
    resample_spectrum,
    write_spectrum_csv,
)

__version__ = "1.0.0"

__all__ = [
    "AnalysisSettings",
    "AnalyzedSpectrum",
    "ColorimetricResult",
    "ILLUMINANT_OPTIONS",
    "MARKER_SEQUENCE",
    "OBSERVER_OPTIONS",
    "PlotStyle",
    "PointAppearance",
    "SpectrumData",
    "ZoomRegion",
    "analyze_fdtd_file",
    "analyze_fdtd_files",
    "analyze_reflectance_spectrum",
    "convert_batch",
    "convert_fdtd_txt_to_csv",
    "export_results_csv",
    "load_spectrum_txt",
    "plot_cie1976_main",
    "plot_cie1976_zoom",
    "plot_spectra",
    "resample_spectrum",
    "write_spectrum_csv",
    "xyz_to_uv1976",
    "xyz_to_xy",
]
