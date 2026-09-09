"""Colour-science engine for FDTD reflectance spectra.

The public API in this module is intentionally independent from the Tkinter GUI.
It converts FDTD reflectance spectra into CIE XYZ, CIE 1931 xy, and CIE 1976
u'v' coordinates under a selectable illuminant.

The FDTD reader and linear resampler are imported from
``txt_to_interpolated_csv_gui`` so the same parsing behaviour is shared by the
CSV converter and the CIE application.
"""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from .txt_to_interpolated_csv_gui import (
    SpectrumData,
    load_spectrum_txt,
    parse_numeric_pair,
    resample_spectrum,
)


CIE_START_NM = 360.0
CIE_END_NM = 800.0
CIE_STEP_NM = 1.0

OBSERVER_CIE1931_2 = "CIE 1931 2 Degree Standard Observer"
OBSERVER_CIE1964_10 = "CIE 1964 10 Degree Standard Observer"

OBSERVER_OPTIONS: tuple[str, ...] = (
    OBSERVER_CIE1931_2,
    OBSERVER_CIE1964_10,
)

ILLUMINANT_EQUAL_ENERGY = "Equal-energy E"
ILLUMINANT_D65 = "CIE D65 - Daylight"
ILLUMINANT_D50 = "CIE D50 - Daylight / 5000 K"
ILLUMINANT_A = "CIE A - Tungsten / Halogen"
ILLUMINANT_BB3000 = "Warm white approximation - Blackbody 3000 K"
ILLUMINANT_BB4000 = "Neutral white approximation - Blackbody 4000 K"
ILLUMINANT_BB6500 = "Cool white approximation - Blackbody 6500 K"
ILLUMINANT_CUSTOM = "Custom illuminant spectrum"

ILLUMINANT_OPTIONS: tuple[str, ...] = (
    ILLUMINANT_EQUAL_ENERGY,
    ILLUMINANT_D65,
    ILLUMINANT_D50,
    ILLUMINANT_A,
    ILLUMINANT_BB3000,
    ILLUMINANT_BB4000,
    ILLUMINANT_BB6500,
    ILLUMINANT_CUSTOM,
)

_STANDARD_ILLUMINANT_KEYS = {
    ILLUMINANT_EQUAL_ENERGY: "E",
    ILLUMINANT_D65: "D65",
    ILLUMINANT_D50: "D50",
    ILLUMINANT_A: "A",
}

_BLACKBODY_TEMPERATURES = {
    ILLUMINANT_BB3000: 3000.0,
    ILLUMINANT_BB4000: 4000.0,
    ILLUMINANT_BB6500: 6500.0,
}


@dataclass(frozen=True)
class AnalysisSettings:
    """Settings used for colourimetric analysis."""

    observer: str = OBSERVER_CIE1931_2
    illuminant: str = ILLUMINANT_EQUAL_ENERGY
    custom_illuminant_path: Path | None = None
    integration_start_nm: float = CIE_START_NM
    integration_end_nm: float = CIE_END_NM
    integration_step_nm: float = CIE_STEP_NM


@dataclass(frozen=True)
class ColorimetricResult:
    """Colourimetric result for one reflectance spectrum."""

    source_path: Path
    label: str
    observer: str
    illuminant: str
    X: float
    Y: float
    Z: float
    x: float
    y: float
    u_prime: float
    v_prime: float
    rgb: tuple[float, float, float]
    hex_color: str
    out_of_srgb_gamut: bool
    reflectance_min: float
    reflectance_max: float
    source_min_nm: float
    source_max_nm: float
    warning: str = ""

    @property
    def xyz(self) -> tuple[float, float, float]:
        return self.X, self.Y, self.Z

    @property
    def xy(self) -> tuple[float, float]:
        return self.x, self.y

    @property
    def uv(self) -> tuple[float, float]:
        return self.u_prime, self.v_prime


@dataclass(frozen=True)
class AnalyzedSpectrum:
    """A result plus the spectra used by the plotting layer."""

    result: ColorimetricResult
    reflectance: SpectrumData
    illuminated_relative: SpectrumData


class ColourDependencyError(RuntimeError):
    """Raised when the optional scientific colour package is unavailable."""


def _import_colour():
    try:
        import colour  # type: ignore
    except ImportError as exc:  # pragma: no cover - depends on local environment
        raise ColourDependencyError(
            "The 'colour-science' package is required for CIE calculations. "
            "Install the project dependencies with 'pip install -e .'."
        ) from exc
    return colour


def xyz_to_xy(XYZ: Sequence[float]) -> tuple[float, float]:
    """Convert CIE XYZ to CIE 1931 x,y chromaticity coordinates."""

    X, Y, Z = (float(v) for v in XYZ)
    denominator = X + Y + Z
    if not math.isfinite(denominator) or abs(denominator) < 1e-15:
        return float("nan"), float("nan")
    return X / denominator, Y / denominator


def xyz_to_uv1976(XYZ: Sequence[float]) -> tuple[float, float]:
    """Convert CIE XYZ to CIE 1976 UCS u',v' coordinates."""

    X, Y, Z = (float(v) for v in XYZ)
    denominator = X + 15.0 * Y + 3.0 * Z
    if not math.isfinite(denominator) or abs(denominator) < 1e-15:
        return float("nan"), float("nan")
    return 4.0 * X / denominator, 9.0 * Y / denominator


def _format_hex(rgb: Sequence[float]) -> str:
    channels = [int(round(max(0.0, min(1.0, float(c))) * 255.0)) for c in rgb]
    return "#{:02X}{:02X}{:02X}".format(*channels)


def readable_text_rgb(rgb: Sequence[float], max_luminance: float = 0.55) -> tuple[float, float, float]:
    """Darken a bright display RGB colour for readable legend text on white."""

    r, g, b = (max(0.0, min(1.0, float(v))) for v in rgb)
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    if luminance <= max_luminance or luminance <= 0:
        return r, g, b
    factor = max_luminance / luminance
    return r * factor, g * factor, b * factor


def _spectrum_to_colour_sd(spectrum: SpectrumData, name: str):
    colour = _import_colour()
    data = dict(zip(spectrum.wavelength_nm, spectrum.values))
    return colour.SpectralDistribution(data, name=name)


def _spectral_shape(settings: AnalysisSettings):
    colour = _import_colour()
    return colour.SpectralShape(
        settings.integration_start_nm,
        settings.integration_end_nm,
        settings.integration_step_nm,
    )


def _load_custom_two_column_spectrum(path: str | Path) -> SpectrumData:
    """Load a custom illuminant spectrum with automatic m/nm wavelength units.

    The first two numeric fields on each row are used. If all wavelengths are
    smaller than 0.01, the first column is assumed to be in metres and is
    multiplied by 1e9; otherwise it is treated as nanometres.
    """

    file_path = Path(path)
    if not file_path.is_file():
        raise FileNotFoundError(f"Custom illuminant file does not exist: {file_path}")

    raw_pairs: list[tuple[float, float]] = []
    with file_path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for line in handle:
            parsed = parse_numeric_pair(line)
            if parsed is not None:
                raw_pairs.append(parsed)

    if len(raw_pairs) < 2:
        raise ValueError("Custom illuminant file must contain at least two numeric data rows.")

    raw_wavelengths = [pair[0] for pair in raw_pairs]
    scale = 1e9 if max(abs(v) for v in raw_wavelengths) < 0.01 else 1.0
    pairs = sorted((w * scale, value) for w, value in raw_pairs)

    wavelengths: list[float] = []
    values: list[float] = []
    index = 0
    while index < len(pairs):
        wavelength = pairs[index][0]
        total = pairs[index][1]
        count = 1
        index += 1
        while index < len(pairs) and pairs[index][0] == wavelength:
            total += pairs[index][1]
            count += 1
            index += 1
        wavelengths.append(wavelength)
        values.append(total / count)

    return SpectrumData(tuple(wavelengths), tuple(values))


def get_illuminant_sd(settings: AnalysisSettings):
    """Return a Colour ``SpectralDistribution`` for the selected illuminant."""

    colour = _import_colour()
    shape = _spectral_shape(settings)

    if settings.illuminant in _STANDARD_ILLUMINANT_KEYS:
        key = _STANDARD_ILLUMINANT_KEYS[settings.illuminant]
        return colour.SDS_ILLUMINANTS[key].copy().align(shape)

    if settings.illuminant in _BLACKBODY_TEMPERATURES:
        temperature = _BLACKBODY_TEMPERATURES[settings.illuminant]
        sd = colour.sd_blackbody(temperature, shape=shape)
        # Absolute scaling is irrelevant here because reflected-object XYZ is
        # normalised to the selected illuminant white during colorimetry.
        sd.name = f"Blackbody {temperature:g} K"
        return sd

    if settings.illuminant == ILLUMINANT_CUSTOM:
        if settings.custom_illuminant_path is None:
            raise ValueError("Select a custom illuminant spectrum file first.")
        custom = _load_custom_two_column_spectrum(settings.custom_illuminant_path)
        custom = resample_spectrum(
            custom,
            start_nm=settings.integration_start_nm,
            end_nm=settings.integration_end_nm,
            step_nm=settings.integration_step_nm,
        )
        sd = _spectrum_to_colour_sd(custom, settings.custom_illuminant_path.stem)
        return sd.align(shape)

    raise ValueError(f"Unknown illuminant: {settings.illuminant}")


def get_observer_cmfs(observer: str):
    """Return the selected CIE standard observer colour matching functions."""

    if observer not in OBSERVER_OPTIONS:
        raise ValueError(f"Unknown observer: {observer}")
    colour = _import_colour()
    return colour.MSDS_CMFS[observer].copy()


def _white_xyz_for_illuminant(cmfs, illuminant, shape) -> np.ndarray:
    colour = _import_colour()
    perfect_reflector = colour.sd_ones(shape)
    XYZ_w = np.asarray(colour.sd_to_XYZ(perfect_reflector, cmfs, illuminant, method="Integration"), dtype=float)
    if XYZ_w[1] <= 0:
        raise ValueError("Selected illuminant produced an invalid reference white.")
    return XYZ_w / XYZ_w[1]


def xyz_to_display_rgb(XYZ: Sequence[float], cmfs, illuminant, shape) -> tuple[tuple[float, float, float], bool]:
    """Approximate an object colour as sRGB for plotting markers and legends.

    The selected illuminant white is chromatically adapted to the sRGB D65
    display white using the Bradford transform. This is a display aid, not a
    camera simulation.
    """

    colour = _import_colour()
    XYZ_array = np.asarray(XYZ, dtype=float) / 100.0
    white_XYZ = _white_xyz_for_illuminant(cmfs, illuminant, shape)
    source_xy = np.asarray(colour.XYZ_to_xy(white_XYZ), dtype=float)

    raw_rgb = np.asarray(
        colour.XYZ_to_sRGB(
            XYZ_array,
            illuminant=source_xy,
            chromatic_adaptation_transform="Bradford",
            apply_cctf_encoding=True,
        ),
        dtype=float,
    )
    out_of_gamut = bool(np.any(raw_rgb < -1e-9) or np.any(raw_rgb > 1.0 + 1e-9))
    clipped = np.clip(raw_rgb, 0.0, 1.0)
    return (float(clipped[0]), float(clipped[1]), float(clipped[2])), out_of_gamut


def analyze_reflectance_spectrum(
    spectrum: SpectrumData,
    *,
    source_path: str | Path,
    label: str | None = None,
    settings: AnalysisSettings = AnalysisSettings(),
) -> AnalyzedSpectrum:
    """Calculate XYZ, xy, u'v', and an approximate display colour."""

    colour = _import_colour()
    shape = _spectral_shape(settings)
    cmfs = get_observer_cmfs(settings.observer).align(shape)
    illuminant = get_illuminant_sd(settings).align(shape)

    resampled = resample_spectrum(
        spectrum,
        start_nm=settings.integration_start_nm,
        end_nm=settings.integration_end_nm,
        step_nm=settings.integration_step_nm,
    )

    reflectance_values = np.asarray(resampled.values, dtype=float)
    reflectance_sd = _spectrum_to_colour_sd(resampled, label or Path(source_path).stem).align(shape)
    XYZ = np.asarray(colour.sd_to_XYZ(reflectance_sd, cmfs, illuminant, method="Integration"), dtype=float)

    x, y = xyz_to_xy(XYZ)
    u_prime, v_prime = xyz_to_uv1976(XYZ)
    rgb, out_of_gamut = xyz_to_display_rgb(XYZ, cmfs, illuminant, shape)

    illumination_values = np.asarray(illuminant.values, dtype=float)
    illumination_max = float(np.nanmax(illumination_values))
    if illumination_max > 0:
        illumination_relative = illumination_values / illumination_max
    else:
        illumination_relative = illumination_values
    illuminated_values = reflectance_values * illumination_relative
    illuminated = SpectrumData(resampled.wavelength_nm, tuple(float(v) for v in illuminated_values))

    reflectance_min = float(np.nanmin(reflectance_values))
    reflectance_max = float(np.nanmax(reflectance_values))
    warnings: list[str] = []
    if reflectance_min < -1e-6 or reflectance_max > 1.0 + 1e-6:
        warnings.append(
            f"reflectance outside [0, 1] ({reflectance_min:.4g} to {reflectance_max:.4g})"
        )
    if out_of_gamut:
        warnings.append("display colour is outside sRGB gamut and was clipped")

    path = Path(source_path)
    result = ColorimetricResult(
        source_path=path,
        label=label or path.stem,
        observer=settings.observer,
        illuminant=settings.illuminant,
        X=float(XYZ[0]),
        Y=float(XYZ[1]),
        Z=float(XYZ[2]),
        x=float(x),
        y=float(y),
        u_prime=float(u_prime),
        v_prime=float(v_prime),
        rgb=rgb,
        hex_color=_format_hex(rgb),
        out_of_srgb_gamut=out_of_gamut,
        reflectance_min=reflectance_min,
        reflectance_max=reflectance_max,
        source_min_nm=spectrum.min_wavelength_nm,
        source_max_nm=spectrum.max_wavelength_nm,
        warning="; ".join(warnings),
    )
    return AnalyzedSpectrum(result=result, reflectance=resampled, illuminated_relative=illuminated)


def analyze_fdtd_file(
    txt_path: str | Path,
    *,
    label: str | None = None,
    settings: AnalysisSettings = AnalysisSettings(),
) -> AnalyzedSpectrum:
    """Load one FDTD TXT file and calculate its CIE coordinates."""

    path = Path(txt_path)
    source = load_spectrum_txt(path)
    return analyze_reflectance_spectrum(
        source,
        source_path=path,
        label=label,
        settings=settings,
    )


def analyze_fdtd_files(
    txt_paths: Iterable[str | Path],
    *,
    settings: AnalysisSettings = AnalysisSettings(),
) -> list[AnalyzedSpectrum]:
    """Analyze multiple FDTD spectra in input order."""

    return [analyze_fdtd_file(path, settings=settings) for path in txt_paths]


def export_results_csv(path: str | Path, analyses: Iterable[AnalyzedSpectrum]) -> Path:
    """Export a results table with XYZ, xy, u'v', and approximate display RGB."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "file",
                "label",
                "observer",
                "illuminant",
                "X",
                "Y",
                "Z",
                "x",
                "y",
                "u_prime",
                "v_prime",
                "sRGB_R",
                "sRGB_G",
                "sRGB_B",
                "HEX",
                "out_of_sRGB_gamut",
                "reflectance_min",
                "reflectance_max",
                "warning",
            ]
        )
        for analysis in analyses:
            r = analysis.result
            writer.writerow(
                [
                    str(r.source_path),
                    r.label,
                    r.observer,
                    r.illuminant,
                    f"{r.X:.12g}",
                    f"{r.Y:.12g}",
                    f"{r.Z:.12g}",
                    f"{r.x:.12g}",
                    f"{r.y:.12g}",
                    f"{r.u_prime:.12g}",
                    f"{r.v_prime:.12g}",
                    f"{r.rgb[0]:.12g}",
                    f"{r.rgb[1]:.12g}",
                    f"{r.rgb[2]:.12g}",
                    r.hex_color,
                    int(r.out_of_srgb_gamut),
                    f"{r.reflectance_min:.12g}",
                    f"{r.reflectance_max:.12g}",
                    r.warning,
                ]
            )
    return output


__all__ = [
    "AnalysisSettings",
    "AnalyzedSpectrum",
    "ColorimetricResult",
    "ColourDependencyError",
    "CIE_START_NM",
    "CIE_END_NM",
    "CIE_STEP_NM",
    "OBSERVER_CIE1931_2",
    "OBSERVER_CIE1964_10",
    "OBSERVER_OPTIONS",
    "ILLUMINANT_EQUAL_ENERGY",
    "ILLUMINANT_D65",
    "ILLUMINANT_D50",
    "ILLUMINANT_A",
    "ILLUMINANT_BB3000",
    "ILLUMINANT_BB4000",
    "ILLUMINANT_BB6500",
    "ILLUMINANT_CUSTOM",
    "ILLUMINANT_OPTIONS",
    "analyze_fdtd_file",
    "analyze_fdtd_files",
    "analyze_reflectance_spectrum",
    "export_results_csv",
    "get_illuminant_sd",
    "get_observer_cmfs",
    "readable_text_rgb",
    "xyz_to_display_rgb",
    "xyz_to_uv1976",
    "xyz_to_xy",
]
