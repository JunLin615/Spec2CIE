"""Matplotlib plotting utilities for SPEC2CIE.

The plotting layer is independent from Tkinter. It consumes ``AnalyzedSpectrum``
objects from :mod:`SPEC2CIE.colorimetry` and returns ordinary Matplotlib figures,
so the same functions can be used by the desktop GUI, notebooks, scripts, or
future batch pipelines.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from .colorimetry import (
    AnalyzedSpectrum,
    OBSERVER_CIE1931_2,
    readable_text_rgb,
)


MARKER_SEQUENCE: tuple[str, ...] = (
    "o", "s", "^", "v", "<", ">", "D", "d", "p", "h", "H", "P", "X", "*",
    "8", "+", "x", "1", "2", "3", "4", "|", "_", ".", ",", "$A$", "$B$",
    "$C$", "$D$", "$E$", "$F$", "$G$", "$H$", "$I$", "$J$", "$K$", "$L$",
    "$M$", "$N$", "$O$", "$P$", "$Q$", "$R$", "$S$", "$T$", "$U$", "$V$",
    "$W$", "$X$", "$Y$", "$Z$",
)


@dataclass(frozen=True)
class PlotStyle:
    """Common appearance settings for CIE, zoom, and spectrum figures."""

    font_family: str = "Arial"
    axis_font_size: float = 12.0
    tick_font_size: float = 10.0
    legend_font_size: float = 10.0
    marker_size: float = 8.0
    axis_line_width: float = 1.2
    spectral_locus_line_width: float = 1.0
    figure_width_in: float = 6.5
    figure_height_in: float = 6.0
    dpi: int = 300
    show_grid: bool = False
    show_spectral_locus: bool = True
    show_diagram_colours: bool = True
    show_point_labels: bool = False
    show_legend: bool = True
    show_zoom_rectangle: bool = True
    match_legend_text_color: bool = True
    improve_light_text_readability: bool = True
    # Lock the physical Axes rectangle rather than the full Figure canvas.
    # Main CIE uses its visible data-range ratio to preserve equal u'/v' scale.
    lock_main_cie_box_aspect: bool = False
    lock_spectra_box_aspect: bool = False
    spectra_box_aspect: float = 0.75  # Axes height / width.


@dataclass(frozen=True)
class CIERange:
    """Visible range of a CIE 1976 u'v' diagram.

    The defaults are intentionally tighter than the colour-science plotting
    defaults while still covering the full spectral locus with a small margin.
    """

    u_min: float = 0.00
    u_max: float = 0.63
    v_min: float = 0.00
    v_max: float = 0.60

    def normalized(self) -> "CIERange":
        return CIERange(
            min(self.u_min, self.u_max),
            max(self.u_min, self.u_max),
            min(self.v_min, self.v_max),
            max(self.v_min, self.v_max),
        )

    @property
    def width(self) -> float:
        region = self.normalized()
        return region.u_max - region.u_min

    @property
    def height(self) -> float:
        region = self.normalized()
        return region.v_max - region.v_min


@dataclass(frozen=True)
class ZoomRegion:
    u_min: float = 0.12
    u_max: float = 0.28
    v_min: float = 0.35
    v_max: float = 0.60

    def normalized(self) -> "ZoomRegion":
        return ZoomRegion(
            min(self.u_min, self.u_max),
            max(self.u_min, self.u_max),
            min(self.v_min, self.v_max),
            max(self.v_min, self.v_max),
        )

    @property
    def width(self) -> float:
        region = self.normalized()
        return region.u_max - region.u_min

    @property
    def height(self) -> float:
        region = self.normalized()
        return region.v_max - region.v_min

    def contains(self, u: float, v: float) -> bool:
        region = self.normalized()
        return region.u_min <= u <= region.u_max and region.v_min <= v <= region.v_max


@dataclass(frozen=True)
class PointAppearance:
    marker: str
    label: str


def marker_for_index(index: int) -> str:
    return MARKER_SEQUENCE[index % len(MARKER_SEQUENCE)]


def _rc_params(style: PlotStyle) -> dict:
    return {
        "font.family": style.font_family,
        "font.size": style.tick_font_size,
        "axes.labelsize": style.axis_font_size,
        "axes.titlesize": style.axis_font_size,
        "xtick.labelsize": style.tick_font_size,
        "ytick.labelsize": style.tick_font_size,
        "legend.fontsize": style.legend_font_size,
        "axes.linewidth": style.axis_line_width,
        "savefig.dpi": style.dpi,
    }


def _import_colour_plotting():
    try:
        from colour.plotting import plot_chromaticity_diagram_CIE1976UCS  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "The 'colour-science' package is required for CIE diagram plotting."
        ) from exc
    return plot_chromaticity_diagram_CIE1976UCS


def _style_axes(ax, style: PlotStyle, *, x_label: str | None = None, y_label: str | None = None) -> None:
    if x_label is not None:
        ax.set_xlabel(x_label, fontsize=style.axis_font_size)
    if y_label is not None:
        ax.set_ylabel(y_label, fontsize=style.axis_font_size)
    ax.tick_params(width=style.axis_line_width, labelsize=style.tick_font_size)
    for spine in ax.spines.values():
        spine.set_linewidth(style.axis_line_width)
    ax.grid(style.show_grid, linewidth=0.5, alpha=0.35)


def _new_figure(style: PlotStyle) -> tuple[Figure, object]:
    fig, ax = plt.subplots(
        figsize=(style.figure_width_in, style.figure_height_in),
        dpi=max(72, min(style.dpi, 200)),
    )
    return fig, ax


def _draw_cie_background(ax, style: PlotStyle, observer: str, bounding_box=None) -> None:
    plot_cie = _import_colour_plotting()
    kwargs = {
        "cmfs": observer,
        "show_diagram_colours": style.show_diagram_colours,
        "show_spectral_locus": style.show_spectral_locus,
        "axes": ax,
        "show": False,
        "title": None,
        "x_label": "u'",
        "y_label": "v'",
    }
    if bounding_box is not None:
        kwargs["bounding_box"] = bounding_box
    plot_cie(**kwargs)

    # Colour's default spectral locus style is intentionally overridden only in
    # width, preserving its wavelength-colour treatment.
    if style.show_spectral_locus:
        for line in ax.lines:
            line.set_linewidth(style.spectral_locus_line_width)


def cie1976_position_rgb(u_prime: float, v_prime: float) -> tuple[float, float, float]:
    """Return the same display RGB mapping used by Colour's CIE 1976 background.

    The colour is derived only from the CIE 1976 chromaticity position, not from
    the measured Y tristimulus value and not from the selected illuminant white.
    Consequently, a marker placed on the diagram visually matches the diagram
    colour at that position. When a new observer or illuminant is analysed, the
    point coordinates move and the plot colour follows the new coordinates.

    This deliberately mirrors ``colour.plotting``: CIE 1976 u'v' -> XYZ at a
    nominal luminance -> plotting sRGB -> per-colour maximum normalisation.
    """

    try:
        import numpy as np
        from colour.algebra import normalise_maximum  # type: ignore
        from colour.plotting import (  # type: ignore
            CONSTANTS_COLOUR_STYLE,
            XYZ_to_plotting_colourspace,
        )
        from colour.plotting.diagrams import METHODS_CHROMATICITY_DIAGRAM  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "The 'colour-science' package is required for CIE plot colours."
        ) from exc

    ij = np.asarray([float(u_prime), float(v_prime)], dtype=float)
    if not np.all(np.isfinite(ij)):
        return (0.5, 0.5, 0.5)

    illuminant = CONSTANTS_COLOUR_STYLE.colour.colourspace.whitepoint
    ij_to_XYZ = METHODS_CHROMATICITY_DIAGRAM["CIE 1976 UCS"]["ij_to_XYZ"]
    XYZ = np.asarray(ij_to_XYZ(ij, illuminant), dtype=float)
    RGB = np.asarray(XYZ_to_plotting_colourspace(XYZ, illuminant), dtype=float)
    RGB = np.asarray(normalise_maximum(RGB, axis=-1), dtype=float)
    RGB = np.nan_to_num(RGB, nan=0.0, posinf=1.0, neginf=0.0)
    RGB = np.clip(RGB, 0.0, 1.0)
    return float(RGB[0]), float(RGB[1]), float(RGB[2])


def _plot_colour_map(
    analyses: Sequence[AnalyzedSpectrum],
) -> dict[str, tuple[float, float, float]]:
    """Return stable CIE-position colours for all current analyses."""

    return {
        str(a.result.source_path): cie1976_position_rgb(
            a.result.u_prime, a.result.v_prime
        )
        for a in analyses
    }


def _appearance_map(
    analyses: Sequence[AnalyzedSpectrum],
    appearances: Mapping[str, PointAppearance] | None,
) -> dict[str, PointAppearance]:
    mapping: dict[str, PointAppearance] = {}
    for index, analysis in enumerate(analyses):
        key = str(analysis.result.source_path)
        if appearances is not None and key in appearances:
            mapping[key] = appearances[key]
        else:
            mapping[key] = PointAppearance(marker_for_index(index), analysis.result.label)
    return mapping


def _plot_points(
    ax,
    analyses: Sequence[AnalyzedSpectrum],
    style: PlotStyle,
    appearances: Mapping[str, PointAppearance] | None = None,
) -> None:
    mapping = _appearance_map(analyses, appearances)
    plot_colours = _plot_colour_map(analyses)
    legend_handles = []
    legend_labels = []
    legend_text_colours = []

    for analysis in analyses:
        result = analysis.result
        appearance = mapping[str(result.source_path)]
        marker = appearance.marker
        colour = plot_colours[str(result.source_path)]

        handle = ax.plot(
            result.u_prime,
            result.v_prime,
            linestyle="None",
            marker=marker,
            markersize=style.marker_size,
            markerfacecolor=colour,
            markeredgecolor="black",
            markeredgewidth=0.6,
            zorder=20,
        )[0]

        if style.show_point_labels:
            ax.annotate(
                appearance.label,
                (result.u_prime, result.v_prime),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=style.tick_font_size,
                color=readable_text_rgb(colour),
                zorder=21,
            )

        legend_handles.append(handle)
        legend_labels.append(appearance.label)
        if style.improve_light_text_readability:
            legend_text_colours.append(readable_text_rgb(colour))
        else:
            legend_text_colours.append(colour)

    if style.show_legend and legend_handles:
        legend = ax.legend(
            legend_handles,
            legend_labels,
            loc="best",
            frameon=True,
            fontsize=style.legend_font_size,
        )
        if style.match_legend_text_color:
            for text, colour in zip(legend.get_texts(), legend_text_colours):
                text.set_color(colour)


def plot_cie1976_main(
    analyses: Iterable[AnalyzedSpectrum],
    *,
    style: PlotStyle = PlotStyle(),
    observer: str = OBSERVER_CIE1931_2,
    main_range: CIERange | None = CIERange(),
    zoom_region: ZoomRegion | None = None,
    appearances: Mapping[str, PointAppearance] | None = None,
) -> Figure:
    """Create a publication-ready main CIE 1976 u'v' figure.

    ``main_range`` controls only the visible range of the main CIE diagram.
    Pass ``None`` to fall back to colour-science's native plotting bounds.
    """

    items = list(analyses)
    region = main_range.normalized() if main_range is not None else None
    if region is not None and (region.width <= 0 or region.height <= 0):
        raise ValueError("Main CIE range must have non-zero width and height.")

    with plt.rc_context(_rc_params(style)):
        fig, ax = _new_figure(style)
        bounding_box = None
        if region is not None:
            bounding_box = (region.u_min, region.u_max, region.v_min, region.v_max)
        _draw_cie_background(ax, style, observer, bounding_box=bounding_box)
        if region is not None:
            ax.set_xlim(region.u_min, region.u_max)
            ax.set_ylim(region.v_min, region.v_max)
        _plot_points(ax, items, style, appearances)
        _style_axes(ax, style, x_label="u'", y_label="v'")

        if style.lock_main_cie_box_aspect:
            if region is not None:
                # Height/width of the physical axes rectangle. This value also
                # preserves equal physical scale for one u' and one v' unit.
                ax.set_box_aspect(region.height / region.width)
            else:
                # Native Colour bounds are normally square-ish; use the actual
                # resolved limits so locking never changes the plotted limits.
                x0, x1 = ax.get_xlim()
                y0, y1 = ax.get_ylim()
                if x1 != x0 and y1 != y0:
                    ax.set_box_aspect(abs((y1 - y0) / (x1 - x0)))

        if zoom_region is not None and style.show_zoom_rectangle:
            region = zoom_region.normalized()
            if region.width > 0 and region.height > 0:
                ax.add_patch(
                    Rectangle(
                        (region.u_min, region.v_min),
                        region.width,
                        region.height,
                        fill=False,
                        linewidth=max(1.0, style.axis_line_width),
                        linestyle="--",
                        edgecolor="black",
                        zorder=25,
                    )
                )

        fig.tight_layout()
        return fig


def plot_cie1976_zoom(
    analyses: Iterable[AnalyzedSpectrum],
    region: ZoomRegion,
    *,
    style: PlotStyle = PlotStyle(),
    observer: str = OBSERVER_CIE1931_2,
    appearances: Mapping[str, PointAppearance] | None = None,
    include_only_points_inside: bool = True,
) -> Figure:
    """Create a zoomed CIE 1976 figure for a selected u'v' region."""

    region = region.normalized()
    if region.width <= 0 or region.height <= 0:
        raise ValueError("Zoom region must have non-zero width and height.")

    items = list(analyses)
    if include_only_points_inside:
        points = [a for a in items if region.contains(a.result.u_prime, a.result.v_prime)]
    else:
        points = items

    with plt.rc_context(_rc_params(style)):
        fig, ax = _new_figure(style)
        _draw_cie_background(
            ax,
            style,
            observer,
            bounding_box=(region.u_min, region.u_max, region.v_min, region.v_max),
        )
        ax.set_xlim(region.u_min, region.u_max)
        ax.set_ylim(region.v_min, region.v_max)
        _plot_points(ax, points, style, appearances)
        _style_axes(ax, style, x_label="u'", y_label="v'")
        fig.tight_layout()
        return fig


def plot_spectra(
    analyses: Iterable[AnalyzedSpectrum],
    *,
    style: PlotStyle = PlotStyle(),
    mode: str = "Raw reflectance",
    appearances: Mapping[str, PointAppearance] | None = None,
) -> Figure:
    """Plot raw reflectance or illumination-weighted relative spectra."""

    items = list(analyses)
    mapping = _appearance_map(items, appearances)
    plot_colours = _plot_colour_map(items)
    with plt.rc_context(_rc_params(style)):
        fig, ax = _new_figure(style)

        for analysis in items:
            appearance = mapping[str(analysis.result.source_path)]
            spectrum = (
                analysis.illuminated_relative
                if mode == "Illuminated spectrum S(lambda) * R(lambda)"
                else analysis.reflectance
            )
            ax.plot(
                spectrum.wavelength_nm,
                spectrum.values,
                label=appearance.label,
                linewidth=1.5,
                color=plot_colours[str(analysis.result.source_path)],
            )

        ax.set_xlabel("Wavelength (nm)", fontsize=style.axis_font_size)
        if mode == "Illuminated spectrum S(lambda) * R(lambda)":
            ax.set_ylabel("Relative illuminated reflectance", fontsize=style.axis_font_size)
        else:
            ax.set_ylabel("Reflectance", fontsize=style.axis_font_size)
        _style_axes(ax, style)
        if style.show_legend and items:
            legend = ax.legend(loc="best", fontsize=style.legend_font_size)
            if style.match_legend_text_color:
                for text, analysis in zip(legend.get_texts(), items):
                    colour = plot_colours[str(analysis.result.source_path)]
                    if style.improve_light_text_readability:
                        colour = readable_text_rgb(colour)
                    text.set_color(colour)
        if style.lock_spectra_box_aspect:
            ax.set_box_aspect(style.spectra_box_aspect)
        fig.tight_layout()
        return fig


def export_figure(figure: Figure, path: str, *, dpi: int = 300, transparent_background: bool | None = None) -> None:
    """Export a Matplotlib figure with robust PNG transparency.

    For transparent export, the Figure patch, every Axes patch, and legend frame
    are temporarily forced to alpha=0 before rendering. Their original appearance
    is restored immediately afterwards, so the on-screen GUI is unchanged.

    PNG files are transparent by default. A post-export alpha check is performed
    for PNG so an unexpectedly opaque file is reported instead of silently saved.
    """

    suffix = Path(path).suffix.lower()
    if transparent_background is None:
        transparent_background = suffix == ".png"

    # Save patch state so transparent export does not alter the interactive figure.
    patch_states = []

    def make_patch_transparent(patch) -> None:
        if patch is None:
            return
        patch_states.append((patch, patch.get_facecolor(), patch.get_edgecolor(), patch.get_alpha()))
        patch.set_facecolor("none")
        patch.set_edgecolor("none")
        patch.set_alpha(0.0)

    if transparent_background:
        make_patch_transparent(figure.patch)
        for ax in figure.axes:
            make_patch_transparent(ax.patch)
            legend = ax.get_legend()
            if legend is not None:
                make_patch_transparent(legend.get_frame())

    save_kwargs = {
        "dpi": dpi,
        "bbox_inches": "tight",
        "transparent": bool(transparent_background),
    }
    if transparent_background:
        save_kwargs["facecolor"] = (0.0, 0.0, 0.0, 0.0)
        save_kwargs["edgecolor"] = (0.0, 0.0, 0.0, 0.0)

    try:
        figure.savefig(path, **save_kwargs)
    finally:
        for patch, facecolor, edgecolor, alpha in reversed(patch_states):
            patch.set_facecolor(facecolor)
            patch.set_edgecolor(edgecolor)
            patch.set_alpha(alpha)

        # Refresh an interactive canvas after restoring its visible background.
        canvas = getattr(figure, "canvas", None)
        if canvas is not None:
            try:
                canvas.draw_idle()
            except Exception:
                pass

    # Validate the actual written PNG, not just Matplotlib's requested settings.
    if transparent_background and suffix == ".png":
        try:
            from PIL import Image

            with Image.open(path) as image:
                if "A" not in image.getbands():
                    raise RuntimeError(
                        "PNG export produced no alpha channel even though transparent export was requested."
                    )
                alpha_min, alpha_max = image.getchannel("A").getextrema()
                if alpha_min == 255 and alpha_max == 255:
                    raise RuntimeError(
                        "PNG export produced an alpha channel, but every pixel is fully opaque. "
                        "The transparent background was not applied by the active Matplotlib backend."
                    )
        except ImportError:
            # Pillow is normally present because Matplotlib depends on it. If it is
            # unavailable, the file is still kept; only validation is skipped.
            pass


__all__ = [
    "MARKER_SEQUENCE",
    "CIERange",
    "PlotStyle",
    "PointAppearance",
    "ZoomRegion",
    "cie1976_position_rgb",
    "export_figure",
    "marker_for_index",
    "plot_cie1976_main",
    "plot_cie1976_zoom",
    "plot_spectra",
]
