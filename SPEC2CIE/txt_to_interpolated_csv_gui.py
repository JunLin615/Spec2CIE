#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FDTD reflectance TXT -> CIELUV-compatible CSV converter.

This module can be used in two ways:

1. As an importable library for other Python projects.
2. As a standalone Tkinter GUI application.

Expected input data
-------------------
The input text file should contain at least two numeric columns per data row:

    wavelength_in_meters, value

A header such as ``lambda(m), Y`` is allowed and is ignored automatically.
Comma, semicolon, tab, or whitespace delimiters are accepted.

Output CSV
----------
The output contains no header and exactly two columns:

    wavelength_nm, interpolated_value

The default wavelength grid is 380-800 nm at 1 nm intervals.
Linear interpolation is used and extrapolation is intentionally disabled.

Only the Python standard library is required.
"""

from __future__ import annotations

import bisect
import csv
import math
import re
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence


__version__ = "2.0.0"

APP_TITLE = "FDTD Spectrum to CIELUV CSV"
DEFAULT_START_NM = 380.0
DEFAULT_END_NM = 800.0
DEFAULT_STEP_NM = 1.0
DEFAULT_WAVELENGTH_SCALE = 1e9  # meters -> nanometers


# ---------------------------------------------------------------------------
# Reusable data structures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SpectrumData:
    """A wavelength/value spectrum with wavelengths stored in nanometers."""

    wavelength_nm: tuple[float, ...]
    values: tuple[float, ...]

    def __post_init__(self) -> None:
        if len(self.wavelength_nm) != len(self.values):
            raise ValueError("Wavelength and value arrays must have the same length.")
        if len(self.wavelength_nm) < 2:
            raise ValueError("At least two spectrum points are required.")

    @property
    def min_wavelength_nm(self) -> float:
        return self.wavelength_nm[0]

    @property
    def max_wavelength_nm(self) -> float:
        return self.wavelength_nm[-1]

    def __len__(self) -> int:
        return len(self.wavelength_nm)


@dataclass(frozen=True)
class ConversionResult:
    """Metadata returned after one successful file conversion."""

    input_path: Path
    output_path: Path
    point_count: int
    source_min_nm: float
    source_max_nm: float


@dataclass(frozen=True)
class ConversionFailure:
    """Information about one failed file conversion."""

    input_path: Path
    error: str


@dataclass(frozen=True)
class BatchConversionResult:
    """Summary returned by :func:`convert_batch`."""

    successes: tuple[ConversionResult, ...]
    failures: tuple[ConversionFailure, ...]

    @property
    def success_count(self) -> int:
        return len(self.successes)

    @property
    def failure_count(self) -> int:
        return len(self.failures)

    @property
    def total_count(self) -> int:
        return self.success_count + self.failure_count


# ---------------------------------------------------------------------------
# Parsing and spectrum preparation
# ---------------------------------------------------------------------------


_NUMBER_SPLIT_RE = re.compile(r"[,;\s]+")


def parse_numeric_pair(line: str) -> tuple[float, float] | None:
    """Parse the first two numeric fields from one text line.

    Returns ``None`` for blank lines, headers, comments, or malformed rows.
    Delimiters may be commas, semicolons, tabs, or whitespace.
    """

    parts = [part for part in _NUMBER_SPLIT_RE.split(line.strip()) if part]
    if len(parts) < 2:
        return None

    try:
        first = float(parts[0])
        second = float(parts[1])
    except ValueError:
        return None

    if not (math.isfinite(first) and math.isfinite(second)):
        return None

    return first, second


def normalize_spectrum_pairs(
    pairs: Iterable[tuple[float, float]],
) -> SpectrumData:
    """Sort wavelength/value pairs and average exact duplicate wavelengths."""

    sorted_pairs = sorted(pairs, key=lambda pair: pair[0])
    if len(sorted_pairs) < 2:
        raise ValueError("At least two valid numeric data points are required.")

    wavelengths: list[float] = []
    values: list[float] = []

    index = 0
    while index < len(sorted_pairs):
        wavelength = sorted_pairs[index][0]
        value_sum = sorted_pairs[index][1]
        count = 1
        index += 1

        while index < len(sorted_pairs) and sorted_pairs[index][0] == wavelength:
            value_sum += sorted_pairs[index][1]
            count += 1
            index += 1

        wavelengths.append(wavelength)
        values.append(value_sum / count)

    if len(wavelengths) < 2:
        raise ValueError("Only one unique wavelength is available; interpolation is impossible.")

    return SpectrumData(tuple(wavelengths), tuple(values))


def load_spectrum_txt(
    txt_path: str | Path,
    *,
    wavelength_scale: float = DEFAULT_WAVELENGTH_SCALE,
    encoding: str = "utf-8-sig",
) -> SpectrumData:
    """Load an FDTD-style spectrum text file.

    Parameters
    ----------
    txt_path:
        Input text file.
    wavelength_scale:
        Multiplier applied to the first numeric column. The default ``1e9``
        converts wavelengths from meters to nanometers.
    encoding:
        Text encoding used to read the file.

    Returns
    -------
    SpectrumData
        Sorted spectrum in nanometers with duplicate wavelengths averaged.
    """

    path = Path(txt_path)
    if not path.is_file():
        raise FileNotFoundError(f"Input file does not exist: {path}")
    if not math.isfinite(wavelength_scale) or wavelength_scale == 0:
        raise ValueError("wavelength_scale must be a finite, non-zero number.")

    pairs: list[tuple[float, float]] = []
    with path.open("r", encoding=encoding, errors="replace") as handle:
        for line in handle:
            parsed = parse_numeric_pair(line)
            if parsed is None:
                continue
            wavelength_raw, value = parsed
            pairs.append((wavelength_raw * wavelength_scale, value))

    return normalize_spectrum_pairs(pairs)


# ---------------------------------------------------------------------------
# Grid generation and interpolation
# ---------------------------------------------------------------------------


def validate_resampling_parameters(
    start_nm: float,
    end_nm: float,
    step_nm: float,
) -> None:
    """Validate wavelength-grid parameters."""

    if not all(math.isfinite(v) for v in (start_nm, end_nm, step_nm)):
        raise ValueError("Start, end, and step must be finite numbers.")
    if step_nm <= 0:
        raise ValueError("Step must be greater than 0 nm.")
    if end_nm < start_nm:
        raise ValueError("End wavelength must be greater than or equal to start wavelength.")


def build_wavelength_grid(
    start_nm: float = DEFAULT_START_NM,
    end_nm: float = DEFAULT_END_NM,
    step_nm: float = DEFAULT_STEP_NM,
) -> list[float]:
    """Create an evenly spaced wavelength grid in nanometers.

    The endpoint is included when it is exactly reachable within normal
    floating-point tolerance.
    """

    validate_resampling_parameters(start_nm, end_nm, step_nm)

    span = end_nm - start_nm
    count = int(math.floor(span / step_nm + 1e-12))
    grid = [start_nm + i * step_nm for i in range(count + 1)]

    tolerance = max(1e-10, abs(end_nm) * 1e-12)
    if grid and abs(grid[-1] - end_nm) <= tolerance:
        grid[-1] = end_nm

    return grid


def interpolate_value(
    wavelengths_nm: Sequence[float],
    values: Sequence[float],
    wavelength_nm: float,
) -> float:
    """Linearly interpolate one point from an ascending wavelength array.

    Extrapolation is not allowed.
    """

    if len(wavelengths_nm) != len(values) or len(wavelengths_nm) < 2:
        raise ValueError("Interpolation requires matching arrays with at least two points.")

    minimum = wavelengths_nm[0]
    maximum = wavelengths_nm[-1]
    tolerance = max(
        1e-10,
        max(abs(minimum), abs(maximum), abs(wavelength_nm)) * 1e-12,
    )

    if wavelength_nm < minimum - tolerance or wavelength_nm > maximum + tolerance:
        raise ValueError(
            f"Target wavelength {wavelength_nm:g} nm is outside the source range "
            f"[{minimum:g}, {maximum:g}] nm. Extrapolation is disabled."
        )

    if abs(wavelength_nm - minimum) <= tolerance:
        return values[0]
    if abs(wavelength_nm - maximum) <= tolerance:
        return values[-1]

    index = bisect.bisect_left(wavelengths_nm, wavelength_nm)
    if index < len(wavelengths_nm) and abs(wavelengths_nm[index] - wavelength_nm) <= tolerance:
        return values[index]

    x0, x1 = wavelengths_nm[index - 1], wavelengths_nm[index]
    y0, y1 = values[index - 1], values[index]

    fraction = (wavelength_nm - x0) / (x1 - x0)
    return y0 + fraction * (y1 - y0)


def resample_spectrum(
    spectrum: SpectrumData,
    *,
    start_nm: float = DEFAULT_START_NM,
    end_nm: float = DEFAULT_END_NM,
    step_nm: float = DEFAULT_STEP_NM,
) -> SpectrumData:
    """Resample a spectrum onto an evenly spaced wavelength grid."""

    grid = build_wavelength_grid(start_nm, end_nm, step_nm)
    if not grid:
        raise ValueError("No target wavelength points were generated.")

    source_min = spectrum.min_wavelength_nm
    source_max = spectrum.max_wavelength_nm
    tolerance = max(1e-10, max(abs(source_min), abs(source_max)) * 1e-12)

    if grid[0] < source_min - tolerance or grid[-1] > source_max + tolerance:
        raise ValueError(
            f"Requested range [{grid[0]:g}, {grid[-1]:g}] nm is outside the source range "
            f"[{source_min:g}, {source_max:g}] nm. Extrapolation is disabled."
        )

    interpolated = tuple(
        interpolate_value(spectrum.wavelength_nm, spectrum.values, wavelength)
        for wavelength in grid
    )
    return SpectrumData(tuple(grid), interpolated)


# Backward-compatible names from v1.
make_grid = build_wavelength_grid
linear_interpolate = interpolate_value


def read_txt_data(txt_path: str | Path) -> tuple[list[float], list[float]]:
    """Backward-compatible wrapper returning separate wavelength/value lists."""

    spectrum = load_spectrum_txt(txt_path)
    return list(spectrum.wavelength_nm), list(spectrum.values)


# ---------------------------------------------------------------------------
# CSV writing and file conversion
# ---------------------------------------------------------------------------


def _format_wavelength(value: float) -> str:
    if abs(value - round(value)) < 1e-10:
        return str(int(round(value)))
    return f"{value:.12g}"


def _format_value(value: float) -> str:
    return f"{value:.12g}"


def write_spectrum_csv(
    csv_path: str | Path,
    spectrum: SpectrumData,
    *,
    create_parent: bool = True,
) -> Path:
    """Write a headerless two-column CSV: wavelength_nm, value."""

    path = Path(csv_path)
    if create_parent:
        path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        for wavelength, value in zip(spectrum.wavelength_nm, spectrum.values):
            writer.writerow([_format_wavelength(wavelength), _format_value(value)])

    return path


def convert_fdtd_txt_to_csv(
    txt_path: str | Path,
    csv_path: str | Path,
    *,
    start_nm: float = DEFAULT_START_NM,
    end_nm: float = DEFAULT_END_NM,
    step_nm: float = DEFAULT_STEP_NM,
    wavelength_scale: float = DEFAULT_WAVELENGTH_SCALE,
) -> ConversionResult:
    """Convert one FDTD-exported TXT spectrum to the target CSV format."""

    input_path = Path(txt_path)
    output_path = Path(csv_path)

    source = load_spectrum_txt(input_path, wavelength_scale=wavelength_scale)
    resampled = resample_spectrum(
        source,
        start_nm=start_nm,
        end_nm=end_nm,
        step_nm=step_nm,
    )
    write_spectrum_csv(output_path, resampled)

    return ConversionResult(
        input_path=input_path,
        output_path=output_path,
        point_count=len(resampled),
        source_min_nm=source.min_wavelength_nm,
        source_max_nm=source.max_wavelength_nm,
    )


def convert_one_file(
    txt_path: str | Path,
    csv_path: str | Path,
    start_nm: float,
    end_nm: float,
    step_nm: float,
) -> tuple[int, float, float]:
    """Backward-compatible wrapper matching the v1 function signature."""

    result = convert_fdtd_txt_to_csv(
        txt_path,
        csv_path,
        start_nm=start_nm,
        end_nm=end_nm,
        step_nm=step_nm,
    )
    return result.point_count, result.source_min_nm, result.source_max_nm


def find_txt_files(folder: str | Path, *, recursive: bool = True) -> list[Path]:
    """Return sorted TXT files from a folder."""

    root = Path(folder)
    if not root.is_dir():
        raise NotADirectoryError(f"Input folder does not exist: {root}")

    iterator = root.rglob("*") if recursive else root.glob("*")
    return sorted(
        path
        for path in iterator
        if path.is_file() and path.suffix.lower() == ".txt"
    )


def output_path_for_input(
    txt_path: str | Path,
    *,
    input_root: str | Path,
    output_root: str | Path,
) -> Path:
    """Map an input TXT path to its CSV path while preserving subfolders."""

    txt = Path(txt_path)
    input_base = Path(input_root)
    output_base = Path(output_root)
    relative = txt.relative_to(input_base)
    return output_base / relative.with_suffix(".csv")


def convert_batch(
    input_folder: str | Path,
    output_folder: str | Path,
    *,
    start_nm: float = DEFAULT_START_NM,
    end_nm: float = DEFAULT_END_NM,
    step_nm: float = DEFAULT_STEP_NM,
    wavelength_scale: float = DEFAULT_WAVELENGTH_SCALE,
    recursive: bool = True,
    progress_callback: Callable[[int, int, Path, ConversionResult | ConversionFailure], None] | None = None,
) -> BatchConversionResult:
    """Convert every TXT file in a folder while preserving subfolder structure.

    The optional ``progress_callback`` receives
    ``(current_index, total_count, input_path, result_or_failure)``.
    """

    input_root = Path(input_folder)
    output_root = Path(output_folder)
    files = find_txt_files(input_root, recursive=recursive)

    successes: list[ConversionResult] = []
    failures: list[ConversionFailure] = []
    total = len(files)

    for index, txt_path in enumerate(files, start=1):
        csv_path = output_path_for_input(
            txt_path,
            input_root=input_root,
            output_root=output_root,
        )

        try:
            result: ConversionResult | ConversionFailure = convert_fdtd_txt_to_csv(
                txt_path,
                csv_path,
                start_nm=start_nm,
                end_nm=end_nm,
                step_nm=step_nm,
                wavelength_scale=wavelength_scale,
            )
            successes.append(result)
        except Exception as exc:  # batch processing should continue after one bad file
            result = ConversionFailure(txt_path, str(exc))
            failures.append(result)

        if progress_callback is not None:
            progress_callback(index, total, txt_path, result)

    return BatchConversionResult(tuple(successes), tuple(failures))


# ---------------------------------------------------------------------------
# Optional Tkinter GUI
# ---------------------------------------------------------------------------


def launch_gui() -> None:
    """Launch the standalone Tkinter user interface.

    Tkinter is imported lazily so importing this module for data processing does
    not create GUI dependencies at import time.
    """

    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    class ConverterGUI(tk.Tk):
        def __init__(self) -> None:
            super().__init__()
            self.title(APP_TITLE)
            self.geometry("820x560")
            self.minsize(740, 510)

            self.input_mode: str | None = None
            self.input_path = tk.StringVar()
            self.output_path = tk.StringVar()
            self.start_nm = tk.StringVar(value=f"{DEFAULT_START_NM:g}")
            self.end_nm = tk.StringVar(value=f"{DEFAULT_END_NM:g}")
            self.step_nm = tk.StringVar(value=f"{DEFAULT_STEP_NM:g}")
            self.status_var = tk.StringVar(value="Select a TXT file or a folder to begin.")

            self._build_ui()

        def _build_ui(self) -> None:
            outer = ttk.Frame(self, padding=14)
            outer.pack(fill="both", expand=True)
            outer.columnconfigure(1, weight=1)

            ttk.Label(outer, text="Input").grid(row=0, column=0, sticky="w", pady=(0, 6))
            ttk.Entry(outer, textvariable=self.input_path, state="readonly").grid(
                row=0, column=1, columnspan=2, sticky="ew", padx=(8, 8), pady=(0, 6)
            )
            ttk.Button(outer, text="Select TXT File", command=self.choose_file).grid(
                row=0, column=3, padx=(0, 6), pady=(0, 6)
            )
            ttk.Button(outer, text="Select Folder", command=self.choose_folder).grid(
                row=0, column=4, pady=(0, 6)
            )

            ttk.Label(outer, text="Output Folder").grid(row=1, column=0, sticky="w", pady=6)
            ttk.Entry(outer, textvariable=self.output_path, state="readonly").grid(
                row=1, column=1, columnspan=3, sticky="ew", padx=(8, 8), pady=6
            )
            ttk.Button(outer, text="Select Output", command=self.choose_output).grid(
                row=1, column=4, pady=6
            )

            params = ttk.LabelFrame(outer, text="Interpolation Settings (nm)", padding=10)
            params.grid(row=2, column=0, columnspan=5, sticky="ew", pady=(12, 8))
            for column in range(6):
                params.columnconfigure(column, weight=1 if column in (1, 3, 5) else 0)

            ttk.Label(params, text="Start").grid(row=0, column=0, sticky="e", padx=(0, 6))
            ttk.Entry(params, textvariable=self.start_nm, width=12).grid(
                row=0, column=1, sticky="ew"
            )
            ttk.Label(params, text="End").grid(row=0, column=2, sticky="e", padx=(16, 6))
            ttk.Entry(params, textvariable=self.end_nm, width=12).grid(
                row=0, column=3, sticky="ew"
            )
            ttk.Label(params, text="Step").grid(row=0, column=4, sticky="e", padx=(16, 6))
            ttk.Entry(params, textvariable=self.step_nm, width=12).grid(
                row=0, column=5, sticky="ew"
            )

            note = (
                "Output CSV: no header, two columns (wavelength in nm, interpolated value).\n"
                "Folder mode scans TXT files recursively and preserves the original subfolder structure."
            )
            ttk.Label(outer, text=note, justify="left").grid(
                row=3, column=0, columnspan=5, sticky="w", pady=(4, 10)
            )

            self.convert_button = ttk.Button(
                outer,
                text="Convert",
                command=self.start_conversion,
            )
            self.convert_button.grid(
                row=4,
                column=0,
                columnspan=5,
                sticky="ew",
                pady=(4, 10),
            )

            self.progress = ttk.Progressbar(outer, mode="determinate", maximum=100)
            self.progress.grid(row=5, column=0, columnspan=5, sticky="ew", pady=(0, 8))

            ttk.Label(outer, textvariable=self.status_var).grid(
                row=6, column=0, columnspan=5, sticky="w", pady=(0, 6)
            )

            log_frame = ttk.LabelFrame(outer, text="Conversion Log", padding=6)
            log_frame.grid(row=7, column=0, columnspan=5, sticky="nsew")
            outer.rowconfigure(7, weight=1)
            log_frame.rowconfigure(0, weight=1)
            log_frame.columnconfigure(0, weight=1)

            self.log = tk.Text(log_frame, height=12, wrap="word", state="disabled")
            self.log.grid(row=0, column=0, sticky="nsew")
            scrollbar = ttk.Scrollbar(log_frame, orient="vertical", command=self.log.yview)
            scrollbar.grid(row=0, column=1, sticky="ns")
            self.log.configure(yscrollcommand=scrollbar.set)

        def choose_file(self) -> None:
            path = filedialog.askopenfilename(
                title="Select FDTD TXT File",
                filetypes=[("TXT files", "*.txt"), ("All files", "*.*")],
            )
            if path:
                self.input_mode = "file"
                self.input_path.set(path)
                if not self.output_path.get():
                    self.output_path.set(str(Path(path).parent))
                self.status_var.set("Single-file mode selected.")

        def choose_folder(self) -> None:
            path = filedialog.askdirectory(title="Select Folder Containing TXT Files")
            if path:
                self.input_mode = "folder"
                self.input_path.set(path)
                if not self.output_path.get():
                    source = Path(path)
                    self.output_path.set(str(source.parent / f"{source.name}_csv"))
                self.status_var.set("Folder mode selected; TXT files will be scanned recursively.")

        def choose_output(self) -> None:
            path = filedialog.askdirectory(title="Select CSV Output Folder")
            if path:
                self.output_path.set(path)

        def _append_log(self, text: str) -> None:
            def update() -> None:
                self.log.configure(state="normal")
                self.log.insert("end", text + "\n")
                self.log.see("end")
                self.log.configure(state="disabled")

            self.after(0, update)

        def _set_progress(self, current: int, total: int) -> None:
            value = 0 if total <= 0 else current / total * 100
            self.after(0, lambda: self.progress.configure(value=value))
            self.after(0, lambda: self.status_var.set(f"Converting: {current}/{total}"))

        def _read_parameters(self) -> tuple[float, float, float]:
            start_nm = float(self.start_nm.get())
            end_nm = float(self.end_nm.get())
            step_nm = float(self.step_nm.get())
            validate_resampling_parameters(start_nm, end_nm, step_nm)
            return start_nm, end_nm, step_nm

        def start_conversion(self) -> None:
            if not self.input_mode or not self.input_path.get():
                messagebox.showerror("Missing Input", "Select a TXT file or folder first.")
                return
            if not self.output_path.get():
                messagebox.showerror("Missing Output", "Select an output folder first.")
                return

            try:
                start_nm, end_nm, step_nm = self._read_parameters()
            except ValueError as exc:
                messagebox.showerror("Invalid Parameters", str(exc))
                return

            input_path = Path(self.input_path.get())
            output_root = Path(self.output_path.get())

            if self.input_mode == "file":
                files = [input_path]
            else:
                try:
                    files = find_txt_files(input_path, recursive=True)
                except Exception as exc:
                    messagebox.showerror("Input Error", str(exc))
                    return

            if not files:
                messagebox.showwarning("No Files Found", "No .txt files were found in the selected location.")
                return

            self.convert_button.configure(state="disabled")
            self.progress.configure(value=0)
            self.status_var.set(f"Preparing to convert {len(files)} file(s)...")
            self._append_log("=" * 70)
            self._append_log(
                f"Settings: start={start_nm:g} nm, end={end_nm:g} nm, step={step_nm:g} nm"
            )

            worker = threading.Thread(
                target=self._conversion_worker,
                args=(files, input_path, output_root, start_nm, end_nm, step_nm),
                daemon=True,
            )
            worker.start()

        def _conversion_worker(
            self,
            files: Sequence[Path],
            input_path: Path,
            output_root: Path,
            start_nm: float,
            end_nm: float,
            step_nm: float,
        ) -> None:
            successes: list[ConversionResult] = []
            failures: list[ConversionFailure] = []

            for index, txt_path in enumerate(files, start=1):
                if self.input_mode == "file":
                    csv_path = output_root / f"{txt_path.stem}.csv"
                else:
                    csv_path = output_path_for_input(
                        txt_path,
                        input_root=input_path,
                        output_root=output_root,
                    )

                try:
                    result = convert_fdtd_txt_to_csv(
                        txt_path,
                        csv_path,
                        start_nm=start_nm,
                        end_nm=end_nm,
                        step_nm=step_nm,
                    )
                    successes.append(result)
                    self._append_log(
                        f"OK   {txt_path} -> {csv_path} "
                        f"({result.point_count} points; source "
                        f"{result.source_min_nm:g}-{result.source_max_nm:g} nm)"
                    )
                except Exception as exc:
                    failure = ConversionFailure(txt_path, str(exc))
                    failures.append(failure)
                    self._append_log(f"FAIL {txt_path}: {exc}")

                self._set_progress(index, len(files))

            def finish() -> None:
                self.convert_button.configure(state="normal")
                success_count = len(successes)
                failure_count = len(failures)

                if failure_count == 0:
                    self.status_var.set(f"Completed: {success_count} file(s) converted successfully.")
                    messagebox.showinfo(
                        "Conversion Complete",
                        f"Successfully converted {success_count} file(s).",
                    )
                else:
                    self.status_var.set(
                        f"Completed: {success_count} succeeded, {failure_count} failed."
                    )
                    preview = "\n".join(
                        f"{failure.input_path.name}: {failure.error}"
                        for failure in failures[:8]
                    )
                    if failure_count > 8:
                        preview += f"\n... and {failure_count - 8} more error(s). See the log."
                    messagebox.showwarning(
                        "Conversion Complete with Errors",
                        f"Succeeded: {success_count}\nFailed: {failure_count}\n\n{preview}",
                    )

            self.after(0, finish)

    app = ConverterGUI()
    app.mainloop()


if __name__ == "__main__":
    launch_gui()
