"""Tkinter desktop application for SPEC2CIE.

Run as a module::

    python -m SPEC2CIE.cie1976_fdtd_viewer

or install the project and use the ``spec2cie`` console command.
"""

from __future__ import annotations

import threading
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

from .colorimetry import (
    AnalysisSettings,
    AnalyzedSpectrum,
    CIE_END_NM,
    CIE_START_NM,
    CIE_STEP_NM,
    ILLUMINANT_CUSTOM,
    ILLUMINANT_EQUAL_ENERGY,
    ILLUMINANT_OPTIONS,
    OBSERVER_CIE1931_2,
    OBSERVER_OPTIONS,
    analyze_fdtd_file,
    export_results_csv,
)
try:
    # v2 filename is used for side-by-side testing. After the user renames the
    # files back to their canonical names, the fallback import below is used.
    from .plotting import (
        CIERange,
        PlotStyle,
        PointAppearance,
        ZoomRegion,
        export_figure,
        marker_for_index,
        plot_cie1976_main,
        plot_cie1976_zoom,
        plot_spectra,
    )
except ImportError:
    from .plotting import (
        CIERange,
        PlotStyle,
        PointAppearance,
        ZoomRegion,
        export_figure,
        marker_for_index,
        plot_cie1976_main,
        plot_cie1976_zoom,
        plot_spectra,
    )
from .txt_to_interpolated_csv_gui import find_txt_files


APP_TITLE = "SPEC2CIE - FDTD Reflectance to CIE 1976"
__version__ = "1.1.0"


@dataclass
class SpectrumEntry:
    path: Path
    label: str
    marker: str
    analysis: AnalyzedSpectrum | None = None
    error: str = ""


class SPEC2CIEApp:
    """Desktop GUI controller. Tkinter is imported lazily by ``launch_gui``."""

    def __init__(self, root) -> None:
        import tkinter as tk
        from tkinter import ttk

        self.tk = tk
        self.ttk = ttk
        self.root = root
        self.entries: list[SpectrumEntry] = []
        self._analysis_running = False
        self._zoom_selector = None

        self.main_figure = None
        self.zoom_figure = None
        self.spectra_figure = None
        self.main_canvas = None
        self.zoom_canvas = None
        self.spectra_canvas = None
        self.main_toolbar = None
        self.zoom_toolbar = None
        self.spectra_toolbar = None

        root.title(APP_TITLE)
        root.geometry("1500x900")
        root.minsize(1180, 720)

        # Colourimetry controls.
        self.observer_var = tk.StringVar(value=OBSERVER_CIE1931_2)
        self.illuminant_var = tk.StringVar(value=ILLUMINANT_EQUAL_ENERGY)
        self.custom_illuminant_var = tk.StringVar(value="")

        # Plot controls.
        self.font_family_var = tk.StringVar(value="Arial")
        self.axis_font_size_var = tk.StringVar(value="12")
        self.tick_font_size_var = tk.StringVar(value="10")
        self.legend_font_size_var = tk.StringVar(value="10")
        self.marker_size_var = tk.StringVar(value="8")
        self.axis_line_width_var = tk.StringVar(value="1.2")
        self.figure_width_var = tk.StringVar(value="6.5")
        self.figure_height_var = tk.StringVar(value="6.0")
        self.dpi_var = tk.StringVar(value="300")

        self.show_grid_var = tk.BooleanVar(value=False)
        self.show_locus_var = tk.BooleanVar(value=True)
        self.show_background_var = tk.BooleanVar(value=True)
        self.show_labels_var = tk.BooleanVar(value=False)
        self.show_legend_var = tk.BooleanVar(value=True)
        self.show_zoom_rectangle_var = tk.BooleanVar(value=True)
        self.legend_text_color_var = tk.BooleanVar(value=True)
        self.readable_legend_var = tk.BooleanVar(value=True)

        # Main CIE visible-range controls.
        self.use_default_main_range_var = tk.BooleanVar(value=True)
        self.main_u_min_var = tk.StringVar(value="0.00")
        self.main_u_max_var = tk.StringVar(value="0.63")
        self.main_v_min_var = tk.StringVar(value="0.00")
        self.main_v_max_var = tk.StringVar(value="0.60")
        self._main_range_entries = []

        # Zoom controls.
        self.zoom_u_min_var = tk.StringVar(value="0.12")
        self.zoom_u_max_var = tk.StringVar(value="0.28")
        self.zoom_v_min_var = tk.StringVar(value="0.35")
        self.zoom_v_max_var = tk.StringVar(value="0.60")

        self.spectra_mode_var = tk.StringVar(value="Raw reflectance")
        self.status_var = tk.StringVar(value="Add one or more FDTD TXT spectra to begin.")

        self._build_ui()

    # ------------------------------------------------------------------ UI

    def _build_ui(self) -> None:
        ttk = self.ttk

        outer = ttk.Frame(self.root, padding=8)
        outer.pack(fill="both", expand=True)
        outer.rowconfigure(0, weight=1)
        outer.columnconfigure(0, weight=1)

        paned = ttk.Panedwindow(outer, orient="horizontal")
        paned.grid(row=0, column=0, sticky="nsew")

        left = ttk.Frame(paned, padding=(4, 4, 8, 4), width=390)
        right = ttk.Frame(paned, padding=(8, 4, 4, 4))
        paned.add(left, weight=0)
        paned.add(right, weight=1)
        left.grid_propagate(False)

        left.columnconfigure(0, weight=1)
        right.columnconfigure(0, weight=1)
        right.rowconfigure(1, weight=1)

        self._build_input_panel(left)
        self._build_colorimetry_panel(left)
        self._build_plot_settings_panel(left)
        self._build_main_range_panel(left)
        self._build_zoom_panel(left)
        self._build_action_panel(left)

        self._build_right_toolbar(right)
        self._build_plot_notebook(right)

        status = ttk.Label(outer, textvariable=self.status_var, anchor="w")
        status.grid(row=1, column=0, sticky="ew", pady=(6, 0))

    def _build_input_panel(self, parent) -> None:
        ttk = self.ttk
        frame = ttk.LabelFrame(parent, text="Input Spectra", padding=8)
        frame.grid(row=0, column=0, sticky="nsew", pady=(0, 8))
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)

        columns = ("label", "marker", "u", "v")
        self.file_tree = ttk.Treeview(
            frame,
            columns=columns,
            show="headings",
            height=10,
            selectmode="extended",
        )
        self.file_tree.heading("label", text="Name")
        self.file_tree.heading("marker", text="Marker")
        self.file_tree.heading("u", text="u'")
        self.file_tree.heading("v", text="v'")
        self.file_tree.column("label", width=175, anchor="w")
        self.file_tree.column("marker", width=55, anchor="center")
        self.file_tree.column("u", width=65, anchor="e")
        self.file_tree.column("v", width=65, anchor="e")
        self.file_tree.grid(row=0, column=0, columnspan=6, sticky="nsew")

        scrollbar = ttk.Scrollbar(frame, orient="vertical", command=self.file_tree.yview)
        scrollbar.grid(row=0, column=6, sticky="ns")
        self.file_tree.configure(yscrollcommand=scrollbar.set)

        ttk.Button(frame, text="Add Files", command=self.add_files).grid(row=1, column=0, sticky="ew", pady=(6, 0), padx=(0, 3))
        ttk.Button(frame, text="Add Folder", command=self.add_folder).grid(row=1, column=1, sticky="ew", pady=(6, 0), padx=3)
        ttk.Button(frame, text="Remove", command=self.remove_selected).grid(row=1, column=2, sticky="ew", pady=(6, 0), padx=3)
        ttk.Button(frame, text="Rename", command=self.rename_selected).grid(row=1, column=3, sticky="ew", pady=(6, 0), padx=3)
        ttk.Button(frame, text="Up", command=lambda: self.move_selected(-1)).grid(row=1, column=4, sticky="ew", pady=(6, 0), padx=3)
        ttk.Button(frame, text="Down", command=lambda: self.move_selected(1)).grid(row=1, column=5, sticky="ew", pady=(6, 0), padx=(3, 0))

        ttk.Button(frame, text="Clear All", command=self.clear_all).grid(row=2, column=0, columnspan=6, sticky="ew", pady=(5, 0))

    def _build_colorimetry_panel(self, parent) -> None:
        ttk = self.ttk
        frame = ttk.LabelFrame(parent, text="Colorimetry", padding=8)
        frame.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        frame.columnconfigure(1, weight=1)

        ttk.Label(frame, text="Observer").grid(row=0, column=0, sticky="w")
        observer = ttk.Combobox(frame, textvariable=self.observer_var, values=OBSERVER_OPTIONS, state="readonly", width=31)
        observer.grid(row=0, column=1, sticky="ew", padx=(8, 0))
        observer.bind("<<ComboboxSelected>>", lambda _event: self._mark_analysis_stale())

        ttk.Label(frame, text="Illuminant").grid(row=1, column=0, sticky="w", pady=(6, 0))
        illuminant = ttk.Combobox(frame, textvariable=self.illuminant_var, values=ILLUMINANT_OPTIONS, state="readonly", width=31)
        illuminant.grid(row=1, column=1, sticky="ew", padx=(8, 0), pady=(6, 0))
        illuminant.bind("<<ComboboxSelected>>", self._on_illuminant_changed)

        ttk.Label(frame, text="Custom SPD").grid(row=2, column=0, sticky="w", pady=(6, 0))
        custom_row = ttk.Frame(frame)
        custom_row.grid(row=2, column=1, sticky="ew", padx=(8, 0), pady=(6, 0))
        custom_row.columnconfigure(0, weight=1)
        ttk.Entry(custom_row, textvariable=self.custom_illuminant_var, state="readonly").grid(row=0, column=0, sticky="ew")
        ttk.Button(custom_row, text="Browse", command=self.choose_custom_illuminant).grid(row=0, column=1, padx=(5, 0))

        ttk.Label(
            frame,
            text=f"Integration: {CIE_START_NM:g}-{CIE_END_NM:g} nm, {CIE_STEP_NM:g} nm interval",
        ).grid(row=3, column=0, columnspan=2, sticky="w", pady=(7, 0))

    def _build_plot_settings_panel(self, parent) -> None:
        ttk = self.ttk
        frame = ttk.LabelFrame(parent, text="Plot Settings", padding=8)
        frame.grid(row=2, column=0, sticky="ew", pady=(0, 8))
        for col in (1, 3):
            frame.columnconfigure(col, weight=1)

        fields = [
            ("Font", self.font_family_var, "Axis font", self.axis_font_size_var),
            ("Tick font", self.tick_font_size_var, "Legend font", self.legend_font_size_var),
            ("Marker size", self.marker_size_var, "Axis width", self.axis_line_width_var),
            ("Figure width", self.figure_width_var, "Figure height", self.figure_height_var),
            ("DPI", self.dpi_var, None, None),
        ]
        for row, (label1, var1, label2, var2) in enumerate(fields):
            ttk.Label(frame, text=label1).grid(row=row, column=0, sticky="w", pady=2)
            ttk.Entry(frame, textvariable=var1, width=10).grid(row=row, column=1, sticky="ew", padx=(5, 10), pady=2)
            if label2 is not None:
                ttk.Label(frame, text=label2).grid(row=row, column=2, sticky="w", pady=2)
                ttk.Entry(frame, textvariable=var2, width=10).grid(row=row, column=3, sticky="ew", padx=(5, 0), pady=2)

        checks = ttk.Frame(frame)
        checks.grid(row=5, column=0, columnspan=4, sticky="ew", pady=(5, 0))
        ttk.Checkbutton(checks, text="Diagram colours", variable=self.show_background_var).grid(row=0, column=0, sticky="w")
        ttk.Checkbutton(checks, text="Spectral locus", variable=self.show_locus_var).grid(row=0, column=1, sticky="w", padx=(8, 0))
        ttk.Checkbutton(checks, text="Grid", variable=self.show_grid_var).grid(row=1, column=0, sticky="w")
        ttk.Checkbutton(checks, text="Point labels", variable=self.show_labels_var).grid(row=1, column=1, sticky="w", padx=(8, 0))
        ttk.Checkbutton(checks, text="Legend", variable=self.show_legend_var).grid(row=2, column=0, sticky="w")
        ttk.Checkbutton(checks, text="Zoom rectangle", variable=self.show_zoom_rectangle_var).grid(row=2, column=1, sticky="w", padx=(8, 0))
        ttk.Checkbutton(checks, text="Legend text = point colour", variable=self.legend_text_color_var).grid(row=3, column=0, columnspan=2, sticky="w")
        ttk.Checkbutton(checks, text="Darken light legend text", variable=self.readable_legend_var).grid(row=4, column=0, columnspan=2, sticky="w")

    def _build_main_range_panel(self, parent) -> None:
        ttk = self.ttk
        frame = ttk.LabelFrame(parent, text="Main CIE Range", padding=8)
        frame.grid(row=3, column=0, sticky="ew", pady=(0, 8))
        for col in (1, 3):
            frame.columnconfigure(col, weight=1)

        ttk.Checkbutton(
            frame,
            text="Use default full-CIE range",
            variable=self.use_default_main_range_var,
            command=self._toggle_main_range_entries,
        ).grid(row=0, column=0, columnspan=4, sticky="w", pady=(0, 5))

        specs = (
            ("u' min", self.main_u_min_var, 1, 0),
            ("u' max", self.main_u_max_var, 1, 2),
            ("v' min", self.main_v_min_var, 2, 0),
            ("v' max", self.main_v_max_var, 2, 2),
        )
        for label, variable, row, column in specs:
            ttk.Label(frame, text=label).grid(row=row, column=column, sticky="w", pady=2)
            entry = ttk.Entry(frame, textvariable=variable, width=9)
            entry.grid(row=row, column=column + 1, sticky="ew", padx=(5, 10 if column == 0 else 0), pady=2)
            self._main_range_entries.append(entry)

        ttk.Label(
            frame,
            text="Default: u' 0.00-0.63, v' 0.00-0.60",
        ).grid(row=3, column=0, columnspan=4, sticky="w", pady=(4, 0))
        self._toggle_main_range_entries()

    def _toggle_main_range_entries(self) -> None:
        state = "disabled" if self.use_default_main_range_var.get() else "normal"
        for entry in self._main_range_entries:
            entry.configure(state=state)

    def _build_zoom_panel(self, parent) -> None:
        ttk = self.ttk
        frame = ttk.LabelFrame(parent, text="Zoom Region", padding=8)
        frame.grid(row=4, column=0, sticky="ew", pady=(0, 8))
        for col in (1, 3):
            frame.columnconfigure(col, weight=1)

        ttk.Label(frame, text="u' min").grid(row=0, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.zoom_u_min_var, width=9).grid(row=0, column=1, sticky="ew", padx=(5, 10))
        ttk.Label(frame, text="u' max").grid(row=0, column=2, sticky="w")
        ttk.Entry(frame, textvariable=self.zoom_u_max_var, width=9).grid(row=0, column=3, sticky="ew", padx=(5, 0))
        ttk.Label(frame, text="v' min").grid(row=1, column=0, sticky="w", pady=(4, 0))
        ttk.Entry(frame, textvariable=self.zoom_v_min_var, width=9).grid(row=1, column=1, sticky="ew", padx=(5, 10), pady=(4, 0))
        ttk.Label(frame, text="v' max").grid(row=1, column=2, sticky="w", pady=(4, 0))
        ttk.Entry(frame, textvariable=self.zoom_v_max_var, width=9).grid(row=1, column=3, sticky="ew", padx=(5, 0), pady=(4, 0))

        ttk.Button(frame, text="Select Region on Main Plot", command=self.activate_zoom_selector).grid(row=2, column=0, columnspan=4, sticky="ew", pady=(6, 0))

    def _build_action_panel(self, parent) -> None:
        ttk = self.ttk
        frame = ttk.Frame(parent)
        frame.grid(row=5, column=0, sticky="ew")
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)

        self.analyze_button = ttk.Button(frame, text="Analyze / Recalculate", command=self.start_analysis)
        self.analyze_button.grid(row=0, column=0, sticky="ew", padx=(0, 4))
        ttk.Button(frame, text="Update Plot Style", command=self.update_plots).grid(row=0, column=1, sticky="ew", padx=(4, 0))

    def _build_right_toolbar(self, parent) -> None:
        ttk = self.ttk
        bar = ttk.Frame(parent)
        bar.grid(row=0, column=0, sticky="ew", pady=(0, 6))

        ttk.Button(bar, text="Export Main", command=lambda: self.export_plot("main")).pack(side="left")
        ttk.Button(bar, text="Export Zoom", command=lambda: self.export_plot("zoom")).pack(side="left", padx=(5, 0))
        ttk.Button(bar, text="Export Spectra", command=lambda: self.export_plot("spectra")).pack(side="left", padx=(5, 0))
        ttk.Button(bar, text="Export Results CSV", command=self.export_results).pack(side="left", padx=(5, 0))

        ttk.Label(bar, text="Spectra view:").pack(side="right", padx=(8, 4))
        spectra_mode = ttk.Combobox(
            bar,
            textvariable=self.spectra_mode_var,
            state="readonly",
            width=38,
            values=("Raw reflectance", "Illuminated spectrum S(lambda) * R(lambda)"),
        )
        spectra_mode.pack(side="right")
        spectra_mode.bind("<<ComboboxSelected>>", lambda _e: self.update_plots())

    def _build_plot_notebook(self, parent) -> None:
        ttk = self.ttk
        self.notebook = ttk.Notebook(parent)
        self.notebook.grid(row=1, column=0, sticky="nsew")

        self.main_tab = ttk.Frame(self.notebook)
        self.zoom_tab = ttk.Frame(self.notebook)
        self.spectra_tab = ttk.Frame(self.notebook)
        self.results_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.main_tab, text="Main CIE")
        self.notebook.add(self.zoom_tab, text="Zoom")
        self.notebook.add(self.spectra_tab, text="Spectra")
        self.notebook.add(self.results_tab, text="Results")

        for tab in (self.main_tab, self.zoom_tab, self.spectra_tab, self.results_tab):
            tab.rowconfigure(0, weight=1)
            tab.columnconfigure(0, weight=1)

        columns = ("name", "X", "Y", "Z", "x", "y", "u", "v", "hex", "gamut", "warning")
        self.results_tree = ttk.Treeview(self.results_tab, columns=columns, show="headings")
        headers = {
            "name": "Name", "X": "X", "Y": "Y", "Z": "Z", "x": "x", "y": "y",
            "u": "u'", "v": "v'", "hex": "HEX", "gamut": "sRGB", "warning": "Warning",
        }
        for column in columns:
            self.results_tree.heading(column, text=headers[column])
            width = 80
            if column == "name": width = 160
            if column == "warning": width = 280
            self.results_tree.column(column, width=width, anchor="w" if column in ("name", "warning") else "e")
        self.results_tree.grid(row=0, column=0, sticky="nsew")
        ybar = ttk.Scrollbar(self.results_tab, orient="vertical", command=self.results_tree.yview)
        ybar.grid(row=0, column=1, sticky="ns")
        xbar = ttk.Scrollbar(self.results_tab, orient="horizontal", command=self.results_tree.xview)
        xbar.grid(row=1, column=0, sticky="ew")
        self.results_tree.configure(yscrollcommand=ybar.set, xscrollcommand=xbar.set)

    # ------------------------------------------------------------ Input list

    def add_files(self) -> None:
        from tkinter import filedialog

        paths = filedialog.askopenfilenames(
            title="Select FDTD TXT spectra",
            filetypes=[("TXT spectra", "*.txt"), ("All files", "*.*")],
        )
        self._add_paths(Path(path) for path in paths)

    def add_folder(self) -> None:
        from tkinter import filedialog, messagebox

        folder = filedialog.askdirectory(title="Select folder containing FDTD TXT spectra")
        if not folder:
            return
        try:
            paths = find_txt_files(folder, recursive=True)
        except Exception as exc:
            messagebox.showerror("Folder error", str(exc))
            return
        self._add_paths(paths)

    def _add_paths(self, paths: Iterable[Path]) -> None:
        known = {entry.path.resolve() for entry in self.entries}
        added = 0
        for path in paths:
            path = Path(path)
            if not path.is_file() or path.suffix.lower() != ".txt":
                continue
            resolved = path.resolve()
            if resolved in known:
                continue
            self.entries.append(SpectrumEntry(path=path, label=path.stem, marker="o"))
            known.add(resolved)
            added += 1
        self._assign_markers()
        self._refresh_file_tree()
        if added:
            self.status_var.set(f"Added {added} spectrum file(s). Click Analyze / Recalculate.")

    def remove_selected(self) -> None:
        indices = sorted(self._selected_indices(), reverse=True)
        for index in indices:
            del self.entries[index]
        self._assign_markers()
        self._refresh_file_tree()
        self.update_plots()

    def clear_all(self) -> None:
        self.entries.clear()
        self._refresh_file_tree()
        self._clear_plots()
        self._refresh_results_tree()
        self.status_var.set("Input list cleared.")

    def rename_selected(self) -> None:
        from tkinter import simpledialog

        indices = self._selected_indices()
        if len(indices) != 1:
            return
        entry = self.entries[indices[0]]
        new_name = simpledialog.askstring("Rename", "Display name:", initialvalue=entry.label, parent=self.root)
        if not new_name or not new_name.strip():
            return
        entry.label = new_name.strip()
        if entry.analysis is not None:
            entry.analysis = replace(
                entry.analysis,
                result=replace(entry.analysis.result, label=entry.label),
            )
        self._refresh_file_tree()
        self.update_plots()

    def move_selected(self, direction: int) -> None:
        indices = self._selected_indices()
        if len(indices) != 1:
            return
        index = indices[0]
        target = index + direction
        if target < 0 or target >= len(self.entries):
            return
        self.entries[index], self.entries[target] = self.entries[target], self.entries[index]
        self._assign_markers()
        self._refresh_file_tree(select_index=target)
        self.update_plots()

    def _selected_indices(self) -> list[int]:
        result: list[int] = []
        for iid in self.file_tree.selection():
            try:
                result.append(int(iid))
            except ValueError:
                pass
        return sorted(set(result))

    def _assign_markers(self) -> None:
        for index, entry in enumerate(self.entries):
            entry.marker = marker_for_index(index)

    def _refresh_file_tree(self, select_index: int | None = None) -> None:
        for iid in self.file_tree.get_children():
            self.file_tree.delete(iid)
        for index, entry in enumerate(self.entries):
            if entry.analysis is not None:
                u_text = f"{entry.analysis.result.u_prime:.5f}"
                v_text = f"{entry.analysis.result.v_prime:.5f}"
            elif entry.error:
                u_text = "ERR"
                v_text = "ERR"
            else:
                u_text = "-"
                v_text = "-"
            self.file_tree.insert("", "end", iid=str(index), values=(entry.label, entry.marker, u_text, v_text))
        if select_index is not None and 0 <= select_index < len(self.entries):
            self.file_tree.selection_set(str(select_index))

    # --------------------------------------------------------- Colourimetry

    def _on_illuminant_changed(self, _event=None) -> None:
        if self.illuminant_var.get() == ILLUMINANT_CUSTOM and not self.custom_illuminant_var.get():
            self.choose_custom_illuminant()
        self._mark_analysis_stale()

    def choose_custom_illuminant(self) -> None:
        from tkinter import filedialog

        path = filedialog.askopenfilename(
            title="Select custom illuminant spectrum",
            filetypes=[("Spectrum files", "*.txt *.csv"), ("TXT", "*.txt"), ("CSV", "*.csv"), ("All files", "*.*")],
        )
        if path:
            self.custom_illuminant_var.set(path)
            self.illuminant_var.set(ILLUMINANT_CUSTOM)
            self._mark_analysis_stale()

    def _mark_analysis_stale(self) -> None:
        if self.entries:
            self.status_var.set("Colorimetry settings changed. Click Analyze / Recalculate.")

    def _analysis_settings(self) -> AnalysisSettings:
        custom_path = Path(self.custom_illuminant_var.get()) if self.custom_illuminant_var.get() else None
        return AnalysisSettings(
            observer=self.observer_var.get(),
            illuminant=self.illuminant_var.get(),
            custom_illuminant_path=custom_path,
        )

    def start_analysis(self) -> None:
        from tkinter import messagebox

        if self._analysis_running:
            return
        if not self.entries:
            messagebox.showwarning("No input", "Add at least one FDTD TXT spectrum first.")
            return
        if self.illuminant_var.get() == ILLUMINANT_CUSTOM and not self.custom_illuminant_var.get():
            messagebox.showerror("Custom illuminant", "Select a custom illuminant spectrum first.")
            return

        settings = self._analysis_settings()
        self._analysis_running = True
        self.analyze_button.configure(state="disabled")
        self.status_var.set(f"Analyzing {len(self.entries)} spectrum file(s)...")

        thread = threading.Thread(target=self._analysis_worker, args=(settings,), daemon=True)
        thread.start()

    def _analysis_worker(self, settings: AnalysisSettings) -> None:
        success = 0
        failed = 0
        for index, entry in enumerate(self.entries, start=1):
            try:
                entry.analysis = analyze_fdtd_file(entry.path, label=entry.label, settings=settings)
                entry.error = ""
                success += 1
            except Exception as exc:
                entry.analysis = None
                entry.error = str(exc)
                failed += 1
            self.root.after(0, lambda i=index: self.status_var.set(f"Analyzing {i}/{len(self.entries)}..."))

        self.root.after(0, lambda: self._finish_analysis(success, failed))

    def _finish_analysis(self, success: int, failed: int) -> None:
        from tkinter import messagebox

        self._analysis_running = False
        self.analyze_button.configure(state="normal")
        self._refresh_file_tree()
        self._refresh_results_tree()
        self.update_plots()
        if failed:
            self.status_var.set(f"Analysis complete: {success} succeeded, {failed} failed.")
            errors = [f"{entry.path.name}: {entry.error}" for entry in self.entries if entry.error]
            preview = "\n".join(errors[:8])
            if len(errors) > 8:
                preview += f"\n... and {len(errors) - 8} more"
            messagebox.showwarning("Analysis completed with errors", preview)
        else:
            self.status_var.set(f"Analysis complete: {success} spectrum file(s).")

    # -------------------------------------------------------------- Plotting

    def _plot_style(self) -> PlotStyle:
        try:
            style = PlotStyle(
                font_family=self.font_family_var.get().strip() or "Arial",
                axis_font_size=float(self.axis_font_size_var.get()),
                tick_font_size=float(self.tick_font_size_var.get()),
                legend_font_size=float(self.legend_font_size_var.get()),
                marker_size=float(self.marker_size_var.get()),
                axis_line_width=float(self.axis_line_width_var.get()),
                figure_width_in=float(self.figure_width_var.get()),
                figure_height_in=float(self.figure_height_var.get()),
                dpi=int(float(self.dpi_var.get())),
                show_grid=self.show_grid_var.get(),
                show_spectral_locus=self.show_locus_var.get(),
                show_diagram_colours=self.show_background_var.get(),
                show_point_labels=self.show_labels_var.get(),
                show_legend=self.show_legend_var.get(),
                show_zoom_rectangle=self.show_zoom_rectangle_var.get(),
                match_legend_text_color=self.legend_text_color_var.get(),
                improve_light_text_readability=self.readable_legend_var.get(),
            )
        except ValueError as exc:
            raise ValueError("Plot sizes, font sizes, line widths, and DPI must be numeric.") from exc
        if min(style.axis_font_size, style.tick_font_size, style.legend_font_size, style.marker_size, style.axis_line_width, style.figure_width_in, style.figure_height_in, style.dpi) <= 0:
            raise ValueError("Plot numeric settings must be greater than zero.")
        return style

    def _main_cie_range(self) -> CIERange:
        if self.use_default_main_range_var.get():
            return CIERange()
        try:
            region = CIERange(
                float(self.main_u_min_var.get()),
                float(self.main_u_max_var.get()),
                float(self.main_v_min_var.get()),
                float(self.main_v_max_var.get()),
            ).normalized()
        except ValueError as exc:
            raise ValueError("Main CIE range coordinates must be numeric.") from exc
        if region.width <= 0 or region.height <= 0:
            raise ValueError("Main CIE range must have non-zero width and height.")
        return region

    def _zoom_region(self) -> ZoomRegion:
        try:
            region = ZoomRegion(
                float(self.zoom_u_min_var.get()),
                float(self.zoom_u_max_var.get()),
                float(self.zoom_v_min_var.get()),
                float(self.zoom_v_max_var.get()),
            ).normalized()
        except ValueError as exc:
            raise ValueError("Zoom coordinates must be numeric.") from exc
        if region.width <= 0 or region.height <= 0:
            raise ValueError("Zoom region must have non-zero width and height.")
        return region

    def _valid_analyses(self) -> list[AnalyzedSpectrum]:
        return [entry.analysis for entry in self.entries if entry.analysis is not None]

    def _appearance_map(self) -> dict[str, PointAppearance]:
        return {
            str(entry.path): PointAppearance(entry.marker, entry.label)
            for entry in self.entries
            if entry.analysis is not None
        }

    def update_plots(self) -> None:
        from tkinter import messagebox

        analyses = self._valid_analyses()
        if not analyses:
            self._clear_plots()
            return
        try:
            style = self._plot_style()
            main_range = self._main_cie_range()
            region = self._zoom_region()
            appearances = self._appearance_map()
            # Keep the CIE background consistent with the observer that was
            # actually used to calculate the current points. If the user has
            # changed observer/illuminant controls, the status bar asks for a
            # recalculation rather than silently mixing definitions.
            observer = analyses[0].result.observer

            main = plot_cie1976_main(
                analyses,
                style=style,
                observer=observer,
                main_range=main_range,
                zoom_region=region,
                appearances=appearances,
            )
            zoom = plot_cie1976_zoom(
                analyses,
                region,
                style=style,
                observer=observer,
                appearances=appearances,
                include_only_points_inside=True,
            )
            spectra = plot_spectra(
                analyses,
                style=style,
                mode=self.spectra_mode_var.get(),
                appearances=appearances,
            )
        except Exception as exc:
            messagebox.showerror("Plot error", str(exc))
            return

        self._set_figure("main", main)
        self._set_figure("zoom", zoom)
        self._set_figure("spectra", spectra)

    def _set_figure(self, kind: str, figure) -> None:
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        import matplotlib.pyplot as plt

        tab = {"main": self.main_tab, "zoom": self.zoom_tab, "spectra": self.spectra_tab}[kind]
        canvas_attr = f"{kind}_canvas"
        toolbar_attr = f"{kind}_toolbar"
        figure_attr = f"{kind}_figure"

        old_canvas = getattr(self, canvas_attr)
        old_toolbar = getattr(self, toolbar_attr)
        old_figure = getattr(self, figure_attr)
        if old_toolbar is not None:
            old_toolbar.destroy()
        if old_canvas is not None:
            old_canvas.get_tk_widget().destroy()
        if old_figure is not None:
            plt.close(old_figure)

        canvas = FigureCanvasTkAgg(figure, master=tab)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        toolbar_frame = self.ttk.Frame(tab)
        toolbar_frame.grid(row=1, column=0, sticky="ew")
        toolbar = NavigationToolbar2Tk(canvas, toolbar_frame, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side="left", fill="x")

        setattr(self, canvas_attr, canvas)
        setattr(self, toolbar_attr, toolbar_frame)
        setattr(self, figure_attr, figure)

    def _clear_plots(self) -> None:
        import matplotlib.pyplot as plt

        for kind in ("main", "zoom", "spectra"):
            canvas = getattr(self, f"{kind}_canvas")
            toolbar = getattr(self, f"{kind}_toolbar")
            figure = getattr(self, f"{kind}_figure")
            if toolbar is not None:
                toolbar.destroy()
                setattr(self, f"{kind}_toolbar", None)
            if canvas is not None:
                canvas.get_tk_widget().destroy()
                setattr(self, f"{kind}_canvas", None)
            if figure is not None:
                plt.close(figure)
                setattr(self, f"{kind}_figure", None)

    def activate_zoom_selector(self) -> None:
        from tkinter import messagebox
        from matplotlib.widgets import RectangleSelector

        if self.main_figure is None or not self.main_figure.axes:
            messagebox.showinfo("Zoom selection", "Analyze spectra and create the main CIE plot first.")
            return

        ax = self.main_figure.axes[0]
        if self._zoom_selector is not None:
            try:
                self._zoom_selector.set_active(False)
            except Exception:
                pass

        self._zoom_selector = RectangleSelector(
            ax,
            self._on_zoom_selected,
            useblit=True,
            button=[1],
            minspanx=0.001,
            minspany=0.001,
            spancoords="data",
            interactive=True,
        )
        self.status_var.set("Drag a rectangle on the Main CIE plot to define the zoom region.")

    def _on_zoom_selected(self, eclick, erelease) -> None:
        if None in (eclick.xdata, eclick.ydata, erelease.xdata, erelease.ydata):
            return
        u0, u1 = sorted((float(eclick.xdata), float(erelease.xdata)))
        v0, v1 = sorted((float(eclick.ydata), float(erelease.ydata)))
        self.zoom_u_min_var.set(f"{u0:.6g}")
        self.zoom_u_max_var.set(f"{u1:.6g}")
        self.zoom_v_min_var.set(f"{v0:.6g}")
        self.zoom_v_max_var.set(f"{v1:.6g}")
        if self._zoom_selector is not None:
            self._zoom_selector.set_active(False)
            self._zoom_selector = None
        self.status_var.set("Zoom region updated from the Main CIE plot.")
        self.root.after(0, self.update_plots)

    # --------------------------------------------------------------- Results

    def _refresh_results_tree(self) -> None:
        for iid in self.results_tree.get_children():
            self.results_tree.delete(iid)
        for index, entry in enumerate(self.entries):
            if entry.analysis is None:
                if entry.error:
                    self.results_tree.insert("", "end", iid=str(index), values=(entry.label, "", "", "", "", "", "", "", "", "ERROR", entry.error))
                continue
            r = entry.analysis.result
            gamut = "clipped" if r.out_of_srgb_gamut else "in gamut"
            self.results_tree.insert(
                "",
                "end",
                iid=str(index),
                values=(
                    entry.label,
                    f"{r.X:.6g}", f"{r.Y:.6g}", f"{r.Z:.6g}",
                    f"{r.x:.6f}", f"{r.y:.6f}", f"{r.u_prime:.6f}", f"{r.v_prime:.6f}",
                    r.hex_color, gamut, r.warning,
                ),
            )

    # --------------------------------------------------------------- Export

    def export_plot(self, kind: str) -> None:
        from tkinter import filedialog, messagebox

        figure = getattr(self, f"{kind}_figure")
        if figure is None:
            messagebox.showwarning("Nothing to export", f"No {kind} figure is available yet.")
            return
        path = filedialog.asksaveasfilename(
            title=f"Export {kind} figure",
            defaultextension=".png",
            filetypes=[
                ("PNG", "*.png"),
                ("TIFF", "*.tif *.tiff"),
                ("SVG", "*.svg"),
                ("PDF", "*.pdf"),
                ("All files", "*.*"),
            ],
        )
        if not path:
            return
        try:
            dpi = int(float(self.dpi_var.get()))
            export_figure(figure, path, dpi=dpi)
        except Exception as exc:
            messagebox.showerror("Export error", str(exc))
            return
        self.status_var.set(f"Exported {kind} figure: {path}")

    def export_results(self) -> None:
        from tkinter import filedialog, messagebox

        analyses = self._valid_analyses()
        if not analyses:
            messagebox.showwarning("Nothing to export", "No successful analysis results are available.")
            return
        path = filedialog.asksaveasfilename(
            title="Export colourimetric results",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            export_results_csv(path, analyses)
        except Exception as exc:
            messagebox.showerror("Export error", str(exc))
            return
        self.status_var.set(f"Exported results: {path}")


def launch_gui() -> None:
    """Launch the SPEC2CIE desktop application."""

    import tkinter as tk

    root = tk.Tk()
    SPEC2CIEApp(root)
    root.mainloop()


if __name__ == "__main__":
    launch_gui()
