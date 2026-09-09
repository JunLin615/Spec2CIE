
# SPEC2CIE

SPEC2CIE converts FDTD-exported reflectance spectra directly into CIE 1931 XYZ/xy and CIE 1976 `u'v'` coordinates, and provides publication-oriented plotting for multiple spectra.

The FDTD reader assumes the first numeric column is wavelength in **meters** and the second numeric column is reflectance. Internally, wavelengths are converted to nanometers and resampled to **360-830 nm at 1 nm intervals** for colourimetric integration. Extrapolation is disabled.

## Main features

- Directly reads FDTD TXT spectra; no intermediate CSV is required.
- Multiple TXT files or a recursively scanned folder can be analyzed together.
- CIE 1931 2° observer (default) and CIE 1964 10° observer.
- Selectable illumination:
  - Equal-energy illuminant E (default; closest to treating the reflectance spectrum itself as the spectral weighting)
  - CIE D65 daylight
  - CIE D50 daylight
  - CIE A tungsten/halogen illumination
  - 3000 K / 4000 K / 6500 K blackbody approximations for generic warm/neutral/cool white lighting
  - custom two-column illuminant spectrum
- Main CIE 1976 `u'v'` diagram with multiple points.
- Long automatic marker sequence; markers are assigned by file order.
- Marker colours approximate the calculated object colour on an sRGB display.
- Legend text can follow the point colour, with optional automatic darkening for readability.
- Adjustable Arial font, font sizes, marker size, axis width, figure size, DPI, grid, labels, legend, spectral locus, and chromaticity background.
- Interactive zoom-region selection on the main CIE diagram plus an independent zoom figure.
- Separate **Main CIE**, **Zoom**, **Spectra**, and **Results** tabs.
- Spectra tab can show either raw reflectance or relative `S(lambda) * R(lambda)`.
- Main, zoom, and spectrum figures can be exported to PNG, TIFF, SVG, or PDF.
- Results can be exported to CSV with XYZ, xy, u'v', approximate sRGB/HEX values, gamut status, and warnings.

## Installation

From the project root:

```bash
pip install -e .
```

Then launch the main application:

```bash
spec2cie
```

or:

```bash
python -m SPEC2CIE.cie1976_fdtd_viewer
```

The original FDTD TXT -> interpolated CSV utility remains available as:

```bash
spec2cie-csv
```

## Basic workflow

1. Click **Add Files** or **Add Folder** and select one or more FDTD TXT spectra.
2. Choose an observer and illuminant. **Equal-energy E** is the default.
3. Click **Analyze / Recalculate**.
4. Inspect the **Main CIE**, **Zoom**, **Spectra**, and **Results** tabs.
5. To define a detailed region, click **Select Region on Main Plot** and drag a rectangle on the main CIE figure.
6. Adjust plot settings and click **Update Plot Style** when needed.
7. Export the main/zoom plot or the results CSV.

## Python API

The GUI is only a front end. Core functions can be imported independently:

```python
from SPEC2CIE import AnalysisSettings, analyze_fdtd_file

result = analyze_fdtd_file(
    "R-30nm.txt",
    settings=AnalysisSettings(),
)

print(result.result.u_prime, result.result.v_prime)
```

The existing spectrum-processing module is reused internally:

```python
from SPEC2CIE import load_spectrum_txt, resample_spectrum

spectrum = load_spectrum_txt("R-30nm.txt")
visible = resample_spectrum(spectrum, start_nm=360, end_nm=830, step_nm=1)
```

## Physical interpretation

For a reflective sample, SPEC2CIE calculates object colour from the selected illuminant spectral power distribution `S(lambda)` and FDTD reflectance `R(lambda)`. The CIE coordinates therefore represent **predicted visual chromaticity under the selected illumination**.

The displayed RGB/HEX colour is an sRGB visualization aid. It uses chromatic adaptation to the sRGB D65 display white and clips out-of-gamut colours when necessary. It is **not** a microscope-camera simulation. Camera prediction would additionally require the camera sensor spectral sensitivities, exposure, white balance, colour matrix, and other imaging-pipeline parameters.

The generic 3000/4000/6500 K options are blackbody approximations and should not be interpreted as exact spectra of real LED lamps. For a specific microscope lamp or room light, use **Custom illuminant spectrum** with a measured SPD.

## Colour-science basis

The project uses the open-source `colour-science` package for standard CIE colour-matching functions, standard illuminants, spectral integration, chromatic adaptation, and the CIE diagram background.

Reference resources:

- CIE 1931 2° colour matching functions: https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer
- CIE 1964 10° colour matching functions: https://cie.co.at/datatable/cie-1964-colour-matching-functions-10-degree-observer
- CIE standard illuminants: https://www.cie.co.at/publications/colorimetry-part-2-cie-standard-illuminants-0
- Colour documentation: https://colour.readthedocs.io/
- Reference online CIE 1976 calculator: https://sciapps.sci-sim.com/CIE1976.html


## txt_to_interpolated_csv_gui
FDTD Reflectance Spectrum to CIELUV CSV Converter

This Python tool converts reflectance spectra exported from FDTD simulations into a CSV format suitable for the CIELUV 1976 (u′v′) Chromaticity Diagram Calculator.

The FDTD TXT file is expected to contain wavelength in meters in the first column and the spectrum value in the second column. The converter changes wavelength to nanometers, sorts the data, performs linear interpolation, and writes a headerless two-column CSV:

380,0.123456
381,0.124321
...

By default, the output range is 380–800 nm with a 1 nm interval. These values can be changed in the GUI. The online CIELUV calculator accepts comma- or tab-delimited spectrum data and recommends wavelength intervals such as 1, 2, or 5 nm.

Usage

Run:

python txt_to_interpolated_csv_gui.py

Then:

Select a single FDTD TXT file or a folder containing multiple TXT files.

Select the output folder.

Adjust Start, End, and Step if necessary.

Click Convert.

Upload the generated CSV file to the CIELUV 1976 online calculator.

For folder input, all TXT files are processed recursively and the original subfolder structure is preserved. Output files keep the original filename and use the .csv extension.

Use as a Python Module

The conversion functions can also be imported without launching the GUI:

from txt_to_interpolated_csv_gui import convert_fdtd_txt_to_csv

convert_fdtd_txt_to_csv(
    "spectrum.txt",
    "spectrum.csv",
    start_nm=380,
    end_nm=800,
    step_nm=1,
)

Note: this utility only performs wavelength-unit conversion, interpolation, and CSV formatting. It does not apply an illuminant spectrum or other colorimetric corrections to reflectance data.