
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