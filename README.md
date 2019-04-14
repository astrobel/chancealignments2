# chancealignments2
Updated version of https://github.com/astrobel/chancealignments, now in Python 3! This collection of scripts is a fully-automated pipeline for close inspection of Kepler long cadence data and target pixel files, aimed at aiding the identification of chance alignments.

## Requirements

* numpy
* scipy
* astropy
* matplotlib
* lightkurve
* argparse

## Main files

### 1_lightcurve.py

Prepares light curve for analysis.

Arguments:
* `-k --kic`: KIC ID, required, takes integer
* `-s --smoothing`: Gaussian smoothing kernel, default 100 days, takes integer
* `-c --clip`: Outlier clipping level, default 3 sigma, takes integer
* `-p --plots`: Show plots, default False, takes bool

### 2_spectrum.py

Computes Lomb-Scargle periodogram of light curve. `1_lightcurve.py` _must_ be run first.

Arguments:
* `-k --kic`: KIC ID, required, takes integer
* `-o --oversampling`: LSP oversampling factor, default 5, takes integer
* `-n --nyquistfactor`: LSP Nyquist factor, default 1, takes float
* `-p --plots`: Show plots, default False, takes bool

### 3_phase.py

Phases light curve on highest-amplitude frequency, or chosen frequency. `1_lightcurve.py` _must_ be run first.

Arguments:
* `-k --kic`: KIC ID, required, takes integer
* `-f --foldfreq`: Frequency for phase folding, in &mu;Hertz, default None runs code on highest-amplitude frequency, takes float
* `-o --oversampling`: LSP oversampling factor, default 5, takes integer
* `-n --nyquistfactor`: LSP Nyquist factor, default 1, takes float
* `-p --plots`: Show plots, default False, takes bool

### 4_pixels.py

Examines light curves and amplitude spectra in each individual pixel for a given quarter. Code will exit if there is no data for chosen quarter.

Arguments:
* `-k --kic`: KIC ID, required, takes integer
* `-q --quarter`: Quarter to analyse, required, takes integer between 0 and 17
* `-s --smoothing`: Gaussian smoothing kernel, default 100 days, takes integer
* `-c --clip`: Outlier clipping level, default 3 sigma, takes integer
* `-o --oversampling`: LSP oversampling factor, default 5, takes integer
* `-n --nyquistfactor`: LSP Nyquist factor, default 1, takes float
* `-m --makepdf`: Option to create PDF that shows a close-up of the light curve and amplitude spectrum for each individual pixel, including a numbered plot at the beginning to identify each pixel, default False, takes bool
* `-e --export`: Export a data file containing light curve data for each individual pixel, default False, takes bool
* `-p --plots`: Show plots, default False, takes bool

### 5_comparison.py

Plots a side-by-side comparison of pixel image for one quarter and a 1' UKIRT image of the same area, the latter of which must be downloaded by the user. Code will exit if there is no data for chosen quarter.

Arguments:
* `-k --kic`: KIC ID, required, takes integer
* `-q --quarter`: Quarter to analyse, required, takes integer between 0 and 17
* `-u --ukirt`: UKIRT image filename, required, takes string
* `-r --refpix`: Plot reference pixel location from FITS header info on Kepler image, default False, takes bool
* `-p --plots`: Show plots, default False, takes bool

### 6_difference.py

Performs difference imaging on a given quarter, to find the pixel source of a signal. `1_lightcurve.py` _must_ be run first. Code will exit if there is no data for chosen quarter.

This code uses only time series data which falls within 5% (`tolerance` variable) either side of the peak and trough of the phase-folded light curve, or 20% of the light curve in total. Trough points are subtracted from peak points for each individual pixel and an average is taken to create the difference image. The plot returned is a side-by-side comparison of an average of all unaltered frames and the difference image.

Arguments:
* `-k --kic`: KIC ID, required, takes integer
* `-q --quarter`: Quarter to analyse, required, takes integer between 0 and 17
* `-f --foldfreq`: Frequency for phase folding, in &mu;Hertz, default None runs code on highest-amplitude frequency, takes float
* `-o --oversampling`: LSP oversampling factor, default 5, takes integer
* `-n --nyquistfactor`: LSP Nyquist factor, default 1, takes float
* `-r --refpix`: Plot reference pixel location from FITS header info on Kepler image, default False, takes bool
* `-p --plots`: Show plots, default False, takes bool

### 7_regionmap.py - NEW!

Maps Gaia DR2 sources to a given quarter's TPF postage stamp, and lists sources for easy identification. Code will exit if there is no data for chosen quarter.

Arguments:
* `-k --kic`: KIC ID, required, takes integer
* `-q --quarter`: Quarter to analyse, required, takes integer between 0 and 17
* `-s --sourceprint`: Print Gaia DR2 sources to console, default False, takes bool
* `-p --plots`: Show plots, default False, takes bool

## Auxiliary files

* nancleaner.py: module to remove NaN values from both SAP and TPF light curve data.
* outliers.py: module to handle clipping of outliers given a sigma value for bounding.
* smoothing.py: contains functions which perform the convolution of a Gaussian, boxcar, or polynomial kernel with time series data, and then divides the time series by the fit to normalise. Only Gaussian smoothing is used by the code, but the other functions are included for use at the user's discretion.
* translate.py: module to translate a time series from one range to another, used to normalise time series in `3_phase.py`.
