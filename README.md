# chancealignments2
Updated version of https://github.com/astrobel/chancealignments, now in Python 3!

## Requirements

* numpy
* scipy
* astropy
* lightkurve

## Main files

### 1_smoothing.py

Prepare light curve for analysis.

Arguments:
* `-k --kic`: KIC ID, required, takes integer
* `-s --smoothing`: Gaussian smoothing kernel, default 100 days, takes integer
* `-c --clip`: Outlier clipping level, default 3 sigma, takes integer
* `-p --plots`: Show plots, default False, takes bool

### 2_spectrum.py

Compute Lomb-Scargle periodogram of light curve. `1_smoothing.py` _must_ be run first.

Arguments:
* `-k --kic`: KIC ID, required, takes integer
* `-o --oversampling`: LSP oversampling factor, default 5, takes integer
* `-n --nyquistfactor`: LSP Nyquist factor, default 1, takes float
* `-p --plots`: Show plots, default False, takes bool

### 3_phase.py

Phase light curve on highest-amplitude frequency, or chosen frequency. `1_smoothing.py` _must_ be run first.

Arguments:
* `-k --kic`: KIC ID, required, takes integer
* `-f --foldfreq`: Frequency for phase folding, in &mu;Hertz, default None runs code on highest-amplitude frequency, takes float
* `-o --oversampling`: LSP oversampling factor, default 5, takes integer
* `-n --nyquistfactor`: LSP Nyquist factor, default 1, takes float
* `-p --plots`: Show plots, default False, takes bool

### 4_pixels.py

Examine light curves and amplitude spectra in each individual pixel for one quarter at a time.

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

Side-by-side comparison of pixel image for one quarter and a 1' UKIRT image of the same area, the latter of which must be downloaded by the user.

Arguments:
* `-k --kic`: KIC ID, required, takes integer
* `-q --quarter`: Quarter to analyse, required, takes integer between 0 and 17
* `-u --ukirt`: UKIRT image filename, required, takes string
* `-r --refpix`: Plot reference pixel location from FITS header info on Kepler image, default False, takes bool
* `-p --plots`: Show plots, default False, takes bool

## Auxiliary files

* nancleaner.py: module to remove NaN values from both SAP and TPF light curve data.
* outliers.py: module to handle clipping of outliers given a sigma value for bounding.
* quarters.py: module to fetch the availability of a quarter for a given target, or a list of all available quarters for that target.
* smoothing.py: contains functions which perform the convolution of a Gaussian, boxcar, or polynomial kernel with time series data, and then divides the time series by the fit to normalise. Only Gaussian smoothing is used by the code, but the other functions are included for use at the user's discretion.
* translate.py: module to translate a time series from one range to another, used to normalise time series in 3_phase.py.