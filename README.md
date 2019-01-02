# chancealignments2
Updated version of https://github.com/astrobel/chancealignments, now in Python 3!

## Requirements

* numpy
* scipy
* astropy
* lightkurve

To run 5_comparison.py, it is suggested that you create separate directories for each target analysed

## Main files

### 1_smoothing.py

Prepare light curve for analysis.

Arguments:
* -k --kic: KIC ID, required, takes integer
* -s --smoothing: Gaussian smoothing kernel, default 100 days, takes integer
* -c --clip: Outlier clipping level, default 3 sigma, takes integer
* -p --plots: Show plots, default False, takes bool

### 2_spectrum.py

Compute Lomb-Scargle periodogram of light curve. `1_smoothing.py` _must_ be run first.

Arguments:
* -k --kic: KIC ID, required, takes integer
* -o --oversampling: LSP oversampling factor, default 5, takes integer
* -n --nyquistfactor: LSP Nyquist factor, default 1, takes float
* -p --plots: Show plots, default False, takes bool

### 3_phase.py

Phase light curve on highest-amplitude frequency, or chosen frequency. `1_smoothing.py` _must_ be run first.

Arguments:
* -k --kic: KIC ID, required, takes integer
* -f --foldfreq: Frequency for phase folding, in &mu;Hertz, default None runs code on highest-amplitude frequency, takes float
* -o --oversampling: LSP oversampling factor, default 5, takes integer
* -n --nyquistfactor: LSP Nyquist factor, default 1, takes float
* -p --plots: Show plots, default False, takes bool

## Auxiliary files

Coming soon!
