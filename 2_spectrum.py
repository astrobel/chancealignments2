import numpy as np
from astropy.stats import LombScargle
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse, warnings, os

warnings.simplefilter('ignore', category=UserWarning) # for font conflicts on my system, at least

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
      r'\usepackage{helvet}',
      r'\usepackage[EULERGREEK]{sansmath}',
      r'\sansmath'
]
mpl.rcParams['axes.formatter.useoffset'] = False
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True

parser = argparse.ArgumentParser(description='Prepare long cadence light curve for further analysis.')
parser.add_argument('-k', '--kic', required=True, type=int, help='KIC ID')
parser.add_argument('-o', '--oversampling', dest='over', default=5, type=int, help='LSP oversampling factor')
parser.add_argument('-n', '--nyqistfactor', dest='nyq', default=1, type=int, help='LSP Nyquist factor')
parser.add_argument('-p', '--plots', dest='show', default=False, type=bool, help='Show plots?')

params = parser.parse_args()

# read in light curve
while True:
   try:
      importblend = np.loadtxt(f'kic{params.kic}_lc.dat')
      break
   except OSError:
      print('Wrong KIC number? Or try running 1_smoothing.py first!')
      sys.exit()
clipped_time = importblend[:,0]
clipped_flux = importblend[:,1]

# fourier transform

ofac = params.over
hifac = params.nyq

frequencies, power_spectrum = LombScargle(np.asarray(clipped_time), np.asarray(clipped_flux)).autopower(method='fast', normalization='psd', samples_per_peak=ofac, nyquist_factor=hifac)
hifac *= (283/11.57)/max(frequencies)
frequencies, power_spectrum = LombScargle(np.asarray(clipped_time), np.asarray(clipped_flux)).autopower(method='fast', normalization='psd', samples_per_peak=ofac, nyquist_factor=hifac)
power_spectrum = power_spectrum * 4 / len(clipped_time)
power_spectrum = np.sqrt(power_spectrum)
power_spectrum *= 1e6
frequencies *= 11.57

# logarithmic and linear power spectra
fig1, (pslog, pslin) = plt.subplots(2, 1)
pslog.plot(frequencies, power_spectrum, 'k-', lw=0.5)
pslog.set_xlim(1, max(frequencies))
pslog.set_ylabel('Amplitude (ppm)')
pslog.set_xscale('log')
pslog.set_yscale('log')
pslog.set_title(f'{params.kic}')

pslin.plot(frequencies, power_spectrum, 'k-', lw=0.5)
pslin.set_xlim(1, max(frequencies))
pslin.set_ylim(ymin = 0)
pslin.set_xlabel('Frequency ($\mu$Hz)')
pslin.set_ylabel('Amplitude (ppm)')

plt.tight_layout()
fig1.savefig(f'kic{params.kic}_spectrum.png')

if params.show == True:
   plt.show()
else:
   pass