import numpy as np
from astropy.timeseries import LombScargle
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse, warnings, sys

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

parser = argparse.ArgumentParser(description='Lomb-Scargle periodogram of light curve.')
parser.add_argument('-k', '--kic', required=True, type=int, help='KIC ID')
parser.add_argument('-t', '--timecadence', default='long', choices=['long', 'short'], type=str, help='Cadence of data to use')
parser.add_argument('-s', '--sampling', dest='over', default=5, type=int, help='LSP oversampling factor')
parser.add_argument('-n', '--nyquistfactor', dest='nyq', default=1, type=float, help='LSP Nyquist factor')
parser.add_argument('-u', '--unitscpd', dest='cpd', default=False, type=bool, help='Use cycles per day instead of microhertz in plot?')
parser.add_argument('-p', '--plots', dest='show', default=False, type=bool, help='Show plots?')

params = parser.parse_args()

kic = params.kic
cadence = params.timecadence

# read in light curve
while True:
   try:
      lc = np.loadtxt(f'kic{kic}_lc_{cadence}.dat')
      break
   except OSError:
      print('Wrong KIC number or cadence? Try running 1_smoothing.py first!')
      sys.exit()
clipped_time = lc[:,0]
clipped_flux = lc[:,1]

# fourier transform

ofac = params.over
hifac = params.nyq

frequencies, power_spectrum = LombScargle(np.asarray(clipped_time), np.asarray(clipped_flux)).autopower(method='fast', normalization='psd', samples_per_peak=ofac, nyquist_factor=1)
if cadence == 'long':
   maxcpd = 24.4598
elif cadence == 'short':
   maxcpd = 734.0535
hifac *= maxcpd/max(frequencies)
frequencies, power_spectrum = LombScargle(np.asarray(clipped_time), np.asarray(clipped_flux)).autopower(method='fast', normalization='psd', samples_per_peak=ofac, nyquist_factor=hifac)
power_spectrum = power_spectrum * 4 / len(clipped_time)
power_spectrum = np.sqrt(power_spectrum)
power_spectrum *= 1e6
xaxismin = 0.1
if params.cpd == False:
   frequencies *= 11.57
   xaxismin = 1

# logarithmic and linear power spectra
fig1, (pslog, pslin) = plt.subplots(2, 1)
pslog.plot(frequencies, power_spectrum, 'k-', lw=0.5)
pslog.set_xlim(xaxismin, max(frequencies))
pslog.set_ylabel('Amplitude (ppm)')
pslog.set_xscale('log')
pslog.set_yscale('log')
pslog.set_title(f'{kic}')

pslin.plot(frequencies, power_spectrum, 'k-', lw=0.5)
pslin.set_xlim(xaxismin, max(frequencies))
pslin.set_ylim(ymin = 0)
if params.cpd == False:
   pslin.set_xlabel('Frequency ($\mu$Hz)')
elif params.cpd == True:
   pslin.set_xlabel('Frequency (cd)')
pslin.set_ylabel('Amplitude (ppm)')

plt.tight_layout()
fig1.savefig(f'kic{kic}_spectrum_{cadence}.png')

if params.show == True:
   plt.show()
else:
   pass