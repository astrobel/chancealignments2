import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits as pyfits
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
from astropy.timeseries import LombScargle
import translate as tr
import argparse, warnings

warnings.simplefilter('ignore', category=UserWarning) # for font conflicts on my system, at least
warnings.simplefilter('ignore', category=RuntimeWarning) # may encounter some divide by zeros in process

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

parser = argparse.ArgumentParser(description='Phase light curve on highest-amplitude frequency, or given frequency.')
parser.add_argument('-k', '--kic', required=True, type=int, help='KIC ID')
parser.add_argument('-t', '--timecadence', default='long', choices=['long', 'short'], type=str, help='Cadence of data to use')
parser.add_argument('-f', '--foldfreq', dest='fold', default=None, type=float, help='Frequency for phase folding in microHertz or cycles per day (use -u True)')
parser.add_argument('-u', '--unitscpd', dest='cpd', default=False, type=bool, help='Use cycles per day instead of microhertz?')
parser.add_argument('-o', '--oversampling', dest='over', default=5, type=int, help='LSP oversampling factor')
parser.add_argument('-n', '--nyquistfactor', dest='nyq', default=1, type=float, help='LSP Nyquist factor')
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
      print('Wrong KIC number? Or try running 1_smoothing.py first!')
      sys.exit()
times = lc[:,0]
ampls = lc[:,1]

if params.fold == None:

   ofac = params.over
   hifac = params.nyq
   
   frequencies, power_spectrum = LombScargle(np.asarray(times), np.asarray(ampls)).autopower(method='fast', normalization='psd', samples_per_peak=ofac, nyquist_factor=hifac)
   hifac *= (283/11.57)/max(frequencies)
   frequencies, power_spectrum = LombScargle(np.asarray(times), np.asarray(ampls)).autopower(method='fast', normalization='psd', samples_per_peak=ofac, nyquist_factor=hifac)
   power_spectrum = power_spectrum * 4. / len(times)
   power_spectrum = np.sqrt(power_spectrum)
   power_spectrum *= 1e6
   frequencies *= 11.57

   foldfreq = frequencies[power_spectrum.argmax()]

else:
   if params.cpd == False:
      foldfreq = params.fold
   elif params.cpd == True:
      foldfreq = params.fold * 11.57

tohz = foldfreq * 1e-6
tos = 1/tohz
foldper = tos/86400

binnum = 100
binnum2 = 1000
binsize = foldper / binnum
binsize2 = foldper / binnum2

phasedtimearray = np.zeros(len(times))
finalampls = np.zeros(binnum)
amplcounts = np.zeros(binnum)
finalampls2 = np.zeros(binnum2)
amplcounts2 = np.zeros(binnum2)

for i, val in enumerate(times):
   phasedtime = val % foldper
   newphasedtime = tr.translate(phasedtime, 0, foldper, 0, 1)
   phasedtimearray[i] = newphasedtime
   bindex = int((phasedtime - (phasedtime % binsize)) / binsize - 1)
   finalampls[bindex] += ampls[i]
   amplcounts[bindex] += 1
   bindex2 = int((phasedtime - (phasedtime % binsize2)) / binsize2 - 1)
   finalampls2[bindex2] += ampls[i]
   amplcounts2[bindex2] += 1
   
finalampls = np.divide(finalampls, amplcounts)
finalampls2 = np.divide(finalampls2, amplcounts2)

finaltimes = np.histogram(phasedtimearray, bins=binnum-1, range=(0,1))
finaltimes2 = np.histogram(phasedtimearray, bins=binnum2-1, range=(0,1))
# finalampls = np.histogram(ampls, bins=binnum)

plt.figure(1)

plt.plot(finaltimes2[1], finalampls2, 'kx', markersize=3, alpha=0.3)
plt.plot(finaltimes2[1]+1, finalampls2, 'kx', markersize=3, alpha=0.3)
plt.plot(finaltimes[1], finalampls, 'ko', markersize=2)
plt.plot(finaltimes[1]+1, finalampls, 'ko', markersize=2)
plt.xlim(0, max(finaltimes2[1]+1))
plt.xlabel(f'Normalised Time mod {foldper:.2} days')
plt.ylabel('Fractional Intensity')
plt.title(f'{kic}')
plt.savefig(f'kic{kic}_phase_{cadence}.png')

if params.show == True:
   plt.show()
else:
   pass
