import numpy as np
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
import smoothing 
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.pyplot as plt
from lightkurve import search_lightcurvefile, LightkurveWarning
import nancleaner as nc
import argparse, warnings

warnings.simplefilter('ignore', category=UserWarning) # for font conflicts on my system, at least
warnings.simplefilter('ignore', category=LightkurveWarning) # if it's trying to download an empty quarter, this will be skipped in code

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
parser.add_argument('-t', '--timecadence', default='long', choices=['long', 'short'], type=str, help='Cadence of data to download')
parser.add_argument('-s', '--smoothing', dest='kern', default=100, type=int, help='Gaussian smoothing kernel, in days')
parser.add_argument('-c', '--clip', dest='inp', default=3, type=float, help='Outlier clipping level, in sigma')
parser.add_argument('-p', '--plots', dest='show', default=False, type=bool, help='Show plots?')

params = parser.parse_args()

kic = params.kic
cadence = params.timecadence

time = np.zeros(0)
sap_flux = np.zeros(0)

for q in np.arange(0,18):

   lc = search_lightcurvefile(f'KIC {kic}', quarter=q, cadence=cadence).download()
   
   if lc != None:
      table = lc.hdu[1].data

      sap_flux_temp = table['SAP_FLUX']
      time_temp = table['TIME']

      sap_flux_1, time_1 = nc.nancleaner2d(sap_flux_temp, time_temp)

      sap_flux_2, smth_flux = smoothing.gausssmooth(time_1, sap_flux_1, params.kern)

      time = np.append(time, time_1)
      sap_flux = np.append(sap_flux, sap_flux_2)
      
      continue

# scan for outliers
clip = params.inp * np.std(sap_flux)
meanflux = np.mean(sap_flux)

upperbound = meanflux + clip
lowerbound = meanflux - clip

colours = np.zeros(sap_flux.size)

for i, flux in enumerate(sap_flux):
   if flux < upperbound and flux > lowerbound:
      colours[i] = 1

clipped_flux = []
clipped_time = []
discarded_flux = []
discarded_time = []

# smoothed data

plt.figure(1)

plt.xlabel('Time (d)')
plt.ylabel('Fractional Intensity')
plt.title(f'{kic}')

clipped_flux = [sap_flux[i] for i in range(len(colours)) if colours[i] == 1]
clipped_time = [time[i] for i in range(len(colours)) if colours[i] == 1]
discarded_flux = [sap_flux[i] for i in range(len(colours)) if colours[i] == 0]
discarded_time = [time[i] for i in range(len(colours)) if colours[i] == 0]

plt.plot(discarded_time, discarded_flux, 'kx', alpha=0.1, markersize=3)
plt.plot(clipped_time, clipped_flux, 'ko', markersize=1)
plt.savefig(f'kic{kic}_lc_{cadence}.png')

# export smoothed and clipped data as .dat file
exportblend = np.array([clipped_time, clipped_flux])
exportblend = np.transpose(exportblend)
np.savetxt(f'kic{kic}_lc_{cadence}.dat', exportblend, delimiter=' ', header=f'Smoothed and clipped {cadence} cadence light curve for KIC{kic}')

if params.show == True:
   plt.show()
else:
   pass