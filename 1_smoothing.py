import numpy as np
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
import smoothing 
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.pyplot as plt
from lightkurve import KeplerLightCurveFile
import nancleaner as nc
import quarters as qs
import argparse, os, warnings

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
parser.add_argument('-s', '--smoothing', dest='kern', default=100, type=int, help='Gaussian smoothing kernel')
parser.add_argument('-c', '--clip', dest='inp', default=3, type=int, help='Outlier clipping level')
parser.add_argument('-p', '--plots', dest='show', default=False, type=bool, help='Show plots?')

params = parser.parse_args()

quarterlist = qs.getallquarters(params.kic)

time = np.zeros(0)
sap_flux = np.zeros(0)

for q in quarterlist:
   lc = KeplerLightCurveFile.from_archive(params.kic, quarter=q)

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
plt.title(f'{params.kic}')

clipped_flux = [sap_flux[i] for i in range(len(colours)) if colours[i] == 1]
clipped_time = [time[i] for i in range(len(colours)) if colours[i] == 1]
discarded_flux = [sap_flux[i] for i in range(len(colours)) if colours[i] == 0]
discarded_time = [time[i] for i in range(len(colours)) if colours[i] == 0]

plt.plot(discarded_time, discarded_flux, 'kx', alpha=0.1, markersize=3)
plt.plot(clipped_time, clipped_flux, 'ko', markersize=1)
plt.savefig(f'kic{params.kic}_smooth.png')

# export smoothed and clipped data as .dat file
exportblend = np.array([clipped_time, clipped_flux])
exportblend = np.transpose(exportblend)
np.savetxt(f'kic{params.kic}_lc.dat', exportblend, delimiter=' ', header=f'Smoothed and clipped light curve for KIC{params.kic}')

if params.show == True:
   plt.show()
else:
   pass