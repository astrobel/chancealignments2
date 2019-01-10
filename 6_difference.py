import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy import wcs
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
from astropy.utils.exceptions import AstropyWarning
from astropy.stats import LombScargle
from lightkurve import KeplerTargetPixelFile
import smoothing
import translate as tr
import quarters as qs
import nancleaner as nc
import argparse, warnings

warnings.simplefilter('ignore', category=AstropyWarning)
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

parser = argparse.ArgumentParser(description='Prepare long cadence light curve for further analysis.')
parser.add_argument('-k', '--kic', required=True, type=int, help='KIC ID')
parser.add_argument('-q', '--quarter', required=True, type=int, choices=range(0,18), help='Quarter to analyse')
parser.add_argument('-f', '--foldfreq', dest='fold', default=None, type=float, help='Frequency for phase folding in microHertz')
parser.add_argument('-o', '--oversampling', dest='over', default=5, type=int, help='LSP oversampling factor')
parser.add_argument('-n', '--nyquistfactor', dest='nyq', default=1, type=float, help='LSP Nyquist factor')
parser.add_argument('-r', '--refpix', default=False, type=bool, help='Plot location of reference pixel on image?')
parser.add_argument('-p', '--plots', dest='show', default=False, type=bool, help='Show plots?')

params = parser.parse_args()

q = params.quarter
kic = params.kic

tpf = KeplerTargetPixelFile.from_archive(kic, quarter=q)

channel = tpf.channel
obj_ra = tpf.ra
obj_dec = tpf.dec
table = tpf.hdu[1].data
flux1 = table['FLUX']
time1 = table['TIME']
hd1 = tpf.hdu[1].header
ysize = hd1['NAXIS2']
table2 = tpf.hdu[2].data
hd2 = tpf.hdu[2].header
x = hd2['NAXIS1']
y = hd2['NAXIS2']
refx = hd2['CRPIX1']
refy = hd2['CRPIX2']
xsize = x * y
temp2d = np.zeros((x, y))
w = wcs.WCS(hd1, keysel=['binary'])

if (channel%2) == 0:
   eo = 0
else:
   eo = 1

fluxnew, timenew = nc.nancleaner3d(flux1, time1)

avgflux = np.nanmean(fluxnew, axis=0)

# read in light curve
lc = np.loadtxt(f'kic{kic}_lc.dat')
times = lc[:,0]
ampls = lc[:,1]

if params.fold == None:

   ofac = params.over
   hifac = params.nyq
   
   frequencies, power_spectrum = LombScargle(np.asarray(times), np.asarray(ampls)).autopower(method='fast', normalization='psd', samples_per_peak=ofac, nyquist_factor=hifac)
   hifac *= (283/11.57)/max(frequencies)
   frequencies, power_spectrum = LombScargle(np.asarray(times), np.asarray(ampls)).autopower(method='fast', normalization='psd', samples_per_peak=ofac, nyquist_factor=hifac)
   power_spectrum = power_spectrum * 4 / len(times)
   power_spectrum = np.sqrt(power_spectrum)
   power_spectrum *= 1e6
   frequencies *= 11.57

   foldfreq = frequencies[power_spectrum.argmax()]

else:
   foldfreq = params.fold

tohz = foldfreq * 1e-6
tos = 1/tohz
foldper = tos/86400

binnum = 100
binsize = foldper / binnum

phasedtimearray = np.zeros(len(times))
finalampls = np.zeros(binnum)
amplcounts = np.zeros(binnum)

for i, val in enumerate(times):
   phasedtime = val % foldper
   newphasedtime = tr.translate(phasedtime, 0, foldper, 0, 1)
   phasedtimearray[i] = newphasedtime
   bindex = int((phasedtime - (phasedtime % binsize)) / binsize - 1)
   finalampls[bindex] += ampls[i]
   amplcounts[bindex] += 1
   
finalampls = np.divide(finalampls, amplcounts)

finaltimes = np.histogram(phasedtimearray, bins=binnum-1, range=(0,1))

nanlist = []
for i, val in enumerate(finalampls):
   if np.isnan(val) == True:
      nanlist.append(i)

phasetime = np.delete(finaltimes[1], nanlist)
phaseflux = np.delete(finalampls, nanlist)

maxpos = np.argmax(phaseflux) / np.float64(len(phasetime))
minpos = np.argmin(phaseflux) / np.float64(len(phasetime))

# folding the pixel fluxes
time1 = time1 % foldper
for i, val in enumerate(time1):
   time1[i] = tr.translate(val, 0, foldper, 0, 1)

# find the fluxes for differencing
timeflags = np.zeros(len(time1))
tolerance = 0.05

for i, time in enumerate(time1):
   if minpos - tolerance <= 0:
      if time < minpos + tolerance or time > (minpos - tolerance)%1:
         timeflags[i] = -1
   elif maxpos + tolerance >= 1:
      if time < (minpos + tolerance)%1 or time > minpos - tolerance:
         timeflags[i] = -1
   else:
      if time < minpos + tolerance and time > minpos - tolerance:
         timeflags[i] = -1

   if maxpos - tolerance <= 0:
      if time < maxpos + tolerance or time > (maxpos - tolerance)%1:
         timeflags[i] = 1
   elif maxpos + tolerance >= 1:
      if time < (maxpos + tolerance)%1 or time > maxpos - tolerance:
         timeflags[i] = 1
   else:
      if time < maxpos + tolerance and time > maxpos - tolerance:
         timeflags[i] = 1

# establish number of frames to be used around the maxima and minima
highflags = np.where(timeflags>0)[0]
lowflags = np.where(timeflags<0)[0]
highdim = len(highflags)
lowdim = len(lowflags)

# 3d flux arrays to fill with frames
fluxhigh = np.zeros((highdim, y, x))
fluxlow = np.zeros((lowdim, y, x))

# get rid of nans
highnans = []
lownans = []
for i in range(highdim):
   fluxhigh[i,:] = flux1[highflags[i]]
   if np.isnan(fluxhigh[i,:]).all() == True:
      highnans.append(i)
for i in range(lowdim):
   fluxlow[i,:] = flux1[lowflags[i]]
   if np.isnan(fluxlow[i,:]).all() == True:
      lownans.append(i)

fluxhigh1 = np.delete(fluxhigh, highnans, axis=0)
fluxlow1 = np.delete(fluxlow, lownans, axis=0)

fluxdiff = np.abs(np.nanmean(fluxhigh1, axis=0) - np.nanmean(fluxlow1, axis=0))

imgflux = np.flipud(fluxdiff)
if eo == 0:
   imgflux = np.fliplr(imgflux)
avgflux = np.flipud(avgflux)
if eo == 0:
   avgflux = np.fliplr(avgflux)


### PLOTTING ###

plt.figure(1)

fig, (avgimg, diffimg) = plt.subplots(1, 2) 

crval = w.wcs.crval
north = crval + np.array([0, 6/3600.])
east = crval + np.array([ 6/3600., 0])

ncoords = np.vstack([crval, north])
ecoords = np.vstack([crval, east])
npixels = w.wcs_world2pix(ncoords , 0)
epixels = w.wcs_world2pix(ecoords , 0)
npixels[1, 1] = npixels[0, 1] - (npixels[1, 1] - npixels[0, 1]) # flip ud
epixels[1, 1] = epixels[0, 1] - (epixels[1, 1] - epixels[0, 1])
if eo == 0:
   npixels[1, 0] = npixels[0, 0] - (npixels[1, 0] - npixels[0, 0]) # flip lr
   epixels[1, 0] = epixels[0, 0] - (epixels[1, 0] - epixels[0, 0])
diffimg.plot(npixels[:,0], npixels[:,1], color='#0cb5ed')
diffimg.plot(epixels[:,0], epixels[:,1], '--', color='#0cb5ed')

left = avgimg.imshow(avgflux, cmap='YlOrRd')
avgimg.set_title('Average')
left.set_interpolation('nearest')
avgimg.set_xlim(-0.5, x-0.5)
avgimg.set_ylim(y-0.5, -0.5)
avgimg.plot(npixels[:,0], npixels[:,1], color='#0cb5ed')
avgimg.plot(epixels[:,0], epixels[:,1], '--', color='#0cb5ed')

left.axes.get_xaxis().set_ticklabels([])
left.axes.get_yaxis().set_ticklabels([])
left.axes.get_xaxis().set_ticks([])
left.axes.get_yaxis().set_ticks([])

avgimg.plot([25, 25], [25, 55], '-', color='#0cb5ed')
avgimg.plot([25, 55], [25, 25], '--', color='#0cb5ed')

right = diffimg.imshow(imgflux, cmap='YlOrRd')
diffimg.set_title('Difference')
right.set_interpolation('nearest')
diffimg.set_xlim(-0.5, x-0.5)
diffimg.set_ylim(y-0.5, -0.5)

right.axes.get_xaxis().set_ticklabels([])
right.axes.get_yaxis().set_ticklabels([])
right.axes.get_xaxis().set_ticks([])
right.axes.get_yaxis().set_ticks([])

if params.refpix == True:
   if eo == 1:
      diffimg.plot((x - refx) - 0.5, (y - refy) - 0.5, '*', color='#0cb5ed', ms=10)
      avgimg.plot((x - refx) - 0.5, (y - refy) - 0.5, '*', color='#0cb5ed', ms=10)
   elif eo == 0:
      diffimg.plot(refx - 1.5, (y - refy) - 0.5, '*', color='#0cb5ed', ms=10)
      avgimg.plot(refx - 1.5, (y - refy) - 0.5, '*', color='#0cb5ed', ms=10)

plt.tight_layout()
fig.set_size_inches(7.5, 4.5)
plt.savefig(f'kic{kic}q{q}imgs.png')

if params.show == True:
   plt.show()
else:
   pass
