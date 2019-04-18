import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy import wcs
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
from astropy.utils.exceptions import AstropyWarning
from astropy.stats import LombScargle
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.gaia import Gaia
from lightkurve import search_targetpixelfile, LightkurveWarning
import nancleaner as nc
import argparse, warnings

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=UserWarning) # for font conflicts on my system, at least
warnings.simplefilter('ignore', category=RuntimeWarning) # may encounter some divide by zeros in process
warnings.simplefilter('ignore', category=LightkurveWarning) # if it's trying to download an empty quarter, code will quit

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

parser = argparse.ArgumentParser(description='Map out the region using Gaia DR2.')
parser.add_argument('-k', '--kic', required=True, type=int, help='KIC ID')
parser.add_argument('-q', '--quarter', required=True, type=int, choices=range(0,18), help='Quarter to analyse')
parser.add_argument('-s', '--sourceprint', default=False, type=bool, help='Print all Gaia DR2 sources?')
parser.add_argument('-p', '--plots', dest='show', default=False, type=bool, help='Show plots?')

params = parser.parse_args()

q = params.quarter
kic = params.kic

while True:
   tpf = search_targetpixelfile(kic, quarter=q).download()
   if tpf == None:
      print('No data for this quarter.')
      sys.exit()
   else:
      break

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
avgflux = np.flipud(avgflux)
if eo == 0:
   avgflux = np.fliplr(avgflux)

# find stars in the neighbourhood

target = SkyCoord(ra=obj_ra, dec=obj_dec, unit=(u.degree, u.degree), frame='icrs')
neighbours = Gaia.query_object_async(coordinate=target, width=u.Quantity(0.001*x, u.deg), height=u.Quantity(0.001*y, u.deg))
ras = neighbours['ra']
decs = neighbours['dec']
sources = neighbours['source_id']
mags = neighbours['phot_g_mean_mag']


### PLOTTING ###

fig, ax = plt.subplots(1) 

crval = w.wcs.crval
north = crval + np.array([0, 6/3600.])
east = crval + np.array([ 6/3600., 0])

ncoords = np.vstack([crval, north])
ecoords = np.vstack([crval, east])
npixels = w.wcs_world2pix(ncoords, 0)
epixels = w.wcs_world2pix(ecoords, 0)
npixels[1, 1] = npixels[0, 1] - (npixels[1, 1] - npixels[0, 1]) # flip ud
epixels[1, 1] = epixels[0, 1] - (epixels[1, 1] - epixels[0, 1])
if eo == 0:
   npixels[1, 0] = npixels[0, 0] - (npixels[1, 0] - npixels[0, 0]) # flip lr
   epixels[1, 0] = epixels[0, 0] - (epixels[1, 0] - epixels[0, 0])

left = ax.imshow(avgflux, cmap='YlOrRd')
left.set_interpolation('nearest')
ax.set_xlim(-0.5, x-0.5)
ax.set_ylim(y-0.5, -0.5)
ax.plot(npixels[:,0], npixels[:,1], color='#0cb5ed')
ax.plot(epixels[:,0], epixels[:,1], '--', color='#0cb5ed')

left.axes.get_xaxis().set_ticklabels([])
left.axes.get_yaxis().set_ticklabels([])
left.axes.get_xaxis().set_ticks([])
left.axes.get_yaxis().set_ticks([])

ax.plot([25, 25], [25, 55], '-', color='#0cb5ed')
ax.plot([25, 55], [25, 25], '--', color='#0cb5ed')

# https://stackoverflow.com/questions/17285163/clean-way-to-use-words-as-markers-in-matplotlib-and-make-font-size-and-color
for i in range(len(ras)):
   coords = w.wcs_world2pix(ras[i], decs[i], 0)
   if eo == 1:
      ax.scatter(coords[0], y - coords[1] - 1, c='k', linewidth=0.5, marker=r'$ {:.2f} $'.format(mags[i]), s=300, label=mags[i])
   elif eo == 0:
      ax.scatter(x - coords[0] - 1, y - coords[1] - 1, c='k', linewidth=0.5, marker=r'$ {:.2f} $'.format(mags[i]), s=300, label=mags[i])

   if params.sourceprint == True:
      print(f'Source mag {mags[i]:.2f}: {sources[i]}')

# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width*0.6, box.height])
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.savefig(f'kic{kic}q{q}map.png')

if params.show == True:
   plt.show()
else:
   pass
