import numpy as np
from astropy.io import fits as pyfits
from astropy import wcs
from astropy.utils.exceptions import AstropyWarning
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from lightkurve import KeplerTargetPixelFile
import quarters as qs
import argparse, os, warnings

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=UserWarning) # for font conflicts on my system, at least

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
      r'\usepackage{helvet}',
      r'\usepackage[EULERGREEK]{sansmath}',
      r'\sansmath'
]

# to add: flag for filename of ukirt image; flag for quarter of choice

parser = argparse.ArgumentParser(description='Prepare long cadence light curve for further analysis.')
parser.add_argument('-k', '--kic', required=True, type=int, help='KIC ID')
parser.add_argument('-q', '--quarter', required=True, type=int, choices=range(1,18), help='Quarter to analyse')
parser.add_argument('-u', '--ukirt', required=True, type=str, help='UKIRT image file string, including fits extension')
parser.add_argument('-r', '--refpix', default=False, type=bool, help='Plot location of reference pixel on image?')
parser.add_argument('-p', '--plots', dest='show', default=False, type=bool, help='Show plots?')

params = parser.parse_args()

q = params.quarter
kic = params.kic

### IMAGE 1: ONE KEPLER PIXEL IMAGE, Q_ ###

tpf = KeplerTargetPixelFile.from_archive(kic, quarter=q)

channel = tpf.channel
obj_ra = tpf.ra
obj_dec = tpf.dec
table = tpf.hdu[1].data
flux = table['FLUX']
time = table['TIME']
hd1 = tpf.hdu[1].header
ysize = hd1['NAXIS2']
table2 = tpf.hdu[2].data
hd2 = tpf.hdu[2].header
x = hd2['NAXIS1']
y = hd2['NAXIS2']
xsize = x * y
temp2d = np.zeros((x, y))
w = wcs.WCS(hd1, keysel=['binary'])

if (channel%2) == 0:
   eo = 0
else:
   eo = 1

imgflux = np.flipud(flux[0])
if eo == 0:
   imgflux = np.fliplr(imgflux)


### IMAGE 2: UKIRT IMAGE ###

hdulist = pyfits.open(params.ukirt)

flux2 = hdulist[1].data
parameters = hdulist[1].header
cam = parameters['CAMNUM']

rotby = 360 - 90*cam
rotby /= 90
flux2 = np.rot90(flux2, rotby)


### PLOTTING ###

fig, (kepler, ukirt) = plt.subplots(1, 2) 

fig.suptitle(f'{kic}', fontsize=20)

left = kepler.imshow(imgflux, cmap='YlOrRd')
kepler.set_title('Kepler Aperture')
left.set_interpolation('nearest')
kepler.set_xlim(-0.5, x-0.5)
kepler.set_ylim(y-0.5, -0.5)

left.axes.get_xaxis().set_ticklabels([])
left.axes.get_yaxis().set_ticklabels([])
left.axes.get_xaxis().set_ticks([])
left.axes.get_yaxis().set_ticks([])

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
kepler.plot(npixels[:,0], npixels[:,1], color='#00ff8c')
kepler.plot(epixels[:,0], epixels[:,1], '--', color='#00ff8c')

if params.refpix == True:
   print(' ')
   if eo == 1:
      print('--> Reference pixel:')
      print(f'-----> RA = {refra}')
      print(f'-----> DEC = {refdec}')
      kepler.plot((x - refx) - 0.5, (y - refy) - 0.5, '*', color='#00ff8c', ms=10)
   elif eo == 0:
      print('--> Reference pixel:')
      print(f'-----> RA = {refra}')
      print(f'-----> DEC = {refdec}')
      kepler.plot(refx - 1.5, (y - refy) - 0.5, '*', color='#00ff8c', ms=10)
else:
   pass

right = ukirt.imshow(flux2, cmap='YlOrRd', norm=LogNorm())
ukirt.set_title('UKIRT Image')
right.set_interpolation('bilinear')
ukirt.set_xlim(300, 0)
ukirt.set_ylim(0, 300)

right.axes.get_xaxis().set_ticklabels([])
right.axes.get_yaxis().set_ticklabels([])
right.axes.get_xaxis().set_ticks([])
right.axes.get_yaxis().set_ticks([])

ukirt.plot([25, 25], [25, 55], '-', color='#00ff8c')
ukirt.plot([25, 55], [25, 25], '--', color='#00ff8c')

fig.set_size_inches(7.5, 4.5)
plt.savefig(f'kic{kic}img.png')

if params.show == True:
   plt.show()
else:
   pass
