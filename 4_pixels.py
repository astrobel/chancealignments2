import numpy as np
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
from astropy.stats import LombScargle
import smoothing 
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.pyplot as plt
from lightkurve import KeplerTargetPixelFile
import nancleaner as nc
import quarters as qs
from matplotlib.backends.backend_pdf import PdfPages as pdf
import argparse, warnings

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
parser.add_argument('-q', '--quarter', required=True, type=int, choices=range(0,18), help='Quarter to analyse')
parser.add_argument('-s', '--smoothing', dest='kern', default=100, type=int, help='Gaussian smoothing kernel, in days')
parser.add_argument('-c', '--clip', dest='inp', default=3, type=int, help='Outlier clipping level, in sigma')
parser.add_argument('-o', '--oversampling', dest='over', default=5, type=int, help='LSP oversampling factor')
parser.add_argument('-n', '--nyquistfactor', dest='nyq', default=1, type=float, help='LSP Nyquist factor')
parser.add_argument('-m', '--makepdf', default=False, type=bool, help='Make PDF of pixel close-ups?')
parser.add_argument('-e', '--export', default=False, type=bool, help='Export data for each pixel?')
parser.add_argument('-p', '--plots', dest='show', default=False, type=bool, help='Show plots?')

params = parser.parse_args()

q = params.quarter
kic = params.kic

# flag for quarter of choice

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

if (channel%2) == 0:
   eo = 0
else:
   eo = 1

# dynamic variable names
for (j, k), img in np.ndenumerate(temp2d):
   index = (k + 1) * j + (x - j) * k
   exec("pixel%d_flux = np.array(None)" % index)
   exec("pixel%d_time = np.array(None)" % index)

# filling the flux array
second_flux = np.zeros([xsize, ysize])
for (i, j, k), val in np.ndenumerate(flux):
   index = (j + 1) * k + (x - k) * j
   second_flux[index, i] = val

for (j, k), img in np.ndenumerate(table2):
   index = (j + 1) * k + (x - k) * j
   if img == 0:
      pass
   else:
      tempflux1, temptime = nc.nancleaner2d(second_flux[index,:], time)
      tempflux, smth_flux = smoothing.gausssmooth(temptime, tempflux1, params.kern)

      clip = params.inp * np.std(tempflux)
      meanflux = np.mean(tempflux)
 
      upperbound = meanflux + clip
      lowerbound = meanflux - clip

      colours = np.zeros(tempflux.size)

      for i, flux in enumerate(tempflux):
         if flux < upperbound and flux > lowerbound:
            colours[i] = 1

      clipped_flux = []
      clipped_time = []
      for i, colour in enumerate(colours):
         if colour == 1:
            clipped_flux.append(tempflux[i])
            clipped_time.append(temptime[i])

      exec("pixel%d_flux = clipped_flux" % index)
      exec("pixel%d_time = clipped_time" % index)
 
      # export smoothed and clipped data as .dat file
      if params.export == True:
         exportblend = np.array([clipped_time, clipped_flux])
         exportblend = np.transpose(exportblend)
         np.savetxt(f'kic{kic}_pixel{index+1}_lc.dat', exportblend, delimiter=' ', header=f'Smoothed and clipped light curve for KIC{kic} TPF')
      else:
         pass

      # fourier transform

      ofac = params.over
      hifac = params.nyq

      frequencies, power_spectrum = LombScargle(np.asarray(clipped_time), np.asarray(clipped_flux)).autopower(method='fast', normalization='psd', samples_per_peak=ofac, nyquist_factor=hifac)
      hifac *= (283/11.57)/max(frequencies)
      frequencies, power_spectrum = LombScargle(np.asarray(clipped_time), np.asarray(clipped_flux)).autopower(method='fast', normalization='psd', samples_per_peak=ofac, nyquist_factor=hifac)
      power_spectrum = power_spectrum * 4. / len(clipped_time)
      power_spectrum = np.sqrt(power_spectrum)
      power_spectrum *= 1e6
      frequencies *= 11.57

      exec("pixel%d_freq = frequencies" % index)
      exec("pixel%d_ps = power_spectrum" % index)


### PLOTTING ###

# light curves
fig = plt.figure(1)
gs = gridspec.GridSpec(y, x, wspace=0, hspace=0)
plt.title(f'{kic}')
plt.xlabel('Time (d)')
plt.ylabel('Fractional Intensity')

ax1 = plt.gca()
ax1.get_xaxis().set_ticks([])
ax1.get_yaxis().set_ticks([])

for (j, k), img in np.ndenumerate(table2):
   index = (j + 1) * k + (x - k) * j
   if img == 0:
      if eo == 0:
         ax = fig.add_subplot(gs[y - j - 1, x - k - 1])
      elif eo == 1:
         ax = fig.add_subplot(gs[y - j - 1, k])
      ax.set_xticklabels('')
      ax.set_yticklabels('')
   else:
      exec("flux = pixel%d_flux" % index)
      exec("time = pixel%d_time" % index)
      if eo == 0:
         ax = fig.add_subplot(gs[y - j - 1, x - k - 1])
      elif eo == 1:
         ax = fig.add_subplot(gs[y - j - 1, k])
      ax.set_xticklabels('')
      ax.set_yticklabels('')
      if img == np.amax(table2):
         lower = np.zeros(len(time))
         upper = lower + max(flux)
         ax.fill_between(time, lower, upper, facecolor='#ff99a3')
         plt.plot(time, flux, 'k.', ms=0.1)
         plt.ylim(min(flux), max(flux)) #ymin=0)
         plt.xlim(min(time), max(time))
      else:
         plt.plot(time, flux, 'k.', ms=0.1)
         plt.ylim(min(flux), max(flux)) #ymin=0)
         plt.xlim(min(time), max(time))

plt.savefig(f'kic{kic}q{q}pixelslc.png')

if params.makepdf == True:
   outplot = pdf(f'kic{kic}_q{q}.pdf')
else:
   pass

# power spectra
fig = plt.figure(2)
gs = gridspec.GridSpec(y, x, wspace=0, hspace=0)
plt.title(f'{kic}')
plt.xlabel('Frequency ($\mu$Hz)')
plt.ylabel('Amplitude (ppm)')

ax0 = plt.gca()
ax0.get_xaxis().set_ticks([])
ax0.get_yaxis().set_ticks([])

for (j, k), img in np.ndenumerate(table2):
   index = (j + 1) * k + (x - k) * j
   if img == 0:
      if eo == 0:
         ax = fig.add_subplot(gs[y - j - 1, x - k - 1])
      elif eo == 1:
         ax = fig.add_subplot(gs[y - j - 1, k])
      ax.set_xticklabels('')
      ax.set_yticklabels('')
   else:
      exec("freq = pixel%d_freq" % index)
      exec("ps = pixel%d_ps" % index)
      if eo == 0:
         ax = fig.add_subplot(gs[y - j - 1, x - k - 1])
      elif eo == 1:
         ax = fig.add_subplot(gs[y - j - 1, k])
      ax.set_xticklabels('')
      ax.set_yticklabels('')
      if img == np.amax(table2):
         lower = np.zeros(len(freq))
         upper = lower + max(ps)
         ax.fill_between(freq, lower, upper, facecolor='#ff99a3')
         plt.plot(freq, ps, 'k-', lw=0.5)
         plt.ylim(0, max(ps)) #ymin=0)
         plt.xlim(0, max(freq))
      else:
         plt.plot(freq, ps, 'k-', lw=0.5)
         plt.ylim(0, max(ps)) #ymin=0)
         plt.xlim(0, max(freq))

fig.set_size_inches(14,10)
plt.savefig(f'kic{kic}q{q}pixels.png')

if params.makepdf == True:
   # power spectra 2
   fig = plt.figure(3)
   gs = gridspec.GridSpec(y, x, wspace=0, hspace=0)
   plt.title(f'{kic}')
   plt.xlabel('Frequency ($\mu$Hz)')
   plt.ylabel('Amplitude (ppm)')

   ax0 = plt.gca()
   ax0.get_xaxis().set_ticks([])
   ax0.get_yaxis().set_ticks([])

   for (j, k), img in np.ndenumerate(table2):
      index = (j + 1) * k + (x - k) * j
      if img == 0:
         if eo == 0:
            ax = fig.add_subplot(gs[y - j - 1, x - k - 1])
         elif eo == 1:
            ax = fig.add_subplot(gs[y - j - 1, k])
         ax.set_xticklabels('')
         ax.set_yticklabels('')
      else:
         exec("freq = pixel%d_freq" % index)
         exec("ps = pixel%d_ps" % index)
         if eo == 0:
            ax = fig.add_subplot(gs[y - j - 1, x - k - 1])
         elif eo == 1:
            ax = fig.add_subplot(gs[y - j - 1, k])
         ax.set_xticklabels('')
         ax.set_yticklabels('')
         if img == np.amax(table2):
            lower = np.zeros(len(freq))
            upper = lower + max(ps)
            ax.fill_between(freq, lower, upper, facecolor='#ff99a3')
            plt.plot(freq, ps, 'k-', lw=0.5)
            plt.annotate(index,xy=(max(freq)/2, max(ps)/2))
            plt.ylim(0, max(ps)) #ymin=0)
            plt.xlim(0, max(freq))
         else:
            plt.plot(freq, ps, 'k-', lw=0.5)
            plt.annotate(index,xy=(max(freq)/2, max(ps)/2))
            plt.ylim(0, max(ps)) #ymin=0)
            plt.xlim(0, max(freq))
   outplot.savefig(fig)

   # pixels
   for (j, k), img in np.ndenumerate(table2):
      index = (j + 1) * k + (x - k) * j
      if img == 0:
         pass
      else:
         exec("flux = pixel%d_flux" % index)
         exec("time = pixel%d_time" % index)
         exec("freq = pixel%d_freq" % index)
         exec("ps = pixel%d_ps" % index)
         if img == np.amax(table2):
            fig, ax = plt.subplots(2,1)
            plt.title(index)
            ax[0].plot(time, flux, 'k.', ms=1)
            ax[0].set_ylim(min(flux), max(flux)) #ymin=0)
            ax[0].set_xlim(min(time), max(time))
            ax[1].plot(freq, ps, 'k-', lw=0.5)
            ax[1].set_ylim(0, max(ps)) #ymin=0)
            ax[1].set_xlim(0, max(freq))
            ax[0].set_xlabel('time (d)')
            ax[0].set_ylabel('fractional intensity')
            ax[1].set_xlabel('freq ($\mu$Hz)')
            ax[1].set_ylabel('amplitude (ppm)')
            plt.tight_layout()
            outplot.savefig(fig)
         else:
            fig, ax = plt.subplots(2,1)
            plt.title(index)
            ax[0].plot(time, flux, 'k.', ms=1, alpha=0.5)
            ax[0].set_ylim(min(flux), max(flux)) #ymin=0)
            ax[0].set_xlim(min(time), max(time))
            ax[1].plot(freq, ps, 'k-', lw=0.5, alpha=0.5)
            ax[1].set_ylim(0, max(ps)) #ymin=0)
            ax[1].set_xlim(0, max(freq))
            ax[0].set_xlabel('time (d)')
            ax[0].set_ylabel('fractional intensity')
            ax[1].set_xlabel('freq ($\mu$Hz)')
            ax[1].set_ylabel('amplitude (ppm)')
            plt.tight_layout()
            outplot.savefig(fig)
   outplot.close()
else:
   pass

if params.show == True:
   plt.show()
else:
   pass