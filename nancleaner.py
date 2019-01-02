import numpy as np

def nancleaner3d(flux, time):
   dim = len(flux)
   nans = []

   for i in range(dim):
      if np.isnan(flux[i,:]).all() == True:
         nans.append(i)

   fluxnew = np.delete(flux, nans, axis=0)
   timenew = np.delete(time, nans)
   return fluxnew, timenew

def nancleaner2d(flux, time):
   blend = np.array([time, flux])
   blend = np.transpose(blend)
   blend2 = np.ma.compress_rows(np.ma.fix_invalid(blend))

   timenew = blend2[:,0]
   fluxnew = blend2[:,1]
   return fluxnew, timenew