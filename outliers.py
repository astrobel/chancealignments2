import numpy as np

def clip(time, flux, bounds):

   sigma = np.std(flux)
   avg = np.mean(flux)
   outliers = []

   for i, val in enumerate(flux):
      if val > avg + sigma*bounds or val < avg - sigma*bounds:
         outliers.append(i)

   fluxnew = np.delete(flux, outliers)
   timenew = np.delete(time, outliers)

   return timenew, fluxnew