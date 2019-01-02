import numpy as np
from astropy.io import fits as pyfits
from astropy import wcs
from lightkurve import KeplerTargetPixelFile, KeplerLightCurveFile
import os

def getallquarters(kic):

   # pick a quarter
   quarterlist = []

   for q in np.arange(0,18):
      try:
         tpf = KeplerLightCurveFile.from_archive(kic, quarter=q)
         quarterlist.append(q)
      except:
         continue

   return quarterlist

def getquarter(kic):
   
   quarterlist = getallquarters(kic)

   repquarter = 18 # placeholder, since there is 0 but there is no 18
   while repquarter not in quarterlist:
      while True:
         try:
            repquarter = int(input('--> Which quarter? (0-17) '))
            break
         except ValueError:
            print('--> Please enter an integer')
      if repquarter in quarterlist:
         break
      else:
         print('--> No data for this quarter')
         continue
   print(' ')

   return repquarter