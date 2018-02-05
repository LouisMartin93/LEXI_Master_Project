'''
Calculates the rotation for which the dispersion is minimised for a given zenith angle and saves as a fits file.
'''

#contains funcitons
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.patches import Rectangle
path = '/home/martin/Desktop/LEXI ADC project/atmospheric model/'
sys.path.insert(1,path)
import AtmosModel as AM
import PrismFunctions as PF
import time

Lambda = 0.7 #microns
re = 6371.0*10**3 # mean radius of the earth (m)
h0 = 2324 #altitude of observatory on la palma (m)
Mma = 0.0289644 # mean molar mass of dry air (kg/mol)
Lapse_rate = 6.49/1000. # lapse rate of  atmosphere (K/m)
g = 9.81 #acceleration due to gravity (m/s^2)
T0 = 4+273.15 # mean temperature at 00:00 in November at WHT (K)
R = 8.31447 # universal gas constant (J/(mol*K))
D = 1. # 4.26 #telescope diameter.

q = (g*Mma/(R*Lapse_rate))
H0 = T0/Lapse_rate
a = (q-1.)/H0
b = (q-1.)*(q-2.)/(2.*H0**2)
r0 = re + h0


#loading the prism data.
from astropy.io import ascii
table = ascii.read('glass tables/Glass_worksheet1.txt')
glass_by_shape = table.group_by('Prism Type')
prisms_same = glass_by_shape.groups[0]
limit = 1./10.* Lambda*10**-6 / D *180/np.pi *3600 *462#arcsec 
print 'limit is {} arcseconds'.format(np.round(limit,4))
rotate = np.linspace(110,180,601)
zenith = np.linspace(0,60,61)
prism1 = prisms_same[68]
prism2 = prisms_same[103]
l0= 0.7
bw =0.2
l1 = l0 - bw/2.
l2 = l0 + bw/2.
optimal_vals =[]
for z in zenith:
	disp4phi =[]
	dev_ref = AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l0)
	for phi in rotate:
		#for the smallest wavelenght
		offset1 = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l1)-dev_ref)*462
		ray1 = PF.CorrectionSingle(prism1,prism2,l1,phi,offset=offset1)
		
		#for the central wavelength
		ray_ref = PF.CorrectionSingle(prism1,prism2,l0,phi)
		
		#for the longest wavelength
		offset2 = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l2)-dev_ref)*462
		ray2 = PF.CorrectionSingle(prism1,prism2,l2,phi,offset=offset2)

		#calculating the dispersion
		dispersion = np.max([ray1,ray2,ray_ref])- np.min([ray1,ray2,ray_ref])
		disp4phi.append(dispersion*3600)
	
	disp4phi = np.array(disp4phi) 
	index, = np.where(disp4phi == np.min(disp4phi))
	optimal_vals.append(rotate[index])
	

import astropy.io.fits as fits
fits.writeto('OptRot4Zen.fits',np.array(zip(zenith,optimal_vals)),overwrite=True)
