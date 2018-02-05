'''
This is the test script before the library of functions were made (see OffsetModel.py) to calcualte the offset of a given prism model (Deviation from optic axis) as function of zenith as well as the function to calculate the rotation needed for a given zenith. 
'''


#contains funcitons
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.patches import Rectangle
path = '/home/martin/Desktop/LEXI ADC project/atmospheric model/'
sys.path.insert(1,path)
path = '/home/martin/Desktop/LEXI ADC project/prism model/'
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
table = ascii.read('../glass tables/Glass_worksheet1.txt')
glass_by_shape = table.group_by('Prism Type')
prisms_same = glass_by_shape.groups[0]
limit = 1./10.* Lambda*10**-6 / D *180/np.pi *3600 *462#arcsec 
print 'limit is {} arcseconds'.format(np.round(limit,4))

prism1 = prisms_same[68]
prism2 = prisms_same[103]
l0= 0.7
bw =0.2
l1 = l0 - bw/2.
l2 = l0 + bw/2.
print 'Combination of {} and {} glass with wedge angles of {} and {}'.format(prism1['Material'],prism2['Material'],prism1['Wedge angle'],prism2['Wedge angle'])

import astropy.io.fits as fits
data = fits.getdata('../OptRot4Zen.fits')
#data = fits.getdata('test1.fits')
rotate = data[:,1]
zenith = data[:,0]

offsets =[] #note the pural difference
zero_mat = np.matrix([0.0])

for z, phi in zip(zenith,rotate):
    dev_ref = AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l0)
    # for the central wavelength
    ray = PF.CorrectionSingle(prism1,prism2,l0,phi)
    offset = np.max([ray,zero_mat])
    offsets.append(offset) # offset in degrees

# assuming ray is initally centered along axis.
offsets = np.array(offsets)
plt.plot(zenith,offsets) 
plt.xlabel('Zenith (deg)')
plt.ylabel('Offset angle (deg)')
plt.title('Offset for rotation angle needed at given zenith')
plt.show()

plt.plot(zenith,(180-rotate))
plt.xlabel('Zenith (deg)')
plt.ylabel('rotation angle required (deg)')
plt.title('Optimal rotation angle needed for a given zenith (180 is zeroth position)')
plt.show()

# making my fit model for this particular combination
index, = np.where(rotate != 180)
ind = index[0] # so there are no zero values
y = 180 - rotate[ind:]
x =  zenith[ind:]
co_effs = np.polyfit(x, y, 8)
p = np.poly1d(co_effs)
print p(x)
xp = np.linspace(0, 60, 100)
plt.plot(x, y, '.', label='measured')
plt.plot(xp, p(xp), '-',label='model')
plt.xlabel('Zenith (deg)')
plt.ylabel('Rotation needed (deg)')
plt.show()

diff =[]
for i in range(len(x)):
    #print p(x[i]),y[i], p(x[i])-y[i]
    diff.append(p(x[i])-y[i])

abs_diff = np.abs(diff)
print np.mean(abs_diff)
print np.std(abs_diff)
#looking at the dispersion after the prisms correct for the atmosphere (measured values)
disp4phi =[]
for z, phi in zip(zenith,rotate):
    dev_ref = AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l0)
    #for the smallest wavelenght
    offset1 = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l1)-dev_ref)*462
    ray1 = PF.CorrectionSingle(prism1,prism2,l1,phi,offset=offset1)
		    
    #for the central wavelength
    ray_ref = PF.CorrectionSingle(prism1,prism2,l0,phi)
		    
    #for the longest wavelength
    offset2 = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l2)-dev_ref)*462
    ray2 = PF.CorrectionSingle(prism1,prism2,l2,phi,offset=offset2)

    #calculating the dispersion
    dispersion = np.max([ray1,ray2,ray_ref]) - np.min([ray1,ray2,ray_ref])
    disp4phi.append(dispersion*3600)
	
disp4phi = np.array(disp4phi) 
plt.plot(zenith,disp4phi,label='Recorded values')

#looking at the dispersion after the prisms correct for the atmosphere (model values)
disp4phi =[]
for z in zenith:
    dev_ref = AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l0)
    #for the smallest wavelenght
    offset1 = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l1)-dev_ref)*462
    ray1 = PF.CorrectionSingle(prism1,prism2,l1,180 - p(int(z)),offset=offset1)
		    
    #for the central wavelength
    ray_ref = PF.CorrectionSingle(prism1,prism2,l0,180 - p(int(z)))
		    
    #for the longest wavelength
    offset2 = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l2)-dev_ref)*462
    ray2 = PF.CorrectionSingle(prism1,prism2,l2,180 - p(int(z)),offset=offset2)

    #calculating the dispersion
    dispersion = np.max([ray1,ray2,ray_ref])- np.min([ray1,ray2,ray_ref])
    disp4phi.append(dispersion*3600)

disp4phi = np.array(disp4phi) 
plt.plot(zenith,disp4phi,label='Predicted Values')
plt.legend(loc='best')
plt.xlabel('Zenith angle')
plt.ylabel('Dispersion after correction (arcseconds)')
plt.title('Dispersion after correction at given Zenith angles using model')
plt.show()



fname = 'test1'
prisms=(68,103)
Ord = 8

from OffsetModel import AxisOffsetAndModel

#poly_cos = AxisOffsetAndModel(prisms,fname,Ord)
