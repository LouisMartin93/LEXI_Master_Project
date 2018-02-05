'''
Does as it says on the tin..

This will take an input of prisms to use and will then plot out what the dispersion is for a given rotation at a given zenith
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
rotate = np.linspace(0,180,181)
zenith = np.linspace(0,60,7)

prism1 = prisms_same[63]
prism2 = prisms_same[92]
l0= 0.7
bw =0.2
l1 = l0 - bw/2.
l2 = l0 + bw/2.
print 'Combination of {} and {} glass with wedge angles of {} and {}'.format(prism1['Material'],prism2['Material'],prism1['Wedge angle'],prism2['Wedge angle'])
zero_cross=[]
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
	plt.plot(rotate,disp4phi,label=str(z))
	disp4phi = disp4phi - limit
	sign_cross = np.diff(np.sign(disp4phi))
	
	if 2 in np.abs(sign_cross):
		zero_cross.append(1)
	else:
		zero_cross.append(0)
	#print zero_cross
	

if (len(np.unique(zero_cross)) == 1) and (np.unique(zero_cross) == 1):
	print 'good combination'
plt.plot([0,180],[limit,limit],'--k',label='lim')
plt.legend(loc ='best')
plt.xlabel('Rotation of prisms with respect to one another (deg)')
plt.ylabel('Dispersion (arcseconds)')
plt.title('Combination of {} and {} glass with wedge angles of {} and {}'.format(prism1['Material'],prism2['Material'],prism1['Wedge angle'],prism2['Wedge angle']))
plt.show()


