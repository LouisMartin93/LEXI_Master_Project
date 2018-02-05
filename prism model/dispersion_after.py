'''
This calculates what prism combiniations will meet the criteria of correcting dispersion to the required minimum. the combinations are then printed (though should be saved as a fits file..)
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



theta = np.linspace(180,360,181)




#loading the prism data.
from astropy.io import ascii
table = ascii.read('glass tables/Glass_worksheet1.txt')
glass_by_shape = table.group_by('Prism Type')
prisms_same = glass_by_shape.groups[0]
#print prisms_same
combo = []
limit = 1./10.* Lambda*10**-6 / D *180/np.pi *3600 *462.#arcsec 
print 'limit is {} arcseconds'.format(np.round(limit,4))
rotate = np.linspace(0,180,181)
zenith = np.linspace(0,60,7)

for i in range(len(prisms_same)):
	for j in range(len(prisms_same)):
		#start = time.time()
		prism2 = prisms_same[j]
		prism1 = prisms_same[i]

		

		
		if prism1['Wedge angle'] > prism2['Wedge angle']:
			continue
		#if prism1['nd'] > prism['nd']:
			#break
		print 'Combination of {} and {} glass with wedge angles of {} and {}'.format(prism1['Material'],prism2['Material'],prism1['Wedge angle'],prism2['Wedge angle'])
		l0= 0.7
		bw =0.2


		l1 = l0 - bw/2.
		l2 = l0 + bw/2.

		
		zero_cross =[]
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
			#check to see if crosses 0 at any point.
			disp4phi = np.array(disp4phi) 
			sign_cross = np.diff(np.sign(disp4phi-limit))
			
			if 2 in np.abs(sign_cross):
				zero_cross.append(1)
			else:
				zero_cross.append(0)
		#print zero_cross
		if (len(np.unique(zero_cross)) == 1) and (np.unique(zero_cross) == 1):
			combo.append([i,j])
			print 'good combo'
			print i,j
		#print time.time() - start #time is 6.09s for 2775 combos gives total run time of 4.625 hrs
		

print combo	

"""
prism1 = prisms_same[46]
prism2 = prisms_same[2]
l0= 0.7
bw =0.2
l1 = l0 - bw/2.
l2 = l0 + bw/2.
print 'Combination of {} and {} glass with wedge angles of {} and {}'.format(prism1['Material'],prism2['Material'],prism1['Wedge angle'],prism2['Wedge angle'])
zero_cross=[]
for z in zenith:
	disp4phi =[]
	
	dev_ref = AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l0)
	offset2 = (dev_ref- AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l2))*462
	print offset2 *180./np.pi *3600
	offset1 = (dev_ref- AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l1))*462
	print offset1 *180./np.pi *3600
	for phi in rotate:
		ray1 = PF.CorrectionSingle(prism1,prism2,0.6,phi,offset=offset1)
		ray2 = PF.CorrectionSingle(prism1,prism2,0.8,phi,offset=offset2)
		disp4phi.append((ray1-ray2)*3600)
		
	disp4phi = np.array(disp4phi)[:,0,0]
	plt.plot(rotate,disp4phi,label=str(z))
	sign_cross = np.diff(np.sign(disp4phi))
	
	if 2 in np.abs(sign_cross):
		zero_cross.append(1)
	else:
		zero_cross.append(0)
	print zero_cross
plt.plot([180,360],[0,0],'--k',label='lim')
plt.legend(loc ='best')
plt.show()
"""

