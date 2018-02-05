'''

This script takes the all the possiblities of prism combinations where they met the requirements from ray_trace.py and the calculates the angle of rotation spread wrt as the prism combo corrects from 10 to 60 degrees from the zenith. It then prints the top 10 combinations

An adaptation of this would be to also include the criteria such that off the shelf prisms must be included in at least one of the combinations, and then the top 10 be printed.

'''


import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import sys
path = '/home/martin/Desktop/LEXI ADC project/atmospheric model/'
sys.path.insert(1,path)
import AtmosModel as AM
import PrismFunctions as PF
import time
from astropy.io import ascii


# get combo data
combos = fits.getdata('Combo_list.fits')

Lambda = 0.7 # in microns
limit = 1./10.* Lambda*10**-6 / 1. *180/np.pi *3600 *462.#arcsec 
rotate = np.linspace(90,180,(10*180)+1)

z1 = 60 #degrees
z2 = 10 # degrees

#loading the prism data.
table = ascii.read('glass tables/Glass_worksheet1.txt')
glass_by_shape = table.group_by('Prism Type')
prism_same = glass_by_shape.groups[0]

l0= 0.7
bw =0.2
l1 = l0 - bw/2.
l2 = l0 + bw/2.
RotSpread = {}

# for z = 60
dev_ref_60 = AM.delR(AM.h0,z1*np.pi/180.,AM.a,AM.b,l0)
offset1_60 = (AM.delR(AM.h0,z1*np.pi/180.,AM.a,AM.b,l1)-dev_ref_60)*462
offset2_60 = (AM.delR(AM.h0,z1*np.pi/180.,AM.a,AM.b,l2)-dev_ref_60)*462
# for z = 10
dev_ref_10 = AM.delR(AM.h0,z2*np.pi/180.,AM.a,AM.b,l0)
offset1_10 = (AM.delR(AM.h0,z2*np.pi/180.,AM.a,AM.b,l1)-dev_ref_10)*462
offset2_10 = (AM.delR(AM.h0,z2*np.pi/180.,AM.a,AM.b,l2)-dev_ref_10)*462

for c in combos:
	disp4phi60 =[]
	disp4phi10 =[]
	prism1 = prism_same[c[0]]
	prism2 = prism_same[c[1]]
	
	for phi in rotate:
		
		#for z =60
		#for the smallest wavelenght
		ray1 = PF.CorrectionSingle(prism1,prism2,l1,phi,offset=offset1_60)
		#for the central wavelength
		ray_ref = PF.CorrectionSingle(prism1,prism2,l0,phi)
		#for the longest wavelength
		ray2 = PF.CorrectionSingle(prism1,prism2,l2,phi,offset=offset2_60)
		#calculating the dispersion
		dispersion = np.max([ray1,ray2,ray_ref])- np.min([ray1,ray2,ray_ref])
		disp4phi60.append(dispersion*3600) #arcsec
		
		#for z = 10
		#for the smallest wavelenght
		ray1 = PF.CorrectionSingle(prism1,prism2,l1,phi,offset=offset1_10)
		#for the central wavelength
		ray_ref = PF.CorrectionSingle(prism1,prism2,l0,phi)
		#for the longest wavelength
		ray2 = PF.CorrectionSingle(prism1,prism2,l2,phi,offset=offset2_10)
		#calculating the dispersion
		dispersion = np.max([ray1,ray2,ray_ref])- np.min([ray1,ray2,ray_ref])
		disp4phi10.append(dispersion*3600) #arcsec
	
	disp4phi10 = np.array(disp4phi10)
	#plt.plot(rotate,disp4phi10,label='10')
	disp4phi60 = np.array(disp4phi60)
	#plt.plot(rotate,disp4phi60,label='60')
	min_index10 = np.argmin(disp4phi10) +1
	disp4phi10 = disp4phi10[min_index10:] - limit
	new_min10 = np.argmin(abs(disp4phi10))+min_index10
	#print rotate[new_min10]
	min_index60 = np.argmin(disp4phi60) +1
	disp4phi60 = disp4phi60[:min_index60] - limit
	new_min60 = np.argmin(abs(disp4phi60))
	#print rotate[new_min60]
	#plt.plot([90,180],[limit,limit],'--k')
	#plt.legend(loc='best')
	#plt.show()
	max_spread = rotate[new_min10] - rotate[new_min60]
	
	RotSpread[max_spread] = (c[0],c[1])
print RotSpread 

top10 = sorted(RotSpread.keys(),reverse=True)[:10]
print top10
# might be some combinations overlap due to accuracy of decimal places.
for spread in top10:
	c = RotSpread[spread]
	prism1 = prism_same[c[0]]
	prism2 = prism_same[c[1]]
	print 'Combination of {} and {}  with wedge angles of {} and {} gives spread of {} '.format(prism1['Material'],prism2['Material'],prism1['Wedge angle'],prism2['Wedge angle'],np.round(spread,3))

