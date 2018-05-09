from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
import sys
path = '/home/martin/Desktop/LEXI ADC project/atmospheric model/'
sys.path.insert(1,path)
import AtmosModel as AM
from astropy.io import fits
#test to mimic mikes images


#making a simple PSF of aperture.

Dtel=1 #(m)
GridSize = 100
pupil_grid = make_pupil_grid(GridSize)
aperture = circular_aperture(Dtel)
A = aperture(pupil_grid)
#imshow_field(A)
#plt.show()

#pha =  fits.getdata('phase.fits')
#plt.imshow(pha,origin='lower',interpolation='nearest')
#plt.show()
#pha =  np.hstack(pha)


# wavelenghts at in range 0.6 -0.8 um
wlen1 = 0.488E-6 #(m)
n =11
wlen2 = 0.633E-6 #(m)
wavelengths = np.linspace(wlen1,wlen2,n)


##--------------------------##
##--Simple PSF Single Wave--##
##--------------------------##

img_out0 = 0
l0 = wavelengths[5]
F = 150./1000. #(m)

focal_grid = make_focal_grid(pupil_grid,q=10,num_airy=10,focal_length=F, wavelength=l0)

p1 = 0.#np.sqrt(np.square(pupil_grid.x) + np.square(pupil_grid.y))
E1 =  A*np.exp(-1j*(p1))
E1 = np.nan_to_num(E1)

fraunhofer = FraunhoferPropagator(pupil_grid, focal_grid,focal_length=F, wavelength_0=l0)

wf1 = Wavefront(E1, wavelength=l0)
wf1.total_power = 1

img1 = fraunhofer(wf1)
img_out0 = img_out0 +img1.intensity 
imshow_field(np.log10(img_out0 / img_out0.max()), vmin=-6,interpolation='nearest')
plt.colorbar()
plt.show()
plt.close()


'''
##--------------------------##
##--Simple PSF Multi Wave---##
##--------------------------##

img_out1 = 0
l0 = 0.7 (#um)
F = 46.2 #(m)

focal_grid = make_focal_grid(pupil_grid,q=20,num_airy=10,focal_length=F, wavelength=l0*1E-6)
for w in wavelengths:

	p1 = np.sqrt(np.square(pupil_grid.x) + np.square(pupil_grid.y))
	E1 =  A*np.exp(-1j*(p1))
	E1 = np.nan_to_num(E1)
	#focal_grid = make_focal_grid(pupil_grid,q=20,num_airy=25,focal_length=F, wavelength=w)
	fraunhofer = FraunhoferPropagator(pupil_grid, focal_grid,focal_length=F, wavelength_0=l0*1E-6)

	wf1 = Wavefront(E1, wavelength=w)
	wf1.total_power = 1

	img1 = fraunhofer(wf1)
	img_out1 = img_out1 +img1.intensity 

#plt.subplot(1,2,1)
#imshow_field(np.log10(img_out0 / img_out0.max()), vmin=-6,interpolation='nearest')
#plt.colorbar()
#plt.subplot(1,2,2)
#imshow_field(np.log10(img_out1 / img_out1.max()), vmin=-6,interpolation='nearest')
#plt.colorbar()
#plt.show()
#plt.close()

##------------------------------##
##--Simple PSF added aberation--##
##------------------------------##

# The simulation needed to be in units of lambda/D for all angles! that is the only way to make it all work. 
img_out1 = 0
img_out2 = 0
l0 = 0.7
F = 42.6
phi = 45* np.pi/180
devl0 = AM.delR(AM.h0,phi,AM.a,AM.b,l0)
lambda_D = 0.14 /60 * np.pi/10800


focal_grid = make_focal_grid(pupil_grid,q=20,num_airy=10,focal_length=F, wavelength=l0*1E-6)
fraunhofer = FraunhoferPropagator(pupil_grid, focal_grid, focal_length=F, wavelength_0=l0*1E-6)
coronagraph = PerfectCoronagraph(A)

for w in wavelengths:
	refrac = AM.delR(AM.h0,phi,AM.a,AM.b,w*1E6)
	t = (devl0-refrac) / lambda_D
	p2 = 2*np.pi*(pupil_grid.y / Dtel)*t
	E1 =  A*np.exp(1j*p2)
	E1 = np.nan_to_num(E1)
	
	wf1 = Wavefront(E1, wavelength=w)
	wf1.total_power = 1

	wf2 = coronagraph.forward(wf1)
	imgc = fraunhofer(wf2)

	img1 = fraunhofer(wf1)
	img_out1 = img_out1 +img1.intensity 
	img_out2 = img_out2 +imgc.intensity 
	
plt.subplot(1,2,1)
imshow_field(np.log10(img_out1 / img_out1.max()), vmin=-6,interpolation='nearest')
plt.colorbar()
plt.subplot(1,2,2)
imshow_field(np.log10(img_out2 / img_out1.max()), vmin=-6,interpolation='nearest')
plt.colorbar()
plt.show()
plt.close()
'''
