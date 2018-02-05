from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
import sys
path ='/data2/Louis/LEXI ADC project/atmospheric model/'
sys.path.insert(1,path)
import AtmosModel as AM
from astropy.io import fits
#test to mimic mikes images

Dtel=1
GridSize =900
pupil_grid = make_pupil_grid(GridSize)
aperture = circular_aperture(Dtel)
A = aperture(pupil_grid)
#imshow_field(A)
#plt.show()

pha =  fits.getdata('cMWS/phase.fits')
#plt.imshow(pha,origin='lower',interpolation='nearest')
#plt.show()

pha =  np.hstack(pha)
#print pha.shape
wlen1 = 0.6E-6
n =11
wlen2 = 0.8E-6 

# wavelenghts at in range 0.6 -0.8 um
wavelenghts = np.linspace(wlen1,wlen2,n)
img_out0 = 0
l0 = 0.7
phi = 0.0* np.pi/180

devl0 = AM.delR(AM.h0,phi,AM.a,AM.b,l0)

for w in wavelenghts:
	refrac = AM.delR(AM.h0,phi,AM.a,AM.b,w*1E6)
	t = devl0-refrac
	p1 = np.sqrt(np.square(pupil_grid.x) + np.square(pupil_grid.y))/pha
	p2 = (np.sqrt(np.square(pupil_grid.x) + np.square(pupil_grid.y))/w) *t
	E =  A*np.exp(2*np.pi*1j*(p1 + p2))
	E = np.nan_to_num(E)
	focal_grid = make_focal_grid(pupil_grid, 8, 100, wavelength=l0*1E-6)
	fraunhofer = FraunhoferPropagator(pupil_grid, focal_grid, wavelength_0=l0*1E-6)
	wf = Wavefront(E, wavelength=w)
	wf.total_power = 1
	img = fraunhofer(wf)
	img_out0 = img_out0 +img.intensity
	
imshow_field(np.log10(img_out0 / img_out0.max()), vmin=-6,interpolation='nearest')


plt.colorbar()
plt.show()
plt.close()

img_out60 = 0
phi = 60.0* np.pi/180
devl0 = AM.delR(AM.h0,phi,AM.a,AM.b,l0)

for w in wavelenghts:
	refrac = AM.delR(AM.h0,phi,AM.a,AM.b,w*1E6)
	t = devl0-refrac
	p1 = np.sqrt(np.square(pupil_grid.x) + np.square(pupil_grid.y))/pha
	p2 = (np.sqrt(np.square(pupil_grid.x) + np.square(pupil_grid.y))/w) / t
	E =  A*np.exp(2*np.pi*1j*(p1 + p2))
	E = np.nan_to_num(E)
	focal_grid = make_focal_grid(pupil_grid, 8, 100, wavelength=l0*1E-6)
	fraunhofer = FraunhoferPropagator(pupil_grid, focal_grid, wavelength_0=l0*1E-6)
	wf = Wavefront(E, wavelength=w)
	wf.total_power = 1
	img = fraunhofer(wf)
	img_out60 = img_out60 +img.intensity
	
imshow_field(np.log10(img_out60 / img_out60.max()), vmin=-6,interpolation='nearest')

plt.colorbar()
plt.show()
plt.close()


img_outdiff = img_out60 -img_out0
imshow_field((img_outdiff / img_outdiff.max()),vmin=0.0,vmax=0.1)

plt.colorbar()
plt.show()
plt.close()
