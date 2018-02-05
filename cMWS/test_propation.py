from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
import sys
path ='/data2/Louis/LEXI ADC project/atmospheric model/'
sys.path.insert(1,path)
import AtmosModel as AM
Dtel=1
GridSize =512
pupil_grid = make_pupil_grid(GridSize)
aperture = circular_aperture(Dtel)

A = aperture(pupil_grid)
#imshow_field(A)
#plt.show()
n = 11
l0 = 0.7
bw = 0.2

l1 = l0 - bw/2.
l2 = l0 + bw/2.
phi = 60 *np.pi/180.

#onsky refraction caused by atmosphere at zenith angle phi 
devl0 = AM.delR(AM.h0,phi,AM.a,AM.b,l0)
devl1 = AM.delR(AM.h0,phi,AM.a,AM.b,l1)
devl2 = AM.delR(AM.h0,phi,AM.a,AM.b,l2)


# wavelenghts at in range 0.6 -0.8 um
wavelenghts = np.linspace(l1,l2,n)*1E-6
# the corresponding refraction relitive to wavelenght.
theta = np.linspace(devl0-devl1,devl0-devl2,n)

img_out = 0
img_ref = 0
for w,t in zip(wavelenghts,theta):
	refrac = AM.delR(AM.h0,phi,AM.a,AM.b,w*1E6)
	t = devl0-refrac
	phase =  (pupil_grid.x / w )*t #+ (pupil_grid.x / Dtel) 
	E = A * np.exp(2*np.pi *1j * phase)
	focal_grid = make_focal_grid(pupil_grid, 8, 15, wavelength=0.6E-6)
	fraunhofer = FraunhoferPropagator(pupil_grid, focal_grid, wavelength_0=0.6E-6)
	wf = Wavefront(E, wavelength=w)
	wf.total_power = 1


	img = fraunhofer(wf)
	img_out = img_out + img.intensity
	wf = Wavefront(A, wavelength=w)
	wf.total_power = 1


	img = fraunhofer(wf)
	img_ref = img_ref + img.intensity
	
imshow_field(np.log10(img_out / img_out.max()), vmin=-5)

plt.title('On-sky broadband PSF due to atmospheric refraction at {} zenith'.format(phi*180/np.pi))
plt.colorbar()
plt.show()
plt.close()

'''
img_out=0
for w,t in zip(wavelenghts,theta):
	refrac = AM.delR(AM.h0,phi,AM.a,AM.b,w*1E6)
	t1 = 462*(devl0-refrac)
	phase =  (pupil_grid.x / w )*t1 #+ (pupil_grid.x / Dtel) 
	E = A * np.exp(2*np.pi *1j * phase)
	focal_grid = make_focal_grid(pupil_grid, 8, 15, wavelength=0.7E-6)
	fraunhofer = FraunhoferPropagator(pupil_grid, focal_grid, wavelength_0=0.7E-6)
	wf = Wavefront(E, wavelength=w)
	wf.total_power = 1


	img = fraunhofer(wf)
	img_out = img_out + img.intensity
	wf = Wavefront(A, wavelength=w)
	wf.total_power = 1


	img = fraunhofer(wf)
	img_ref = img_ref + img.intensity
	
imshow_field(np.log10(img_out / img_out.max()), vmin=-5)

plt.title('Broadband PSF due to atmospheric refraction at {} zenith'.format(phi*180/np.pi))
plt.colorbar()
plt.show()
plt.close()
'''

