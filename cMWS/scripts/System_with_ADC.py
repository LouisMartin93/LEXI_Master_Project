# importing packages etc... path will need to be changes as to where those folders are.
from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
import sys
path = '/home/martin/Desktop/LEXI ADC project/atmospheric model/'
sys.path.insert(1,path)
path = '/home/martin/Desktop/LEXI ADC project/prism model/'
sys.path.insert(1,path)
import AtmosModel as AM
import PrismFunctions as PF
from astropy.io import fits


#loading the prism data.
from astropy.table import Table
table = Table.read('../../prism model/glass tables/OptimalPrisms.fits')
prism1 = table[0]
prism2 = table[1]
prisms =(prism1,prism2)
print prisms


#atmospheric affects added to mikes PSF!
Dtel=1
GridSize = 900
pupil_grid = make_pupil_grid(GridSize)
aperture = circular_aperture(Dtel)
A = aperture(pupil_grid)
#imshow_field(A)
#plt.show()

pha =  fits.getdata('../cMWS fits/phase.fits')
#plt.imshow(pha,origin='lower',interpolation='nearest')
#plt.show()

pha =  np.hstack(pha)

wlen1 = 0.6E-6
n =11
wlen2 = 0.8E-6 

# wavelenghts in range 0.6 - 0.8 um
wavelenghts = np.linspace(wlen1,wlen2,n)
img_out0 = 0
img_outc = 0
l0 = 0.7
F = 46.2

#zenith and corresponding rotation required
zenith = np.array([0,10,20,30,40,50,60])* np.pi/180
rotAng =np.array([90,95.05,101.1,108.1,116.19,129.11,156.17])* np.pi/180 # corresponding rotation angle for zenith, this is already worked out (not robustly) for the prisms in put.

#indicate the zenith and so rotation working with.
ind = 6 #indicies
z = zenith[ind]
phi = rotAng[ind]

# defining the foacl plane 'mesh', the propagation and the coronagraph
focal_grid = make_focal_grid(pupil_grid,q=10,num_airy=50,focal_length=F, wavelength=l0*1E-6)
fraunhofer = FraunhoferPropagator(pupil_grid, focal_grid, focal_length=F, wavelength_0=l0*1E-6)
coronagraph = PerfectCoronagraph(A)

# calculating the refraction due to the central wavelength, this will be used such that will be realative atmospheric refraction is calculated.
devl0 = AM.delR(AM.h0,z,AM.a,AM.b,l0)

lambda_D = 0.14 /60 * np.pi/10800 #lambda/D in radians

#boolean trigger for ADC
ADCon=False


for w in wavelenghts:
	#on-sky atmospheric dispersion for the 1m aperture
	refrac = AM.delR(AM.h0,z,AM.a,AM.b,w*1E6)
	
	Offset = (devl0-refrac)# assuming central ray is incident on axis
	if ADCon == True:
		
		#relative atmospheric dispersion that needs to be corrected by the angular magnifications
		Offset = Offset*462
		# past through the ADC model (3-D vectorial ray calculation and dispersion)
		vector,rayAng,xyoffset = PF.CompoundDoublet(prism1, prism2, w*1E6, phi=phi, offset=Offset)		
		vector,rayAng,xyoffset = PF.CompoundDoubletMirror(prism2, prism1, w*1E6, phi=-phi, vec_offset=vector, shift=xyoffset)
		# then recorrected assuming perfect optics.
		tx = vector[0,0]/462 
		ty = vector[1,0]/462
		p2 = 2*np.pi*((pupil_grid.x / Dtel  * tx/lambda_D) + (pupil_grid.y / Dtel  * ty/lambda_D))
	else:
		# if ADC is not on then jsut take the refraction.		
		t=Offset
		p2 = 2*np.pi*((pupil_grid.y / Dtel  * t/lambda_D))
	
	#set phase of the cMWS and the E-Field.
	p1 = pha
	E1 =  A*np.exp(1j*(p1+p2))
	E1 = np.nan_to_num(E1)
	
	wf1 = Wavefront(E1, wavelength=w)
	wf1.total_power = 1
	img1 = fraunhofer(wf1)

	wfc1 = coronagraph.forward(wf1)
	imgc1 = fraunhofer(wfc1)

	E2 =  A*np.exp(-1j*(p1-p2))
	E2 = np.nan_to_num(E2)
	wf2 = Wavefront(E2, wavelength=w)
	wf2.total_power = 1
	img2 = fraunhofer(wf2)

	wfc2 = coronagraph.forward(wf2)
	imgc2 = fraunhofer(wfc2)
	
	img_out0 = img_out0 +img1.intensity +img2.intensity
	img_outc = img_outc +imgc1.intensity +imgc2.intensity

reShape = 1000

if ADCon==True:
	titlestr = 'ADC on'
else:
	titlestr = 'ADC off'

fig,axis = plt.subplots(ncols=2,nrows=1)
ax1=axis.flat[0]
im1 = ax1.imshow(np.log10(np.reshape(img_out0,(reShape,reShape)) / img_out0.max()),vmax=0, vmin=-6,interpolation='nearest',origin='lower')
ax1.set_title('PSF at zenith {} for {}'.format(z/np.pi*180,titlestr))

ax2=axis.flat[1]
im2 = ax2.imshow(np.log10(np.reshape(img_outc,(reShape,reShape)) / img_out0.max()),vmax=0, vmin=-6,interpolation='nearest',origin='lower')
ax2.set_title('Coronagraph at zenith {} for {}'.format(z/np.pi*180,titlestr))

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im1, cax=cbar_ax)
plt.show()


