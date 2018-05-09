# this will be used as a test calculation for the doublet model using the prisms as well as functions from the PrismFuncitons

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

# general loading of data.
Lambda = 0.7 #microns
re = 6371.0*10**3 # mean radius of the earth (m)
h0 = 2350 #altitude of observatory on la palma (m)
Mma = 0.0289644 # mean molar mass of dry air (kg/mol)
Lapse_rate = 6.49/1000. # lapse rate of  atmosphere (K/m)
g = 9.81 #acceleration due to gravity (m/s^2)
T0 = 4+273.15 # mean temperature at 00:00 in November at WHT (K)
R = 8.31447 # universal gas constant (J/(mol*K))
D =  1#4.26 #telescope diameter.

q = (g*Mma/(R*Lapse_rate))
H0 = T0/Lapse_rate
a = (q-1.)/H0
b = (q-1.)*(q-2.)/(2.*H0**2)
r0 = re + h0

P = 0.765852


#loading the prism data.
from astropy.table import Table
import astropy.io.fits as fits
table = Table.read('../glass tables/Glass_worksheet2.fits')
glass_by_shape = table.group_by('Prism Type')
prisms_same = glass_by_shape.groups[0]
limit = 1./10.* Lambda*10**-6 / D *180/np.pi *3600 *462#arcsec 
print 'limit is {} arcseconds'.format(np.round(limit,4))

l0 = 0.75
bw = 0.3
l1 = l0 - bw/2.
l2 = l0 + bw/2.
print l1,l2


prism1a = prisms_same[25] #25
prism2a = prisms_same[78] #78
prism1b = prisms_same[26] #25
prism2b = prisms_same[79] #78
prism1a['Material'] = 'NSF11'
prism1b['Material'] = 'NSF11'
n1 = eval('PF.'+prism1a['Material']+'({})'.format(l0))
n2 = eval('PF.'+prism2a['Material']+'({})'.format(l0))

prism1a['Wedge angle'] = 3. #- 5./60.
prism2a['Wedge angle'] = 5.35 #+ 5./60.
prism1b['Wedge angle'] = 3. #- 5./60.
prism2b['Wedge angle'] = 5.35 #+ 5./60.
print 'Combination of {} and {} glass with wedge angles of {} and {}'.format(prism1a['Material'],prism2a['Material'],prism1a['Wedge angle'],prism2a['Wedge angle'])
#going to compare the two values
zenith = np.array([0,10,20,30,40,50,60])

n = 11
mid = (n-1)/2
wavelenghts = np.linspace(l1,l2,n)

deg=180
rotate = np.linspace(80,deg,101) 
#optRot=[90,93.84,97.95,102.67,108.59,116.92,131.16]
#optRot=[90,95.12,100.62,106.99,115.13,127.11,151.28]

OptRot=[]
for z in zenith:	
	offsets = (AM.Refraction(wavelenghts,z*np.pi/180.,h0,T0,P))*462
	offsets =  (offsets[mid] - offsets)
	disp4phi=[]
	
	for phi in rotate:
		#optphi = optphi * np.pi/180.
		phi =  phi*np.pi/180. 
		Vectors=[]
		RayAngles=[]
		XYoffsets=[]
		for Lambda,Offset in zip(wavelenghts,offsets):
			vector,rayAng,xyoffset = PF.CompoundDoublet(prism1a, prism2a, Lambda, phi=phi, offset=Offset)
			
			vector,rayAng,xyoffset = PF.CompoundDoubletMirror(prism2b, prism1b, Lambda, phi=-phi, vec_offset=vector, shift=xyoffset)
			
			Vectors.append(vector)
			RayAngles.append(rayAng[0,0])
			XYoffsets.append(xyoffset)
			
		maxindex = np.argmax(RayAngles)
		minindex = np.argmin(RayAngles)
		#dispersion = RayAngles[maxindex] - RayAngles[minindex]
		dispersion_ = np.arccos(np.dot(Vectors[maxindex].T,Vectors[minindex]))*180/np.pi
		#disp4phi.append(dispersion*3600) #dispersion in arcseconds
		disp4phi.append(dispersion_[0,0]*3600) #dispersion in arcseconds
		
	disp4phi = np.array(disp4phi)
	plt.plot(rotate,disp4phi,label=str(z))
	
	OptRot.append(rotate[np.argmin(disp4phi)])
	
#print np.array(zip(zenith,OptRot))

plt.plot([80,deg],[limit,limit],'--k',label='lim')
plt.legend(loc ='best')
plt.xlabel('Rotation of prisms with respect to one another (deg)')
plt.ylabel('Dispersion (arcseconds)')
plt.title('Combination of {} and {} glass with wedge angles of {} and {}'.format(prism1a['Material'],prism2a['Material'],prism1a['Wedge angle'],prism2a['Wedge angle']))
plt.axis('tight')
plt.show()

print np.array(zip(zenith,OptRot))
print 'Beam Deviation','min Dispersion','rotation'
for z,optphi in zip(zenith,OptRot):	
	offsets = (AM.Refraction(wavelenghts,z*np.pi/180.,h0,T0,P))*462
	offsets =  (offsets[mid] - offsets)
	disp4phi=[]
	
	
	#optphi = optphi * np.pi/180.
	phi =  optphi*np.pi/180. 
	Vectors=[]
	RayAngles=[]
	XYoffsets=[]
	for Lambda,Offset in zip(wavelenghts,offsets):
		vector,rayAng,xyoffset = PF.CompoundDoublet(prism1a, prism2a, Lambda, phi=phi, offset=Offset)
			
		vector,rayAng,xyoffset = PF.CompoundDoubletMirror(prism2b, prism1b, Lambda, phi=-phi, vec_offset=vector, shift=xyoffset)
			
		Vectors.append(vector)
		RayAngles.append(rayAng[0,0])
		XYoffsets.append(xyoffset)
			
	maxindex = np.argmax(RayAngles)
	minindex = np.argmin(RayAngles)
	#dispersion = RayAngles[maxindex] - RayAngles[minindex]
	dispersion_ = np.arccos(np.dot(Vectors[maxindex].T,Vectors[minindex]))*180/np.pi
	#disp4phi.append(dispersion*3600) #dispersion in arcseconds
	disp4phi.append(dispersion_[0,0]*3600) #dispersion in arcseconds
	print RayAngles[mid] *60, disp4phi,optphi-90



z = np.linspace(0,np.pi/3.,100)

a_vals = [0.56,0.5,0.4,0.3,0.2,0.1]

for a in a_vals:
	func =  np.arccos(-a*np.tan(z))*180./np.pi
	plt.plot(z*180./np.pi, func-90.,label='a={}'.format(a))


plt.plot(zenith,(np.array(OptRot)-90),label='My values')
plt.xlabel('zenith (deg)')
plt.ylabel('rotation (deg)')
plt.legend(loc='best')
plt.show()







