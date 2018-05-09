# tests the beam displacement through the system. 



import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.patches import Rectangle
path = '/home/martin/Desktop/LEXI ADC project/atmospheric model/'
sys.path.insert(1,path)
import AtmosModel as AM

Lambda0 = 0.7 #microns
#loading the prism data.
from astropy.table import Table
import astropy.io.fits as fits
table = Table.read('glass tables/Glass_worksheet2.fits')
glass_by_shape = table.group_by('Prism Type')
prisms_same = glass_by_shape.groups[0]
limit = 1./10.* Lambda0*10**-6  *180/np.pi *3600 *462#arcsec 
print 'limit is {} arcseconds'.format(np.round(limit,4))
rotate = np.linspace(0,360,181)
import PrismFunctions as PF
import time
prism1 = prisms_same[25]
prism2 = prisms_same[78]
prism1['Material'] = 'NSF11'

# double singlets
#theoretical calculations
n1 = eval('PF.'+prism1['Material']+'({})'.format(Lambda0))
n2 = eval('PF.'+prism2['Material']+'({})'.format(Lambda0))
prism2['Wedge angle'] = 5.37#prism1['Wedge angle']*np.round(1+(n1/(1.34*n2)),2)
a1 = prism1['Wedge angle']
a2 = prism2['Wedge angle']

print 'Combination of {} and {} glass with wedge angles of {} and {}'.format(prism1['Material'],prism2['Material'],prism1['Wedge angle'],prism2['Wedge angle'])
sys_rotation = 1.*np.pi/180.
#compound doublet
rotxyvar =[]
d = np.array([0,0,9+6.76685]) #mm distance between prisms and thickness of prisms case 1
d = np.array([0,0,11.3+7.83783999]) #mm distance between prisms and thickness of prisms case 2
rotate = [90,93.84,97.95,102.67,108.59,116.92,131.16] # these are for N-SF11 3 and CaF2 0.5
rotate = [90,95,101,108,117,129,157] # these are for N-SF11 3 and CaF2 5.37
zen = [0,10,20,30,40,50,60]
#prism1['T(mm)'] = 0
#prism2['T(mm)'] = 0
for r in rotate:
	vec,ray_ref,xyoff = PF.CompoundDoublet(prism1,prism2,Lambda0,phi= r*np.pi/180.,offset=0.,Print=False)

	rotx1 = PF.RotMat3Dx(-(a2-a1))
	rotz  = PF.RotMat3Dz(-r*np.pi/180.)
	surf_N = np.matrix([0.,0.,-1.]).T
	surf_N1_ = rotx1*surf_N
	surf_N1 = rotz*surf_N1_
	shift = PF.LineSurfIntersec(surf_N1.getA().T[0,:],d,vec.getA().T[0,:],xyoff,epsilon=1e-16)
	

	vec,ray_ref,xyoff = PF.CompoundDoubletMirror(prism2, prism1, Lambda0, phi=-r*np.pi/180., Print=False, shift=shift, vec_offset=vec)
	rotxyvar.append(xyoff)
	print 'Ref ray: {}'.format(ray_ref[0,0]*60)
	print 'offset on the x-y plane:{}'.format(xyoff)
rotxyvar = np.array(rotxyvar)

prism_radius = prism1['Diameter']/2.
beam_rad = 9./2. *0.05
theta = np.linspace(0,2*np.pi)
circ_x = prism_radius * np.cos(theta)
circ_y = prism_radius * np.sin(theta)
beam_x = beam_rad * np.cos(theta)
beam_y = beam_rad * np.sin(theta) 


plt.scatter(rotxyvar[:,0],rotxyvar[:,1],label='rotation 0-180')
for i in range(len(rotxyvar)):
	plt.plot([0,np.max(rotxyvar[i,0])],[0,np.max(rotxyvar[i,1])],'-b')
plt.plot(beam_x,beam_y,'--k',label='limit of deviaton')
plt.title('Beam steering for Combination of {} and {} glass with wedge angles of {} and {}'.format(prism1['Material'],prism2['Material'],prism1['Wedge angle'],prism2['Wedge angle']),size=17)
plt.xlabel('x-coord (mm)',size=16)
plt.ylabel('y-coord (mm)',size=16)
plt.show()

	
rotxyvar = np.array(rotxyvar)
for i in range(len(rotxyvar)):
	plt.plot(beam_x+np.max(rotxyvar[i,0]),beam_y+np.max(rotxyvar[i,1]),'-',label='offset at z: {}'.format(zen[i]))

plt.plot(beam_x,beam_y,'--k',label='Beam enter')
plt.plot(circ_x,circ_y,'-.k',label='Prism Exit surface')
plt.xlabel('x co-ord (mm)',size=16)
plt.ylabel('y co-ord (mm)',size=16)
plt.title('Lateral shift at exit pupil due to Dispersion',size=18)
plt.legend(loc='best')
plt.grid()
plt.show()

