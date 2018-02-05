'''
This was the first test of the ray tracing through a prism for rotations etc.. proper library of functions can be found in PrismFunctions.py generally abreviated as PF

'''

#contains funcitons
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.patches import Rectangle
path = '/data2/Louis/LEXI ADC project/atmospheric model/'
sys.path.insert(1,path)
import atmos_model as am
n_air = 1.0
R = 14.0 *np.pi/180. * (1./3600.)

def VectorForm(psi,phi):
	#must be in radians
	x = np.sin(psi)*np.cos(phi)
	y = np.sin(psi)*np.sin(phi)
	z = np.cos(psi)
	return np.matrix([x,y,z]).T

def snellLaw(s1,surf_N,n1,n2,Print=False):
	"""
	Inputs are s1 which is a matrix vector describing the incident direction of light

	surf_N is a matrix vector describing the normal of the surface of the prism for the light entering the prism.

	n1 and n2 are the refractive indexes of medium before medium and of the prism option of including n3 if medium the light exits is not the same (i.e. air --> prism --> second prism)
	
	Vectors should be normalised.

	"""
	

	# incident angle
	in_angle = np.arccos(np.dot(-s1.T,surf_N))*180/np.pi
	if Print==True:
		print 'incident angle: {}'.format(in_angle)

	# cross product of surface normal and the incident beam
	a = np.cross(surf_N,s1,axis=0)

	# refractant vector.
	s2 = np.dot((n1/n2),np.cross(-surf_N,a,axis=0))-surf_N*np.sqrt(1 -np.dot((n1/n2)**2 , np.dot(a.T,a)))

	# refraction angle in degrees
	angle = np.arccos(np.dot(s2.T,-surf_N))*180/np.pi
	if Print==True:
		print 'refractive angle: {}'.format(angle)
	return s2
	
	
def RotMat3Dx(beta):
	rotation_matrix_3d_x = np.matrix([[1,0,0],[0,np.cos(beta), -np.sin(beta)],[0,np.sin(beta),np.cos(beta)]])
	return rotation_matrix_3d_x
def RotMat3Dy(beta):
	return np.matrix([[np.cos(beta),0,np.sin(beta)],[0,1,0],[-np.sin(beta),0,np.cos(beta)]])
def RotMat3Dz(beta):
	rotation_matrix_3d_z = np.matrix([[np.cos(beta),-np.sin(beta),0],[np.sin(beta),np.cos(beta),0],[0,0,1]])
	return rotation_matrix_3d_z


def check_bounds(beam_point,alpha,s2,surf_N,C,Tedge=0):
	'''
	returns the height at which the beam hits next surface. 
	'''

	Cprime = beam_point/np.sin(alpha)
	#print (Cprime)
	d = (C-Cprime)*np.cos(alpha) +Tedge
	#print d
	refrac_ang = np.arccos(np.dot(s2.T,-surf_N))
	#print refrac_ang*180/np.pi
	gamma = np.pi/2. - alpha - refrac_ang
	#print gamma*180/np.pi
	c = d*np.tan(gamma)
	#print c
	Bprime = beam_point-c
	print 'new height: {}'.format(Bprime) 
	
	return Bprime 

###--------------------------------------------------###
###-------- Funcs for the refractive indcies --------###
###--------------------------------------------------###

def UVFS(lam):
	B1,B2,B3 = (0.6961663,0.4079426,0.8974794)
	C1,C2,C3 = (0.0684043**2,0.1162414**2,9.896161**2)
	n = np.sqrt(1 + (B1*lam**2)/(lam**2 - C1) +(B2*lam**2)/(lam**2 - C2)+(B3*lam**2)/(lam**2 - C3))
	return np.round(n,4)

def NBK7(lam):
	B1,B2,B3 = (1.03961212,0.231792344,1.01046945)
	C1,C2,C3 = (0.00600069867,0.0200179144,103.560653)
	n = np.sqrt(1 + (B1*lam**2)/(lam**2 - C1) +(B2*lam**2)/(lam**2 - C2)+(B3*lam**2)/(lam**2 - C3))
	return np.round(n,4)

def F2(lam):
	B1,B2,B3 = (1.34533359,0.209073176,0.937357162)
	C1,C2,C3 = (0.00997743871,0.0470450767,111.886764)
	n = np.sqrt(1 + (B1*lam**2)/(lam**2 - C1) +(B2*lam**2)/(lam**2 - C2)+(B3*lam**2)/(lam**2 - C3))
	return np.round(n,4)
def NSF11(lam):
	B1,B2,B3 = (1.73759695,0.313747346,1.89878101)
	C1,C2,C3 = (0.013188707,0.0623068142,155.23629)
	n = np.sqrt(1 + (B1*lam**2)/(lam**2 - C1) +(B2*lam**2)/(lam**2 - C2)+(B3*lam**2)/(lam**2 - C3))
	return np.round(n,4)
def CaF2(lam):
	B1,B2,B3 = (0.5675888,0.4710914,3.8484723)
	C1,C2,C3 = (0.050263605,0.1003909,34.649040)
	n = np.sqrt(1 + (B1*lam**2)/(lam**2 - C1**2) +(B2*lam**2)/(lam**2 - C2**2)+(B3*lam**2)/(lam**2 - C3**2))
	return np.round(n,4)
def BaF2(lam):
	B1,B2,B3 = (0.81070,0.19652,4.52469)
	C1,C2,C3 = (0.10065,29.87,53.82)
	n = np.sqrt(1.33973 + (B1*lam**2)/(lam**2 - C1**2) +(B2*lam**2)/(lam**2 - C2**2)+(B3*lam**2)/(lam**2 - C3**2))
	return np.round(n,4)

def MgF2(lam):
	B1,B2,B3 = (0.48755108,0.39875031,2.310353)
	C1,C2,C3 = (0.04338408,0.09461442,23.793604)
	n = np.sqrt(1 + (B1*lam**2)/(lam**2 - C1**2) +(B2*lam**2)/(lam**2 - C2**2)+(B3*lam**2)/(lam**2 - C3**2))
	return np.round(n,4)







###--------------------------------------------------###
###------------------Test of the model---------------###
###--------------------------------------------------###


from astropy.io import ascii
table = ascii.read('Glass_worksheet.txt')
glass_by_shape = table.group_by('Prism Type')
prisms_same = glass_by_shape.groups[0]
print prisms_same


def Correction(prism1,prism2,lambda0,bandwidth,phi):
	triangle_dims1 = (float(prism1['Diameter']), float(prism1['B']), float(prism1['C']))
	triangle_dims2 = (float(prism2['Diameter']), float(prism2['B']), float(prism2['C']))
	prism_angle1 = np.array([90. - float(prism1['Wedge angle']),float(prism1['Wedge angle']),90.])*np.pi/180.
	prism_angle2 = np.array([90. - float(prism2['Wedge angle']),float(prism2['Wedge angle']),90.])*np.pi/180.
	lambda1 = lambda0 - bandwidth/2.
	lambda2 = lambda0 + bandwidth/2.

	alpha,beta,gamma = prism_angle1
	alpha1,beta1,gamma1 = prism_angle2
	surf_N = VectorForm(0,0)
	surf_N = surf_N/np.linalg.norm(surf_N)
	rotation1 = RotMat3Dy(-beta)
	surf_N2 = rotation1*surf_N
	surf_N2 = surf_N2/np.linalg.norm(surf_N2)
	rotation2 = RotMat3Dy(-beta)*RotMat3Dz(phi*np.pi/180.)*RotMat3Dy(-beta1)
	rotation2 = rotation2/np.linalg.norm(rotation2)
	surf_N3 = rotation2*surf_N
	surf_N3 = surf_N3/np.linalg.norm(surf_N3)

	# raytrace for first lambda of the bandwidth
	n1 = eval(prism1['Material']+'({})'.format(lambda1))
	n2 = eval(prism2['Material']+'({})'.format(lambda1))
	inangle = np.pi/2. - alpha
	s1 = VectorForm(inangle,np.pi)
	beam_point = 12.5
	s2 = snellLaw(-s1,surf_N,n_air,n1)
	s3 = snellLaw(s2,surf_N2,n1,n2)
	ray1 = snellLaw(s3,surf_N3,n2,n_air)
	ray1 = np.arccos(np.dot(ray1.T,-surf_N3))*180/np.pi

	# raytrace for first lambda of the bandwidth
	n1 = eval(prism1['Material']+'({})'.format(lambda2))
	n2 = eval(prism2['Material']+'({})'.format(lambda2))
	inangle = np.pi/2. - alpha
	s2 = snellLaw(-s1,surf_N,n_air,n1)
	s3 = snellLaw(s2,surf_N2,n1,n2)
	ray2 = snellLaw(s3,surf_N3,n2,n_air)
	ray2 = np.arccos(np.dot(ray2.T,-surf_N3))*180/np.pi
	return ray1, ray2


minimisation = {}
lambda0 = 0.7
BW = 0.2
z0=60.
D = 1. #D= 4.2 usig aperture stop to shorten this...
Plot = False


for i in range(len(prisms_same)):
	for j in range(len(prisms_same)):
		prism1 = prisms_same[i]
		prism2 = prisms_same[j]
		correction = []
		for phi in np.linspace(0,180,181):
			ray1,ray2 =Correction(prism1,prism2,0.7,0.2,phi)
			correction.append(ray1-ray2)
		correction = np.array(correction)*3600./462
		diff_R = am.Refraction(lambda0,BW,z0,D,num=180/z0,Plot=False)
		correction =correction[:,0,0]
		diff= correction - diff_R[::-1]*206265.
		limit =   1./10.* lambda0*10**-6 / D *180/np.pi *3600 #arcsec
		mask = abs(diff) > limit
		uncorrected=diff[mask]

		minimisation[sum(abs(uncorrected))]= (i,j)
		if Plot:
			#then plotting the two optimal combiniation plots (i.e. with the correction)
			rect = Rectangle((0, -limit), 180, 2*limit, color='g', alpha=0.3,label='limit bounadry')
			fig = plt.figure()
			ax1 = fig.add_subplot(121)
			ax2 = ax1.twiny()
			ax1.add_patch(rect)
			phi2 = np.linspace(0,z0,(180/z0)*z0+1)
			phi1 = np.linspace(0,180,181)
			lns1 = ax1.plot(phi1,correction,'r',label='dispersion correction')
			lns2 = ax2.plot(phi2,diff_R*206265.,label = 'atmospheric dispersion')
			ax1.set_xlabel('Angle of rotation of prisms with repect to one another (degree)')
			ax1.set_ylabel('Compensation to Dispersion (")')
			ax2.set_xlabel('Angle from Zenith')
			ax2.set_xlim(np.max(phi2),0)
			lns =lns1+lns2
			labs = [l.get_label() for l in lns]
			ax1.legend(lns, labs, loc='best')
			rect = Rectangle((0, -limit), 180, 2*limit, color='g', alpha=0.3,label='limit bounadry')
			ax = fig.add_subplot(122)
			ax.plot(phi2,abs(diff[::-1]),'-k',label='correction to dispersion')
			ax.add_patch(rect)
			ax.set_xlabel('Angle of rotation of prisms with repect to one another (degree)')
			ax.set_ylabel('Correction to Dispersion (")')
			plt.suptitle('Prism comb. of {0} and {1} with beam deviation {2} and {3}'.format(prism1['Material'],prism2['Material'], prism1['deviation'], prism2['deviation']),size=16)
			plt.legend(loc='best')
			plt.show()
			plt.close()
		


min_key = np.min(minimisation.keys())
print min_key

combo = minimisation[min_key]
print combo

prism1 = prisms_same[combo[0]]
prism2 = prisms_same[combo[1]]






#plotting the ideal case
#first doing all the calculations
correction = []
for phi in np.linspace(0,180,181):
	ray1,ray2 =Correction(prism1,prism2,0.7,0.2,phi)
	correction.append(ray1-ray2)
correction = np.array(correction)*3600./462
diff_R = am.Refraction(lambda0,BW,z0,D,num=180/z0,Plot=False)
print correction.shape
correction =correction[:,0,0]
diff= correction - diff_R[::-1]*206265. 


#then plotting the two optimal combiniation plots (i.e. with the correction)
limit =   1./10.* lambda0*10**-6 / D *180/np.pi *3600 #arcsec
limitBlue = 1./10.* 0.6*10**-6 / D *180/np.pi *3600 #arcsec
limitRed = 1./10.* 0.8*10**-6 / D *180/np.pi *3600 #arcsec
rect = Rectangle((0, -limit), 180, 2*limit, color='g', alpha=0.3,label='limit bounadry')
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = ax1.twiny()
ax1.add_patch(rect)
phi2 = np.linspace(0,z0,(180/z0)*z0+1)
phi1 = np.linspace(0,180,181)
lns1 = ax1.plot(phi1,correction,'r',label='dispersion correction')
lns2 = ax2.plot(phi2,diff_R*206265.,label = 'atmospheric dispersion')
ax1.set_xlabel('Angle of rotation of prisms with repect to one another (degree)')
ax1.set_ylabel('Compensation to Dispersion (")')
ax2.set_xlabel('Angle from Zenith')
ax2.set_xlim(np.max(phi2),0)
lns =lns1+lns2
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc='best')
rectG = Rectangle((0, -limit), 180, 2*limit, color='g', alpha=0.3,label='limit bounadry')
rectB = Rectangle((0, -limitBlue), 180, 2*limitBlue, color='b', alpha=0.2,label='limit bounadry')
rectR = Rectangle((0, -limitRed), 180, 2*limitRed, color='r', alpha=0.2,label='limit bounadry')
ax = fig.add_subplot(122)
ax.plot(phi2,abs(diff[::-1]),'-k',label='correction to dispersion')
#ax.add_patch(rectR)
ax.add_patch(rectG)
#ax.add_patch(rectB)
ax.set_xlabel('Angle of rotation of prisms with repect to one another (degree)')
ax.set_ylabel('Correction to Dispersion (")')
plt.suptitle('Prism comb. of {0} and {1} with beam deviation {2} and {3}'.format(prism1['Material'],prism2['Material'], prism1['deviation'], prism2['deviation']),size=16)
plt.legend(loc='best')
plt.show()
plt.close()
