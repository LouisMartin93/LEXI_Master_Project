#-----prism functions folder-----#

import numpy as np
import sys
path = '~/Desktop/LEXI ADC project/atmospheric model/'
sys.path.insert(1,path)
import AtmosModel as AM

n_air =1.
###------------------------------------------------###
###-------------Vector functions-------------------###
###------------------------------------------------###


def VectorForm(psi,phi):
	#must be in radians
	x = np.sin(psi)*np.cos(phi)
	y = np.sin(psi)*np.sin(phi)
	z = np.cos(psi)
	return np.matrix([x,y,z]).T
	
	
def RotMat3Dx(beta):
	rotation_matrix_3d_x = np.matrix([[1,0,0],[0,np.cos(beta), -np.sin(beta)],[0,np.sin(beta),np.cos(beta)]])
	return rotation_matrix_3d_x

def RotMat3Dy(beta):
	return np.matrix([[np.cos(beta),0,np.sin(beta)],[0,1,0],[-np.sin(beta),0,np.cos(beta)]])

def RotMat3Dz(beta):
	rotation_matrix_3d_z = np.matrix([[np.cos(beta),-np.sin(beta),0],[np.sin(beta),np.cos(beta),0],[0,0,1]])
	return rotation_matrix_3d_z



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





#-----------------------------------------------#
#-----------------Prism funcitons---------------#
#-----------------------------------------------#

def snellLaw(s1,surf_N,n1,n2,Print=False):
	"""
	Inputs are s1 which is a matrix vector describing the incident direction of light

	surf_N is a matrix vector describing the normal of the surface of the prism for the light entering the prism.

	n1 and n2 are the refractive indexes of medium before medium and of the prism option of including n3 if medium the light exits is not the same (i.e. air --> prism --> second prism)
	
	Vectors should be normalised.

	"""
	if Print == True:
		print surf_N
		print np.linalg.norm(surf_N)
		print s1
		print np.linalg.norm(s1)

	surf_N = surf_N/np.linalg.norm(surf_N)
	s1 = s1/np.linalg.norm(s1)
	if Print == True:
		print surf_N
		print s1
	# incident angle
	in_angle = np.arccos(np.dot(-s1.T,surf_N))*180/np.pi
	if Print==True:
		print 'incident angle: {}'.format(in_angle)

	# cross product of surface normal and the incident beam
	a = np.cross(surf_N,s1,axis=0)
	if Print==True:
		print a 
		print a/np.linalg.norm(a)
	# refractant vector.
	s2 = np.dot((n1/n2),np.cross(-surf_N,a,axis=0))-surf_N*np.sqrt(1 -np.dot((n1/n2)**2 , np.dot(a.T,a)))

	# refraction angle in degrees
	angle = np.arccos(np.dot(s2.T,-surf_N))*180/np.pi
	if Print==True:
		print 'refractive angle: {}'.format(angle)
	return s2


def check_bounds(beam_point,alpha,s2,surf_N,C,Tedge=0):
	'''
	Returns the height at which the beam hits next surface compared to the first surface.
	Inputs are:
	beam_point, the hight at which the centre of the beam hits the suface.
	alpha the wedge angle, 
	the vector of the ray after refraction,
	the surface vector of the first surface,
	the hypotenuse of the wedge,
	and the option to add a thickness.
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

def CorrectionSingle(prism1,prism2,Lambda,phi,offset=0.):
	'''
	Inputs Prism1 and Prism 2 are either tables or dictionaries which must include the following keys 'Material','Diameter','B','C',Wedge angle'
	Lambda is the wavelenght of ray to pass through a Rotational doublet prism.
	Phi is the rotation of the prisms with respect to one another given in degrees (converted to radians)
	Offset is the value at which the angle at which the ray enters the system if 0 the ray enters parrallel to the x axis and at an angle of beta relitive to the normal of the first surface.

	The output is angle of the ray relative to the norm of a surface perpendicular to optical axis given in degrees. 
	'''	
	triangle_dims1 = (float(prism1['Diameter']), float(prism1['B']), float(prism1['C']))
	triangle_dims2 = (float(prism2['Diameter']), float(prism2['B']), float(prism2['C']))
	prism_angle1 = np.array([90. - float(prism1['Wedge angle']),float(prism1['Wedge angle']),90.])*np.pi/180.
	prism_angle2 = np.array([90. - float(prism2['Wedge angle']),float(prism2['Wedge angle']),90.])*np.pi/180.

	alpha,beta,gamma = prism_angle1
	alpha1,beta1,gamma1 = prism_angle2
	
	#all rotations
	rotx1 = RotMat3Dx(alpha)
	rotz = RotMat3Dz(phi*np.pi/180.)
	rotx2 = RotMat3Dx(beta)
	rotx3 = RotMat3Dx(beta1)

	#all surface calculations using co-ords sys of +y is north +z is east
	# rays travel in +z direction
	#first surface
	surf_N = np.matrix([0.,1.,0.]).T
	surf_N = rotx1*surf_N

	#second surface
	surf_N2 = rotx2*surf_N
	
	
	#third surface
	surf_N3 = rotz*rotx3*surf_N2
	
	# raytrace for first lambda of the bandwidth
	n1 = eval(prism1['Material']+'({})'.format(Lambda))
	n2 = eval(prism2['Material']+'({})'.format(Lambda))

	if offset > 0:
		phi1 = np.pi/2.
		psi = abs(offset)
	else:
		phi1 = np.pi*1.5
		psi = abs(offset)
	s1 = VectorForm(psi,phi1)
	beam_point = 12.5
	s2 = snellLaw(-s1,surf_N,n_air,n1)
	s3 = snellLaw(s2,surf_N2,n1,n2)
	ray1 = snellLaw(s3,surf_N3,n2,n_air)
	ray1 = np.arccos(np.dot(ray1.T,-surf_N2))*180/np.pi # deviation of ray wrt flat surface! instead of surface of final prism.
	return ray1

def CorrectionBroad(prism1,prism2,lambda0,bandwidth,phi,z):
	'''
	Runs throught a 2 prism system for a given bandwidth about a central wavelenght. It ruturns the ray angle or vector for the highest and lowest wavelengths for a relative rotation of prisms with respect to one another. the angles are w.r.t the final surface that the light passes from. The incident angles for each wavelenght are assumed to be the same, however will be different, realistically, depending on zenith and wavelenght. 
	inputs are as follows:
	prism1 and prims2 are tables containing at least the folloing data, diameter of wedge (Diameter), wedge angle (Wedge angle), material (Material), thickness (B), wedge hypothenuse (C)
	'''
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
	inangle = np.pi/2. - alpha #- 462*AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,Lambda)
	s1 = VectorForm(inangle,np.pi)
	beam_point = 12.5
	s2 = snellLaw(-s1,surf_N,n_air,n1)
	s3 = snellLaw(s2,surf_N2,n1,n2)
	ray1 = snellLaw(s3,surf_N3,n2,n_air)
	ray1 = np.arccos(np.dot(ray1.T,-surf_N3))*180/np.pi

	# raytrace for second lambda of the bandwidth
	n1 = eval(prism1['Material']+'({})'.format(lambda2))
	n2 = eval(prism2['Material']+'({})'.format(lambda2))
	inangle = np.pi/2. - alpha
	s2 = snellLaw(-s1,surf_N,n_air,n1)
	s3 = snellLaw(s2,surf_N2,n1,n2)
	ray2 = snellLaw(s3,surf_N3,n2,n_air)
	ray2 = np.arccos(np.dot(ray2.T,-surf_N3))*180/np.pi
	return ray1, ray2


def CompoundDoublet(prism1,prism2,Lambda,offset=0.):
	'''
	This function will calculate the deviation of a compound prism made of two materials, i.e. a doublet.

	The function will calculate the ray tracing through vectors as this will be the easiest method for including 
	any rotations necessary.
	
	Inputs will include a lambda, the prism dictionaries, and offset (from the optical axis) 
	Output will be a ray vector (s4) and deviation in degrees relitive to a surface perpendicular to the optical axis (ray1).
	'''
	
	triangle_dims1 = (float(prism1['Diameter']), float(prism1['B']), float(prism1['C']))
	triangle_dims2 = (float(prism2['Diameter']), float(prism2['B']), float(prism2['C']))
	
	# need the angles in radians #
	alpha1 = float(prism1['Wedge angle'])*np.pi/180.
	alpha2 = float(prism2['Wedge angle'])*np.pi/180.
	
	beta1 = 0.0 # this implies that the surface is perpendicular to the optical axis
	gamma1 = alpha1
	gamma2 = gamma1
	beta2 = (gamma2 - alpha2)
	
	#print(alpha1,beta1,gamma1)
	#print(alpha2,beta2,gamma2)
	#all rotations
	rotx1 = RotMat3Dx(beta1)
	rotx2 = RotMat3Dx(-alpha1)
	rotx3 = RotMat3Dx(alpha2)
	#print rotx2,rotx3

	#all surface calculations using co-ords sys of +y is north +z is east
	#rays travel in +z direction
	#first surface
	surf_N = np.matrix([0.,0.,-1.]).T
	surf_N = rotx1*surf_N
	#print surf_N

	#second surface
	surf_N2 = rotx2*surf_N
	#print surf_N2
	
	#Third Surface - exit surface of compound prism
	surf_N3 =  rotx3*surf_N2
	#print surf_N3

	#going to test to make sure deviation is correct here
	# raytrace for first lambda of the bandwidth
	n1 = eval(prism1['Material']+'({})'.format(Lambda))
	n2 = eval(prism2['Material']+'({})'.format(Lambda))

	if offset > 0:
		phi1 = np.pi/2.
		psi = abs(offset)
	else:
		phi1 = np.pi*1.5
		psi = abs(offset)
	s1 = VectorForm(psi,phi1)
	beam_point = 12.5
	s2 = snellLaw(-s1,surf_N,n_air,n1)
	s3 = snellLaw(s2,surf_N2,n1,n2)
	s4 = snellLaw(s3,surf_N3,n2,n_air)
	ray1 = np.arccos(np.dot(s4.T,-surf_N))*180/np.pi # deviation of ray wrt flat surface! instead of surface of final prism.
	return s4, ray1 

def DoubleDoublet(prism1,prism2,Lambda,offset=0.,phi=180.):
	triangle_dims1 = (float(prism1['Diameter']), float(prism1['B']), float(prism1['C']))
	triangle_dims2 = (float(prism2['Diameter']), float(prism2['B']), float(prism2['C']))
	# need the angles in radians #
	alpha1 = float(prism1['Wedge angle'])*np.pi/180.
	alpha2 = float(prism2['Wedge angle'])*np.pi/180.
	
	beta1 = 0.0 # this implies that the surface is perpendicular to the optical axis
	gamma1 = alpha1
	gamma2 = gamma1
	beta2 = (gamma2 - alpha2)
	
	
	alpha3 = alpha2 - alpha1
	#print(alpha3) # radians
	s4,ray1 = CompoundDoublet(prism1,prism2,Lambda,offset=offset)
	
	# surfaces of the next set of compound prisms
	#all rotations
	rotx4 = RotMat3Dx(-alpha3)
	rotx5 = RotMat3Dx(alpha2)
	rotx6 = RotMat3Dx(-alpha1)
	rotz = RotMat3Dz(phi*np.pi/180.)
	#print rotz

	#all surface calculations using co-ords sys of +y is north +z is east
	#rays travel in +z direction
	#first surface
	surf_N = np.matrix([0.,0.,-1.]).T
	surf_N4 = rotz*rotx4*surf_N
	#print surf_N

	#second surface
	surf_N5 = rotx5*surf_N4
	#print surf_N2
	
	#Third Surface - exit surface of compound prism
	surf_N6 =  rotx6*surf_N5
	#print surf_N3
	
	#order has to be reversed for the refractive indeices
	n2 = eval(prism1['Material']+'({})'.format(Lambda))
	n1 = eval(prism2['Material']+'({})'.format(Lambda))
	
	s5 = snellLaw(s4,surf_N4,n_air,n1)
	s6 = snellLaw(s5,surf_N5,n1,n2)
	s7 = snellLaw(s6,surf_N6,n2,n_air)
	ray2 = np.arccos(np.dot(-s7.T,surf_N))*180/np.pi # deviation of ray wrt flat surface! instead of surface of final prism.
	return s7, ray2


def DoubleDoublet2(prism1,prism2,Lambda,offset=0.,phi=180.):
	# need the angles in radians #
	alpha1 = float(prism1['Wedge angle'])*np.pi/180.
	alpha2 = float(prism2['Wedge angle'])*np.pi/180.
	rot1 = alpha2-alpha1
	# surfaces of the next set of compound prisms
	#all rotations
	rotx1 = RotMat3Dx(rot1)
	rotx2 = RotMat3Dx(alpha1)
	rotx3 = RotMat3Dx(-alpha2)
	rotz = RotMat3Dz(phi*np.pi/180.)
	rotx4 = RotMat3Dx(-alpha2)
	rotx5 = RotMat3Dx(alpha1)
	

	#all surface calculations using co-ords sys of +y is north +z is east
	#rays travel in +z direction

	#first surface
	surf_N = np.matrix([0.,0.,-1.]).T
	print surf_N
	surf_N1 = rotx1*surf_N
	print surf_N1
	#second surface
	surf_N2 = rotx2*surf_N1
	print surf_N2
	#Third Surface - exit surface of first compound prism
	surf_N3 =  rotx3*surf_N2
	print surf_N3 #should be flat
	#forth surface
	surf_N4 = rotz*rotx4*surf_N3
	print surf_N4

	#fifth surface
	surf_N5 = rotx5*surf_N4
	print surf_N5
	
	#going to test to make sure deviation is correct here
	# raytrace for first lambda of the bandwidth
	n1 = eval(prism1['Material']+'({})'.format(Lambda))
	n2 = eval(prism2['Material']+'({})'.format(Lambda))

	if offset > 0:
		phi1 = np.pi/2.
		psi = abs(offset)
	else:
		phi1 = np.pi*1.5
		psi = abs(offset)
	s1 = VectorForm(psi,phi1)
	print s1
	s2 = snellLaw(-s1,surf_N1,n_air,n1)
	print s2
	s3 = snellLaw(s2,surf_N2,n1,n2)
	print s3
	s4 = snellLaw(s3,surf_N3,n2,n2)
	print s4
	s5 = snellLaw(s4,surf_N4,n2,n1)
	print s5
	s6 = snellLaw(s5,surf_N5,n1,n_air)
	ray1 = np.arccos(np.dot(s6.T,-surf_N))*180/np.pi # deviation of ray wrt flat surface! instead of surface of final prism.
	print s6, ray1 

	
