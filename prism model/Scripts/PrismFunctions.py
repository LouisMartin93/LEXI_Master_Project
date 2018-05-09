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

def NLLF6(lam):
	B1,B2,B3 = (1.22932376,7.55174062*10**-2,9.88839837*10**-1)
	C1,C2,C3 = (8.77776451*10**-3,4.80429329*10**-2,1.12136098*10**2)
	n = np.sqrt(1 + (B1*lam**2)/(lam**2 - C1) +(B2*lam**2)/(lam**2 - C2)+(B3*lam**2)/(lam**2 - C3))
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
	#surf_N = surf_N/np.linalg.norm(surf_N)
	#s1 = s1/np.linalg.norm(s1)
	if Print == True:
		print surf_N
		print s1
	# incident angle
	in_angle = np.arccos(np.dot(-s1.T,surf_N))*180/np.pi
	if Print==True:
		print 'incident angle: {}'.format(in_angle)

	# cross product of surface normal and the incident beam
	a = np.cross(surf_N,s1,axis=0)
	#if Print==True:
		#print a 
		#print a/np.linalg.norm(a)
	# refractant vector.
	s2 = np.dot((n1/n2),np.cross(-surf_N,a,axis=0))-surf_N*np.sqrt(1 -np.dot((n1/n2)**2 , np.dot(a.T,a)))

	# refraction angle in degrees
	angle = np.arccos(np.dot(s2.T,-surf_N))*180/np.pi
	if Print==True:
		print 'refractive angle: {}'.format(angle)
	return s2 ,angle

def LineSurfIntersec(planeNormal,planePoint,rayDirection,rayPoint,epsilon=1e-16):
	ndotu = planeNormal.dot(rayDirection) 
	if abs(ndotu) < epsilon:
	    print ("no intersection or line is within plane")
	w = rayPoint - planePoint
	si = -planeNormal.dot(w) / ndotu
	Psi = w + si * rayDirection + planePoint
	return Psi

def check_bounds(h0,alpha,s2,surf_N,C,Tedge=0):
	'''
	Returns the height at which the beam hits next surface compared to the first surface.
	Inputs are:
	beam_point, the hight at which the centre of the beam hits the suface.
	alpha the wedge angle, 
	the vector of the ray after refraction,
	the surface vector of the first surface,
	the diameter of the wedge (C)
	and the option to add a thickness.
	'''
	
	Cprime = C-h0
	refrac_ang = np.arccos(np.dot(s2.T,-surf_N))
	#print refrac_ang*180/np.pi
	gamma = np.pi/2. + refrac_ang - alpha
	#print gamma*180/np.pi
	l = np.sin(alpha) * Cprime / np.sin(gamma)
	deltay = l*np.sin(refrac_ang)
	d = l*np.cos(refrac_ang)
	
	S = h0 + deltay
	print 'new height: {}'.format(S) 
	
	return S,d,deltay # the z-y positions and the lateral shift in beam.

def CorrectionSingle(prism1,prism2,Lambda,phi,offset=0.,Print=False,shift=0.):
	'''
	Inputs Prism1 and Prism 2 are either tables or dictionaries which must include the following keys 'Material','Diameter','B','C',Wedge angle'
	Lambda is the wavelenght of ray to pass through a Rotational doublet prism.
	Phi is the rotation of the prisms with respect to one another given in radians
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
	rotz = RotMat3Dz(phi)
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
	C = triangle_dims1[0]
	A1 = triangle_dims1[1]
	A2 = triangle_dims2[1]
	d = A1 - C/2. * np.tan(beta)
	x = A2 - C/2. * np.tan(beta1) #C/2. * np.tan(beta1)/np.cos(beta1)# for the case of a compound prism.
	ray_start = np.array([0,0,0]) +shift
	point_on_surf2 = np.array([0,0,d])
	point_on_surf3 = np.array([0,0,x+d])

	s2,_ = snellLaw(-s1,surf_N,n_air,n1,Print=Print)
	#need to make copies as the matrcies do not handle well with the line segments also need to be transposed to fit function
	s2_copy = s2.getA().T[0,:]
	surf_N2_copy = surf_N2.getA().T[0,:]
	x_y_z1 = LineSurfIntersec(surf_N2_copy,point_on_surf2,s2_copy,ray_start)
	

	s3,_ = snellLaw(s2,surf_N2,n1,n2,Print=Print)
	#similarly need to make copies so that values arent over written.
	s3_copy = s3.getA().T[0,:]
	surf_N3_copy = surf_N3.getA().T[0,:]
	
	x_y_z2 = LineSurfIntersec(surf_N3_copy,point_on_surf3,s3_copy,x_y_z1)
	
	s4,_ = snellLaw(s3,surf_N3,n2,n_air,Print=Print)
	
	ray1 = np.arccos(np.dot(s4.T,-surf_N2))*180/np.pi # deviation of ray wrt flat surface! instead of surface of final prism.
	x_y_offset = np.array([x_y_z2[0],x_y_z2[1]])
	return s4,ray1,x_y_offset



def CompoundDoublet(prism1,prism2,Lambda,phi=0,offset=0.,Print=False,shift=0):
	'''
	This function will calculate the deviation of a compound prism made of two materials, i.e. a doublet.

	The function will calculate the ray tracing through vectors as this will be the easiest method for including 
	any rotations necessary.
	
	Inputs will include a lambda, the prism dictionaries,phi and offset (from the optical axis). phi and offset must be in radians
	Output will be a ray vector (s4) and deviation in degrees relitive to a surface perpendicular to the optical axis (ray1).
	'''
	#prism dimensions
	triangle_dims1 = (float(prism1['Diameter']), float(prism1['B']), float(prism1['C']))
	triangle_dims2 = (float(prism2['Diameter']), float(prism2['B']), float(prism2['C']))
	T_edge1 = float(prism1['T(mm)'])
	T_edge2 = float(prism2['T(mm)'])
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
	rotz = RotMat3Dz(phi)
	#print rotx2,rotx3

	#all surface calculations using co-ords sys of +y is north +z is east
	#rays travel in +z direction
	C = triangle_dims1[0]
	A1 = triangle_dims1[1]
	A2 = triangle_dims2[1]
	d = C/2. * np.tan(alpha1) + T_edge1
	x = C/2. * np.tan(alpha2)/np.cos(alpha2) +T_edge2/np.cos(beta2)# for the case of a compound prism.
	ray_start = np.array([0,0,0]) + shift
	point_on_surf2 = np.array([0,0,d])
	point_on_surf3 = np.array([0,0,x+d])
	
	#first surface
	surf_N = np.matrix([0.,0.,-1.]).T
	surf_N1 = rotx1*surf_N
	#print surf_N

	#second surface
	surf_N2_ = rotx2*surf_N1
	surf_N2 = rotz*surf_N2_
	#print surf_N2
	
	#Third Surface - exit surface of compound prism
	surf_N3_ =  rotx3*surf_N2_
	surf_N3 = rotz*surf_N3_
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
	# need to make copies of the surfaces so values are not overwritten as well as transpose so that they fit the function
	s2,_ = snellLaw(s1,surf_N1,n_air,n1,Print=Print)
	s2_copy = s2.getA().T[0,:]
	surf_N2_copy = surf_N2.getA().T[0,:]
	x_y_z1 = LineSurfIntersec(surf_N2_copy,point_on_surf2,s2_copy,ray_start)
	
	s3,_ = snellLaw(s2,surf_N2,n1,n2,Print=Print)
	s3_copy = s3.getA().T[0,:]
	surf_N3_copy = surf_N3.getA().T[0,:]
	x_y_z2 = LineSurfIntersec(surf_N3_copy,point_on_surf3,s3_copy,x_y_z1)

	s4,_ = snellLaw(s3,surf_N3,n2,n_air,Print=Print)
	ray1 = np.arccos(np.dot(s4.T,-surf_N))*180/np.pi # deviation of ray wrt flat surface! instead of surface of final prism.
	x_y_offset = x_y_z2

	return s4,ray1.getA(),x_y_offset


def CompoundDoubletMirror(prism1,prism2,Lambda,phi=0,offset=0.,Print=False,shift=0,vec_offset=None):
	'''
	This is a mirror function of compoundDoublet and will calculate the deviation of a compound prism made of two materials, i.e. a doublet.

	The function will calculate the ray tracing through vectors as this will be the easiest method for including 
	any rotations necessary.
	
	Inputs will include a lambda, the prism dictionaries,phi and offset (from the optical axis). phi and offset must be in radians
	Output will be a ray vector (s4) and deviation in degrees relitive to a surface perpendicular to the optical axis (ray1).
	'''

	
	#prism dimensions
	T_edge1 = float(prism1['T(mm)'])
	T_edge2 = float(prism2['T(mm)'])
	triangle_dims1 = (float(prism1['Diameter']), float(prism1['B']), float(prism1['C']))
	triangle_dims2 = (float(prism2['Diameter']), float(prism2['B']), float(prism2['C']))
	# need the angles in radians #
	alpha1 = float(prism1['Wedge angle'])*np.pi/180.
	alpha2 = float(prism2['Wedge angle'])*np.pi/180.
	
	beta2 = (alpha1 - alpha2)
	
	
	#all rotations
	rotx1 = RotMat3Dx(-beta2)
	rotx2 = RotMat3Dx(alpha1)
	rotx3 = RotMat3Dx(-alpha2)
	rotz  = RotMat3Dz(phi)
	

	#all surface calculations using co-ords sys of +y is north +z is east
	#rays travel in +z direction
	C = triangle_dims1[0]
	A1 = triangle_dims1[1]
	A2 = triangle_dims2[1]
	d = C/2. * np.tan(alpha2) +T_edge2
	x = C/2. * np.tan(alpha1)/np.cos(alpha1)+T_edge1/np.cos(beta2)# for the case of a compound prism.
	ray_start = np.array([0,0,0]) + shift #generally in x,y,z
	point_on_surf2 = np.array([0,0,x])
	point_on_surf3 = np.array([0,0,x+d])
	
	#first surface
	surf_N = np.matrix([0.,0.,-1.]).T
	surf_N1_ = rotx1*surf_N
	surf_N1 = rotz*surf_N1_
	#print surf_N

	#second surface
	surf_N2_ = rotx2*surf_N1_
	surf_N2 = rotz*surf_N2_
	#print surf_N2
	
	#Third Surface - exit surface of compound prism
	surf_N3_ =  rotx3*surf_N2_
	surf_N3 = surf_N3_
	#print surf_N3
	

	#going to test to make sure deviation is correct here
	# raytrace for first lambda of the bandwidth
	n1 = eval(prism1['Material']+'({})'.format(Lambda))
	n2 = eval(prism2['Material']+'({})'.format(Lambda))
	
	
	s1 = vec_offset
	
	#beam = 9.0 #mm
	#incident = 12.5 #mm
	
	# need to make copies of the surfaces so values are not overwritten as well as transpose so that they fit the function
	s2,_ = snellLaw(s1,surf_N1,n_air,n1,Print=Print)
	s2_copy = s2.getA().T[0,:]
	surf_N2_copy = surf_N2.getA().T[0,:]
	x_y_z1 = LineSurfIntersec(surf_N2_copy,point_on_surf2,s2_copy,ray_start)
	
	s3,_ = snellLaw(s2,surf_N2,n1,n2,Print=Print)
	s3_copy = s3.getA().T[0,:]
	surf_N3_copy = surf_N3.getA().T[0,:]
	x_y_z2 = LineSurfIntersec(surf_N3_copy,point_on_surf3,s3_copy,x_y_z1)

	s4,_ = snellLaw(s3,surf_N3,n2,n_air,Print=Print)
	ray1 = np.arccos(np.dot(s4.T,-surf_N))*180/np.pi # deviation of ray wrt flat surface! instead of surface of final prism.
	x_y_offset = x_y_z2

	return s4,ray1.getA(),x_y_offset


def DoubleDoublet(prisms,Lambda,Bandwidth,n,Plot=True,zenith=np.array([0,10,20,30,40,50]),rotate=np.linspace(0,180,181)):
	'''
	Runs through the entire double doublet model calculating the dispersion, angular offset and lateral displacement
	'''
	l0 = Lambda
	bw = Bandwidth
	l1 = l0 - bw/2.
	l2 = l0 + bw/2.
	mid = (n-1)/2
	wavelenghts = np.linspace(l1,l2,n)
	prism1 = prisms[0]
	prism2 = prisms[1]
	print 'Combination of {} and {} glass with wedge angles of {} and {}'.format(prism1['Material'],prism2['Material'],prism1['Wedge angle'],prism2['Wedge angle'])
	for z in zenith:
		offsets = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,wavelenghts))*462
		offsets =  offsets[mid] - offsets  
		disp4phi=[]
		disp4phi2=[]
		for phi in rotate:
			phi = phi*np.pi/180. 
			Vectors=[]
			RayAngles=[]
			XYoffsets=[]
			for Lambda,Offset in zip(wavelenghts,offsets):
				vector,rayAng,xyoffset = PF.CompoundDoublet(prism1,prism2,Lambda,phi=phi,offset=Offset)
				Vectors.append(vector)
				RayAngles.append(rayAng[0,0])
				XYoffsets.append(xyoffset)
			maxindex = np.argmax(RayAngles)
			minindex = np.argmin(RayAngles)
			dispersion = RayAngles[maxindex] - RayAngles[minindex]
			dispersion_ = np.arccos(np.dot(Vectors[maxindex].T,Vectors[minindex]))*180/np.pi
			disp4phi.append(dispersion_[0,0]*3600) #dispersion in arcseconds
		
		disp4phi = np.array(disp4phi)
		if Plot == True:
			plt.plot(rotate,disp4phi,label=str(z))
		optRot.append(rotate[np.argmin(disp4phi)])
	if Plot == True:
		plt.plot([0,deg],[limit/2,limit/2],'--k',label='lim')
		plt.legend(loc ='best')
		plt.xlabel('Rotation of prisms with respect to one another (deg)')
		plt.ylabel('Dispersion (arcseconds)')
		plt.title('Combination of {} and {} glass with wedge angles of {} and {}'.format(prism1['Material'],prism2['Material'],prism1['Wedge angle'],prism2['Wedge angle']))
		plt.axis('tight')
		plt.show()
		plt.close()

	deg=180
	rotate = np.linspace(0,deg,1*deg+1) 
	OptRot=[]
	for z,optphi in zip(zenith,optRot):
		
		offsets = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,wavelenghts))*462
		offsets =  (offsets[mid] - offsets)
		disp4phi=[]
		
		for phi in rotate:
			optphi = optphi * np.pi/180.
			phi =  optphi - phi*np.pi/180. 
			Vectors=[]
			RayAngles=[]
			XYoffsets=[]
			for Lambda,Offset in zip(wavelenghts,offsets):
				vector,rayAng,xyoffset = PF.CompoundDoublet(prism1,prism2,Lambda,phi=optphi,offset=Offset)
				vector,rayAng,xyoffset = PF.CompoundDoubletMirror(prism2,prism1,Lambda,phi=phi,vec_offset=vector)
				Vectors.append(vector)
				RayAngles.append(rayAng[0,0])
				XYoffsets.append(xyoffset)
			
			maxindex = np.argmax(RayAngles)
			minindex = np.argmin(RayAngles)
			dispersion = RayAngles[maxindex] - RayAngles[minindex]
			dispersion_ = np.arccos(np.dot(Vectors[maxindex].T,Vectors[minindex]))*180/np.pi
			
			disp4phi.append(dispersion_[0,0]*3600) #dispersion in arcseconds
			
		disp4phi = np.array(disp4phi)
		if Plot == True:
			plt.plot(rotate,disp4phi,label=str(z))
		
		OptRot.append(rotate[np.argmin(disp4phi)])
		
	
	if Plot == True:
		plt.plot([0,deg],[limit,limit],'--k',label='lim')
		plt.legend(loc ='best')
		plt.xlabel('Rotation of prisms with respect to one another (deg)')
		plt.ylabel('Dispersion (arcseconds)')
		plt.title('Combination of {} and {} glass with wedge angles of {} and {}'.format(prism1['Material'],prism2['Material'],prism1['Wedge angle'],prism2['Wedge angle']))
		plt.axis('tight')
		plt.show()
	
	return disp4phi 
