import numpy as np
import scipy as sp
from matplotlib import pyplot as plt


###--------------Atmoshperic Model Functions---------------###
###--------------------------------------------------------###


#--------------------------------------------------------#
#------------define the constants and inputs-------------#
#--------------------------------------------------------#

Lambda = 0.633 #microns
re = 6371.0*10**3 # mean radius of the earth (m)
h0 = 2324 #altitude of observatory on la palma (m)
Mma = 0.0289644 # mean molar mass of dry air (kg/mol)
Lapse_rate = 6.49/1000. # lapse rate of  atmosphere (K/m)
g = 9.81 #acceleration due to gravity (m/s^2)
T0 = 4+273.15 # mean temperature at 00:00 in November at WHT (K)
R = 8.31447 # universal gas constant (J/(mol*K))
D = 1. # 4.26 #telescope diameter.


K0 = 238.0185 #(um^-2) 
K1 = 0.05792105 #(um^-2)
K2 = 57.362 #(um^-2)
K3 = 0.00167917 #(um^-2)

#-------------------------------------------------------#
#------------------equation calculation-----------------#
#-------------------------------------------------------#

q = (g*Mma/(R*Lapse_rate))
H0 = T0/Lapse_rate
a = (q-1.)/H0
b = (q-1.)*(q-2.)/(2.*H0**2)
r0 = re + h0


#-------c(lambda)--------#

def c(Lambda):
	return (K1/(K0 - (1/Lambda)**2.) + K3/(K2 - (1/Lambda)**2.))

#------A function----#

def A(a,b,h,phi):
	return (a*r0 - b*r0**2. * np.cos(phi)**2.)*(np.sqrt(2.*h/r0 +np.cos(phi)**2.)- np.cos(phi))

#------B function-------#

def B(b,h,phi):
	return (b*r0**2 /3.)*((2*h/r0 +np.cos(phi)**2.)**1.5 - np.cos(phi)**3.)

#---------dR function----------#

def delR(h,phi,a,b,Lambda):
	'''
	h in meters
	lambda in microns
	phi in radians
	'''
	return c(Lambda)*np.sin(phi)*(A(a,b,h,phi) + B(b,h,phi))*np.pi

#-------Refraction function---------#

def Refraction(Lambda0,Bandwidth,z0,D,Plot=True,num=100):
	'''
	Inputs Lambda0 given as float in microns and is assumed to be the central wavelength.
	Bandwidth is the full bandwidth also given in microns.
	z0 is the zenith angle given in degrees 
	D is the telescope diameter.
	'''
	#phi = np.linspace(0,z0,(num*z0)+1.)*np.pi/180.
	phi = z0*np.pi/180.
	print phi
	factor = .1 #for spectrocopy
	#factor = 1./30. #for coronagraphs
	limit =   factor* Lambda0*10**-6 / D *180./np.pi *3600. #arcsec

	L1 = Lambda0-Bandwidth/2.
	L2 = Lambda0+Bandwidth/2.
	print delR(h0,phi,a,b,Lambda0)*206265.
	R_L1 = delR(h0,phi,a,b,L1)
	print R_L1*206265.
	R_L2 = delR(h0,phi,a,b,L2)
	print R_L2*206265.
	diff_R = R_L1 - R_L2
	if Plot == True:
		plt.plot(phi*180/np.pi,diff_R*206265.,label=str(Bandwidth*1000)+'$nm$')
		plt.plot([phi[0],phi[-1]*180/np.pi],[limit,limit],'k--')
	return diff_R
