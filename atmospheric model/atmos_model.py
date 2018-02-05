import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
#--------------------------------------------------------#
#------------define the constants and inputs-------------#
#--------------------------------------------------------#

Lambda = 0.633 #microns
re = 6371.0*10**3 # mean radius of the earth (m)
h0 = 2350. #4200 #altitude of observatory on la palma (m)
Mma = 0.0289644 # mean molar mass of dry air (kg/mol)
Lapse_rate = 6.5/1000. # 7.91/1000. #  lapse rate of  atmosphere (K/m)
g = 9.81 #acceleration due to gravity (m/s^2)
T0 = 4+273.15 # 319.3 # mean temperature at 00:00 in November at WHT (K)
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
	return c(Lambda)*np.sin(phi)*(A(a,b,h,phi) + B(b,h,phi))*np.pi

#-------Refraction function---------#

def Refraction(Lambda0,Bandwidth,z0,D,Plot=True,num=100):
	phi = np.linspace(0,z0,num*z0+1)*np.pi/180
	factor = .1 #for spectrocopy
	#factor = 1./30. #for coronagraphs
	limit =   factor* Lambda0*10**-6 / D *180/np.pi *3600 #arcsec

	L1 = Lambda0-Bandwidth/2.
	L2 = Lambda0+Bandwidth/2.
	lambda_range = np.linspace(L1,L2, num=5)
	R_L1 = delR(h0,phi,a,b,L1)
	R_L2 = delR(h0,phi,a,b,L2)
	diff_R = R_L1 - R_L2
	if Plot == True:
		plt.plot(phi*180/np.pi,diff_R*206265.,label=str(Bandwidth*1000)+'$nm$')
		plt.plot([phi[0],phi[-1]*180/np.pi],[limit,limit],'k--')
	return diff_R


#print q, H0/1000, a, b, r0/1000, c(Lambda)
wlen2 = 0.8
wlen1 = 0.6
z0 = 60

def calc_refraction(wlen,h0,z0):
	h = np.linspace(h0,110000,1000)
	n = 1. + c(wlen) *(1.+a*h +b*(h**2))
	n0 = n[0]
	zn = np.arcsin(n0*h0*np.sin(z0)/(n*(h)))
	R = delR(h,zn,a,b,wlen)
	deltaR = np.diff(R)
	sumR = np.sum(np.abs(deltaR))
	Rtot1 =  (R[0] + sumR) *180/np.pi *3600
	return Rtot1

def calc_dispersion(wlen1,wlen2,h0,z0):
	'''
	calculates dispersion of atmosphere between 2 wavelengths in microns at a given hight looking at a zenith angle in radians.
	Returned dispersion is in arcseconds.
	'''
	#calculates tropospheres values
	ht =11000
	hs = 100000
	ht_vals = np.linspace(h0,ht,1000)
	hs_vals = np.linspace(ht,hs,1000)

	n = 1. + c(wlen1) *(1.+a*ht_vals +b*(ht_vals**2))
	n0 = n[0]
	nt = n[-1]
	zt_vals = np.arcsin(n0*h0*np.sin(z0)/(n*(ht_vals)))
	zt = zt_vals[-1]
	#calculates stratosphere values
	n = 1. + c(wlen1) *(1.+a*hs_vals +b*(hs_vals**2))
	zs_vals = np.arcsin(nt*ht*np.sin(zt)/(n*(hs_vals)))
	#join arrays
	h = np.concatenate([ht_vals,hs_vals])
	zn = np.concatenate([zt_vals,zs_vals])
	# refraction for wlen1
	R = delR(h,zn,a,b,wlen1)
	deltaR = np.diff(R)
	sumR = np.sum((deltaR))
	Rtot1 =  (R[0]) *180/np.pi *3600

	n = 1. + c(wlen2) *(1.+a*ht_vals +b*(ht_vals**2))
	n0 = n[0]
	nt = n[-1]
	zt_vals = np.arcsin(n0*h0*np.sin(z0)/(n*(ht_vals)))
	zt = zt_vals[-1]
	#calculates stratosphere values
	n = 1. + c(wlen1) *(1.+a*hs_vals +b*(hs_vals**2))
	zs_vals = np.arcsin(nt*ht*np.sin(zt)/(n*(hs_vals)))
	#join arrays
	h = np.concatenate([ht_vals,hs_vals])
	zn = np.concatenate([zt_vals,zs_vals])
	# refraction for wlen2
	R = delR(h,zn,a,b,wlen2)
	deltaR = np.diff(R)
	sumR = np.sum((deltaR))
	Rtot2 =  (R[0]) *180/np.pi *3600
	return Rtot1,Rtot2


for z in [0,10,20,30,40,45,50,60]:
	z0 = z*np.pi/180.
	R1,R2 = delR(h0,z0,a,b,wlen1),delR(h0,z0,a,b,wlen2)
	print z, (R1-R2)*180./np.pi *3600.

Lambda = [0.6,0.7,0.8,0.9,1.0]
phi = np.linspace(0,z0,100*z0+1)
#plot for various Lambda fixed h and L
for i in range(len(Lambda)):
	diff_R = delR(h0,phi,a,b,Lambda[i])
	plt.plot(phi*180/np.pi,diff_R*206265.,label=str(Lambda[i])+'$\mu m$')
	print c(Lambda[i])
	print diff_R[60]*206265.
	print '\n'

plt.ylabel(r'Refraction contribtuion,$\delta R$ (arcsec)',size=17)
plt.xlabel(r'Angle from the zenith,$\phi$ (degrees)',size=17)
plt.legend(loc=2)
plt.title('Comparison of refraction due to the atmosphere for varying angle from zenith',size=19)
plt.show()
plt.close()


lambda0 = 0.7
BW = [0.05,0.1,0.2]
for i in range(len(BW)):
	z0 = 60
	diff_R = Refraction(lambda0,BW[i],z0,D)
	phi = np.linspace(0,z0,100*z0+1)
	limit = 0.1* lambda0*10**-6 / D *180/np.pi *3600 #arcsec
	#correction at 60 deg from zenith
	print(diff_R[-1]*206265 - limit)
	closeness = np.abs(np.array(diff_R*206265.) - limit)
	min_index = np.argmin(closeness)
	max_zenith = phi[min_index]
	print max_zenith
	

plt.ylabel(r'Refraction Difference,$\Delta R$ (arcsec)',size=17)
plt.xlabel(r'Angle from the zenith,$\phi$ (degrees)',size=17)
plt.legend(loc=2)
plt.title('Refraction difference for given bandwidth at '+str(lambda0)+'$\mu m$ vs angle from zenith',size=19)
plt.show()
plt.close()


lambda0 = 0.7
BW = [200./1000.]
for i in range(len(BW)):
	z0 = 60
	diff_R = Refraction(lambda0,BW[i],z0,D)
	phi = np.linspace(0,z0,100*z0+1)
	limit = 0.1* lambda0*10**-6 / D *180/np.pi *3600 #arcsec
	#correction at 60 deg from zenith
	print(diff_R[-1]*206265 - limit)
	closeness = np.abs(np.array(diff_R*206265.) - limit)
	min_index = np.argmin(closeness)
	max_zenith = phi[min_index]
	print max_zenith	

plt.ylabel(r'Refraction Difference,$\Delta R$ (arcsec)',size=17)
plt.xlabel(r'Angle from the zenith,$\phi$ (degrees)',size=17)
plt.legend(loc=2)
plt.title('Refraction difference for Y bandwidth vs angle from zenith',size=19)
plt.show()
plt.close()


"""
bandwidths = np.linspace(1,400,200)
for j in range(len(Lambda)):
	zeniths=[]
	for i in range(len(bandwidths)):
		z0 = 90
		diff_R = Refraction(Lambda[j],bandwidths[i]/1000.,z0,D,Plot=False)
		phi = np.linspace(0,z0,100*z0+1)#*np.pi/180
		closeness = np.abs(np.array(diff_R*206265.) - limit)
		min_index = np.argmin(closeness)
		max_zenith = phi[min_index]
		zeniths.append(max_zenith)

	zeniths = np.array(zeniths)
	plt.plot(bandwidths,zeniths,label=r'$\lambda_{0}$:'+r'{}$\mu m$'.format(Lambda[j]))
plt.xlabel(r'$\Delta \lambda$ full bandwidth (nm)')
plt.ylabel('Zenith angle in degrees')
plt.legend(loc='best')
plt.show()
"""

