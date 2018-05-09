import numpy as np
import scipy as sp
from matplotlib import pyplot as plt





#--------------------------------------------------------#
#------------------define the constants------------------#
#--------------------------------------------------------#

re = 6371.0*10**3 # mean radius of the earth (m)
Mma = 0.0289644 # mean molar mass of dry air (kg/mol)
Lapse_rate = 6.49/1000. # lapse rate of  atmosphere (K/m)
g = 9.81 #acceleration due to gravity (m/s^2)
R = 8.31447 # universal gas constant (J/(mol*K))


###------------------Atmoshperic Models--------------------###
###--------------------------------------------------------###
class SphericallySymmetric_model():
	'''
	The following model is not the same as the concentric spherical shell model and instead follows derivations from the paper 'Atmospheric Refraction 
	Predictions based on Acutal Atomspheric Pressure and Temperature Data' (Nauenberg 20117) https://arxiv.org/abs/1609.08921.
	'''
	
		
	#-------c(lambda)--------#
	@classmethod
	def c(self,Lambda):
		#This is the reflectivity based on Ciddor 1996
		K0 = 238.0185 #(um^-2) 
		K1 = 0.05792105 #(um^-2)
		K2 = 57.362 #(um^-2)
		K3 = 0.00167917 #(um^-2)
		return (K1/(K0 - (1/Lambda)**2.) + K3/(K2 - (1/Lambda)**2.))
	#------A function----#
	@classmethod
	def A(self,a,b,h,phi):
		r0 = re + h
		return (a*r0 - b*r0**2. * np.cos(phi)**2.)*(np.sqrt(2.*h/r0 +np.cos(phi)**2.)- np.cos(phi))

	#------B function-------#
	@classmethod
	def B(self,b,h,phi):
		r0 = re + h
		return (b*r0**2 /3.)*((2*h/r0 +np.cos(phi)**2.)**1.5 - np.cos(phi)**3.)

	#---------dR function----------#
	@classmethod
	def delR(self,h,phi,a,b,Lambda):
		'''
		h in meters
		lambda in microns
		phi in radians
		'''
		return self.c(Lambda)*np.sin(phi)*(self.A(a,b,h,phi) + self.B(b,h,phi))#*np.pi
	
	#-------Refraction function---------#
	def Refraction(self,Lambda,z,h,T):
		'''
		Inputs Lambda given as float in microns.
		z is the zenith angle given in radians 
		h is the heght at which observing
		T is the temperature at the observation given in K
		P is the pressure at the obervation given in atmos
		'''
		
		q = (g*Mma/(R*Lapse_rate))
		H0 = T/Lapse_rate
		a = (q-1.)/H0
		b = (q-1.)*(q-2.)/(2.*H0**2)
		r0 = re + h
		dR = self.delR(h,z,a,b,Lambda)
		return dR


	# this is the correct way to calculate the dispersion... will have to comeback to this later 
	def quick_test(self,wl,zen):
		return np.log(1 + self.c(wl))*np.tan(zen)

class ParallelPlane_model():
    ''' 
    This is the Parallel Plane model used to calculate atmospheric disperions it is calculated
    using the local enviromental parameters. The equation works as refraction = (n-1)tan(z).
    To calculate n, calculations are followed from Ciddor 1996. The refraction is returned in radians.
    '''

    @classmethod
    def n_ax(self,Lambda):
        # refractive index of light due to dry air
        K0 = 238.0185 #(um^-2) 
        K1 = 0.05792105 #(um^-2)
        K2 = 57.362 #(um^-2)
        K3 = 0.00167917 #(um^-2)
        return 1 + (K1/(K0 - (1/Lambda)**2.) + K3/(K2 - (1/Lambda)**2.))

    @classmethod
    def n_ws(self,Lambda):
        # this is the refraction of light in air due to water vapour
        #constants
        w0 = 295.235
        w1 = 2.6422
        w2 = -0.032380
        w3 = 0.004028
        return 1 + 1.022*10**-8 *(w0 + w1*Lambda**-2 + w2*Lambda**-4 +w3*Lambda**-6)

    @classmethod
    def MA(self,carbon):
        #carbon is the CO2 concentraion in ppm
        # in kg/mol is the molar mass of dry air
        return 10**-3 * (28.9635 + 12.011*10**-6 *(carbon-400))
    
	

    @classmethod
    def K2C(self,T):
        #convert kelvin to Celcius
        return T-273.15
    
    @classmethod
    def svp(self,T):
        #saturation vapour pressure of water vapour measured in Pa

        # constants
        A = 1.2378847*10**-5 #K^-2
        B = -1.9121316*10**-2 #K^-1
        C = 33.93711047
        D = -6.3431645*10**3 #K

        return np.exp(A*T**2 + B*T + C +D/T)


    @classmethod
    def f(self,p,T):
        # enhancement factor of water vapour in air
        t = self.K2C(T)
        #constants
        alpha = 1.00062
        beta = 3.14*10**-8 #Pa^-1
        gamma = 5.6*10**-7 #degC^-1

        return alpha + beta*p + gamma*t**2

    @classmethod
    def XW(self,RH,p,T):
        #molar fraction of water vapour in moist air
        return RH *self.f(p,T)*self.svp(T)/p


    @classmethod
    def Zcompress(self,p,T,RH):
        #compresability

        # constants
        a0 = 1.58123*10**-6 # K Pa^-1
        a1 = -2.9331*10**-8 # Pa^-1
        a2 = 1.1043*10**-10 # (KPa) ^-1
        b0 = 5.707*10**-6 # K Pa^-1
        b1 = -2.051*10**-8 # Pa ^-1
        c0 = 1.9898*10**-4 # K Pa^-1
        c1 = -2.376*10**-6 # Pa^-1
        d = 1.83*10**-11 # K^2 Pa^-2
        e = -0.765*10**-8 # K^2 Pa^-2
        # temperature in celcius
        t = self.K2C(T)
        #water vapour content
        xw = self.XW(RH,p,T)

        return 1 - (p/T)*(a0 + a1*t +a2*t**2 + (b0 +b1*t)*xw + (c0+c1*t)*xw**2) + (p/T)**2 *(d + e*xw**2)

    @classmethod
    def BIPM(self,p,T,RH,carbon,method='gen'):
        # atmospheric air density calculation given the fractional water vapour and carbon contents
        Z = self.Zcompress(p,T,RH)
        xw = self.XW(RH,p,T)
        Ma = self.MA(carbon)
        Mw = 0.018015 #kg/mol
        if method=='gen':
            rho = p*Ma/(Z*R*T) * (1 - xw*(1 - Mw/Ma))
        elif method=='water':
            rho = p*Mw*xw/(Z*R*T)
        elif method=='dry':
            rho = p*Ma*(1-xw)/(Z*R*T)
        return rho


    def refraction(self,T,p,RH,wl,xc,zen):
        ''' 
        #-----------Input parameters-------#
        T: temperature given in kelvins (K)
        p: pressure given in Pascals (Pa)
        RH: Relative Humidity 0-1
        wl: wavelenght given in microns (um)
        xc: carbon content in dry air given as parts per million (ppm)
        zen: zenith angle of observer given in degrees

        #-------------Output--------------#
        The Refraction given in radians.
        '''

        #convert degrees to radians
        zen = zen*np.pi/180.

        # calculating densities for standard conditions (scaling values)
        rho_axs = self.BIPM(101325,288.15,0,450)
        rho_ws = self.BIPM(1333,293.15,1,450)

        rho_a = self.BIPM(p,T,RH,450,method='dry') 
        rho_w = self.BIPM(p,T,1,450,method='water')

        refractivity = (rho_a/rho_axs)*(self.n_ax(wl)-1.)+ (rho_w/rho_ws)*(self.n_ws(wl)-1.)
        refraction = refractivity*np.tan(zen) # remember in radians
        return refraction





