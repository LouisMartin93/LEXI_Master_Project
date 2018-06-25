###-----------------------------------------###
###----- Robust Dispersion Calculation -----###
###-----------------------------------------###


###------ import functions & libraries -----###
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.optimize import fmin
from scipy.optimize import minimize
from astropy.table import Table
import astropy.io.fits as fits
path = '/home/martin/Desktop/LEXI ADC project/atmospheric model/' # this dierctory needs to change to where ever the location of the atmospheric model folder is
sys.path.insert(1,path)
path = '/home/martin/Desktop/LEXI ADC project/prism model/Scripts' # this dierctory needs to change to where ever the location of the prism model folder is
sys.path.insert(1,path)

import AtmosModel as AM
import PrismFunctions as PF



def AtmosphericDispersion(Lambda0,Bandwidth,EnviroParam):
	'''
	#---------- Inputs -----------#
	central wavelength (in micron)
	full bandwidth (in micron)
	Enviromental parameters (T,P,RH,xc,zenith) in that order Temp in kelvin,P in pascals, Relative humidity 0-1,xc in ppm, zenith in degrees 
	   
	#---------- Outputs ----------#
	array of length(wavelengths) of the angular offsets in radians about the central wavelength.
	
	### ------------------------###
	calculates the dispersion created by the atmosphere for a particular bandwidth and atmospheric conditions using the parallel plane model. 
	'''
	atmosphere = AM.ParallelPlane_model()
	T,P,RH,xc,zenith = EnviroParam
	l1 = Lambda0 - Bandwidth/2.
	l2 = Lambda0 + Bandwidth/2.
	n = 11
	mid = (n-1)/2
	wavelengths = np.linspace(l1,l2,n)
	
	offsets = (atmosphere.refraction(T,P,RH,wavelengths,xc,zenith))
	offsets =  (offsets[mid] - offsets)

	#Returns the offsets relative the central wavelength in radians and the corresponding wavelengths in microns.	
	return offsets,wavelengths



def Dispersion4rotation(RotAng,Prisms,wavelengths,offsets):
	'''
	#-------------- Inputs ---------# 
	rotation angle for compound prisms (in degrees)
	a table of two prisms with all their required dimension data.
	array of wavelengths
	array of offsets of rays due to atmosphere relative to central wavelenght.
	
	#------------- Outputs ---------#
	The dispersion in Arcseconds.

	###---------------------------###
	Calculates the 3-D ray propagation through a double doublet for given wedge angles and materials and for a given rotation of the prisms w.r.t one another.
	'''

	#rotation angle converted to radians
	RotAng =  RotAng*np.pi/180.
	#empty lists for outputs of prism solver
	Vectors=[]
	RayAngles=[]
	XYoffsets=[]

	#defining the two prisms from the table input.
	prism1 = Prisms[0]
	prism2 = Prisms[1]
		
	
	for Lambda,Offset in zip(wavelengths,offsets):
		vector,rayAng,xyoffset = PF.CompoundDoublet(prism1, prism2, Lambda, phi=RotAng, offset=Offset)
		vector,rayAng,xyoffset = PF.CompoundDoubletMirror(prism2, prism1, Lambda, phi=-RotAng, vec_offset=vector, shift=xyoffset)
		Vectors.append(vector)
		RayAngles.append(rayAng[0,0])
		XYoffsets.append(xyoffset)
		
	maxindex = np.argmax(RayAngles)
	minindex = np.argmin(RayAngles)
	
	dispersion_ = np.arccos(np.dot(Vectors[maxindex].T,Vectors[minindex]))*180/np.pi*3600 # dispersion in arcseconds
	
	# this is only outputting the dispersion in arcseconds	
	return dispersion_[0,0]

#------- Solver function -----------#

def SolveOptimumRotation(Prisms,EnviroParams,Lambda,Bandwidth):
	'''
	###------- Inputs ------ ###
	A table of two prisms with all their required dimension data.
	Enviromental parameters (T,P,RH,xc,zenith) in that order Temp in kelvin,P in pascals, Relative humidity 0-1,xc in ppm, zenith in degrees
	Central wavelength (in micron)
	Full bandwidth (in micron)
	
	###------ Outputs -----###
	The rotation needed for each prism to minimise given the zenith and enviromental parameters.
	The value of the dispersion at this rotation.

	###-------------------###
	Using a powell minimsiation the dispersion created by the atmosphere (in the y axis) is minimised for the system. 
	'''
	offsets,wavelengths = AtmosphericDispersion(Lambda0,Bandwidth,EnviroParams)
	res = minimize(Dispersion4rotation,int(0), args=(Prisms,wavelengths,offsets) ,method='Powell')
	 
	return res.x, res.fun






