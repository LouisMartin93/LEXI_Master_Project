'''
This is a library of functions for robustly calculating a high-order polynomial function for two singlets given such that if a zenith angles is given the out put is the correct angle of rotation to insure that the dispersion is minimised to the criteria given in the calcualtions. 
'''


#contains funcitons
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.patches import Rectangle
path = '/home/martin/Desktop/LEXI ADC project/atmospheric model/'
sys.path.insert(1,path)

path = '/home/martin/Desktop/LEXI ADC project/prism model/'
sys.path.insert(1,path)

import AtmosModel as AM
import PrismFunctions as PF
import time
import astropy.io.fits as fits

Lambda = 0.7 #microns
re = 6371.0*10**3 # mean radius of the earth (m)
h0 = 2324 #altitude of observatory on la palma (m)
Mma = 0.0289644 # mean molar mass of dry air (kg/mol)
Lapse_rate = 6.49/1000. # lapse rate of  atmosphere (K/m)
g = 9.81 #acceleration due to gravity (m/s^2)
T0 = 4+273.15 # mean temperature at 00:00 in November at WHT (K)
R = 8.31447 # universal gas constant (J/(mol*K))
D = 1. # 4.26 #telescope diameter.

q = (g*Mma/(R*Lapse_rate))
H0 = T0/Lapse_rate
a = (q-1.)/H0
b = (q-1.)*(q-2.)/(2.*H0**2)
r0 = re + h0

def OptRotZen(prism1,prism2,fname):
    l0= 0.7
    bw =0.2
    l1 = l0 - bw/2.
    l2 = l0 + bw/2.
    optimal_vals =[]
    rotate = np.linspace(110,180,601)
    zenith = np.linspace(0,60,61)
    for z in zenith:
	    disp4phi =[]
	    dev_ref = AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l0)
	    for phi in rotate:
		    #for the smallest wavelenght
		    offset1 = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l1)-dev_ref)*462
		    ray1 = PF.CorrectionSingle(prism1,prism2,l1,phi,offset=offset1)
		    
		    #for the central wavelength
		    ray_ref = PF.CorrectionSingle(prism1,prism2,l0,phi)
		    
		    #for the longest wavelength
		    offset2 = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l2)-dev_ref)*462
		    ray2 = PF.CorrectionSingle(prism1,prism2,l2,phi,offset=offset2)

		    #calculating the dispersion
		    dispersion = np.max([ray1,ray2,ray_ref])- np.min([ray1,ray2,ray_ref])
		    disp4phi.append(dispersion*3600)
	    
	    disp4phi = np.array(disp4phi) 
	    index, = np.where(disp4phi == np.min(disp4phi))
	    optimal_vals.append(rotate[index])
	    

    fits.writeto(fname+'.fits',np.array(zip(zenith,optimal_vals)),overwrite=True)
    return np.array(zip(zenith,optimal_vals))

def DispersionCalculation(zenith,rotate,prism1,prism2):
    #looking at the dispersion after the prisms correct for the atmosphere (measured values)
    disp4phi =[]
    l0= 0.7
    bw =0.2
    l1 = l0 - bw/2.
    l2 = l0 + bw/2.
    for z, phi in zip(zenith,rotate):
        dev_ref = AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l0)
        #for the smallest wavelenght
        offset1 = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l1)-dev_ref)*462
        ray1 = PF.CorrectionSingle(prism1,prism2,l1,phi,offset=offset1)
		        
        #for the central wavelength
        ray_ref = PF.CorrectionSingle(prism1,prism2,l0,phi)
		        
        #for the longest wavelength
        offset2 = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l2)-dev_ref)*462
        ray2 = PF.CorrectionSingle(prism1,prism2,l2,phi,offset=offset2)

        #calculating the dispersion
        dispersion = np.max([ray1,ray2,ray_ref])- np.min([ray1,ray2,ray_ref])
        disp4phi.append(dispersion*3600)
	    
    disp4phi = np.array(disp4phi) 
    return disp4phi


def AxisOffsetAndModel(prisms,fname,polyOrd,Read=True,Plot=True):


    #loading the prism data.
    from astropy.io import ascii
    table = ascii.read('../glass tables/Glass_worksheet1.txt')
    glass_by_shape = table.group_by('Prism Type')
    prisms_same = glass_by_shape.groups[0]
    prism1 = prisms_same[prisms[0]]
    prism2 = prisms_same[prisms[1]]
    l0= 0.7
    bw =0.2
    l1 = l0 - bw/2.
    l2 = l0 + bw/2.
    if Read == True:
        data = fits.getdata(fname+'.fits')
    else:
        data = OptRotZen(prism1,prism2,fname)
    
    rotate = data[:,1]
    zenith = data[:,0]

    ###------calculates the offsets and makes a plot showing the variation for a given zenith----#
    offsets =[] #note the pural difference
    zero_mat = np.matrix([0.0]) # needed otherwise stuck in matrix format...
    for z, phi in zip(zenith,rotate):
        dev_ref = AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l0)
        # for the central wavelength
        ray = PF.CorrectionSingle(prism1,prism2,l0,phi)
        offset = np.max([ray,zero_mat])
        offsets.append(offset) # offset in degrees

    # assuming ray is initally centered along axis, using central wavelenght of the band so will be.
    offsets = np.array(offsets)
    plt.subplot(1,2,1)
    plt.plot(zenith,offsets) 
    plt.xlabel('Zenith (deg)')
    plt.ylabel('Offset angle (deg)')
    plt.title('Offset for rotation angle needed at given zenith')
    plt.subplot(1,2,2)
    plt.plot(zenith,(180-rotate))
    plt.xlabel('Zenith (deg)')
    plt.ylabel('rotation angle required (deg)')
    plt.title('Optimal rotation angle needed for a given zenith (180 is zeroth position)')
    plt.show()


    #---- Model section ----#
    # making my fit model for this particular combination
    index, = np.where(rotate != 180)
    ind = index[0] # so there are no zero values
    y = 180 - rotate[ind:]
    x =  zenith[ind:]
    co_effs = np.polyfit(x, y, polyOrd)
    p = np.poly1d(co_effs)

    xp = np.linspace(0, 60, 100)
    plt.plot(x, y, '.', label='measured')
    plt.plot(xp, p(xp), '-',label='model')
    plt.xlabel('Zenith (deg)')
    plt.ylabel('Rotation needed')
    plt.title('Function for required dispersion minimsation, poly order: {}'.format(polyOrd))
    plt.show()

    diff =[]
    for i in range(len(x)):
        #print p(x[i]),y[i], p(x[i])-y[i]
        diff.append(p(x[i])-y[i])

    abs_diff = np.abs(diff)
    if Plot:
        #looking at the dispersion after the prisms correct for the atmosphere (measured values) 
        disp4phi = DispersionCalculation(zenith,rotate,prism1,prism2)
        plt.plot(zenith,disp4phi,label='Recorded values')

        #looking at the dispersion after the prisms correct for the atmosphere (model values)
        rotate = 180 - p(zenith)
        disp4phi = DispersionCalculation(zenith,rotate,prism1,prism2)
        plt.plot(zenith,disp4phi,label='Predicted Values')

        plt.legend(loc='best')
        plt.xlabel('Zenith angle')
        plt.ylabel('Dispersion after correction (arcseconds)')
        plt.title('Dispersion after correction at given Zenith angles using model')
        plt.show()
    
    return co_effs


