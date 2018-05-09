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

def OptRotZen(prism1,prism2,fname,model,num=60,zen=[0]):
    l0= 0.7
    bw =0.2
    l1 = l0 - bw/2.
    l2 = l0 + bw/2.
    optimal_vals =[]
    rotate = np.linspace(0,180,601)
    if len(zen) == 1:
        zenith = np.linspace(0,60,num+1)
    else:
        zenith = zen

    for z in zenith:
	    disp4phi =[]
	    dev_ref = AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l0)
	    for phi in rotate:
		    #for the smallest wavelenght
		    offset1 = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l1)-dev_ref)*462
		    _,ray1,off = eval('{}(prism1,prism2,l1,phi*np.pi/180.,offset=offset1)'.format(model))
		    
		    #for the central wavelength
		    _,ray_ref,off = eval('{}(prism1,prism2,l0,phi*np.pi/180.)'.format(model))   
		    #for the longest wavelength
		    offset2 = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l2)-dev_ref)*462
		    _,ray2,off = eval('{}(prism1,prism2,l2,phi*np.pi/180.,offset=offset2)'.format(model)) 

		    #calculating the dispersion
		    dispersion = np.max([ray1,ray2,ray_ref])- np.min([ray1,ray2,ray_ref])
		    disp4phi.append(dispersion*3600)
	    
	    disp4phi = np.array(disp4phi) 
	    index, = np.where(disp4phi == np.min(disp4phi))
	    
	    if len(index) >1:
	        index = index[0]
	    optimal_vals.append(rotate[index])
	    
    
    fits.writeto(fname+'.fits',np.array(zip(zenith,optimal_vals)),overwrite=True)
    return np.array(zip(zenith,optimal_vals))

def DispersionCalculation(zenith,rotate,prism1,prism2,model):
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
        _,ray1,xyoff = eval('{}(prism1,prism2,l1,phi*np.pi/180.,offset=offset1)'.format(model)) #PF.CorrectionSingle(prism1,prism2,l1,phi,offset=offset1)
		        
        #for the central wavelength
        _,ray_ref,xyoff = eval('{}(prism1,prism2,l0,phi*np.pi/180.)'.format(model))#PF.CorrectionSingle(prism1,prism2,l0,phi)
		        
        #for the longest wavelength
        offset2 = (AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l2)-dev_ref)*462
        _,ray2,xyoff = eval('{}(prism1,prism2,l2,phi*np.pi/180.,offset=offset2)'.format(model)) #PF.CorrectionSingle(prism1,prism2,l2,phi,offset=offset2)

        #calculating the dispersion
        dispersion = np.max([ray1,ray2,ray_ref])- np.min([ray1,ray2,ray_ref])
        disp4phi.append(dispersion*3600)
	    
    disp4phi = np.array(disp4phi) 
    return disp4phi


def AxisOffset(prisms,fname,model,Read=True,Plot=True,zen=[0]):
    # this will calculate the offsets for angluar deviations as well as the offsets physically of the beam exiting the system.

    #loading the prism data.
    prism1 = prisms[0]
    prism2 = prisms[1]
    l0= 0.7
    bw =0.2
    l1 = l0 - bw/2.
    l2 = l0 + bw/2.
    if Read == True:
        data = fits.getdata(fname+'.fits')
    else:
        data = OptRotZen(prism1,prism2,fname,model,zen=zen)
    
    rotate = data[:,1]
    zenith = data[:,0]
    x_vals = []
    y_vals = []
    rotxyvar =[]
    ###------calculates the offsets and makes a plot showing the variation for a given zenith----#
    offsets =[] #note the pural difference
    zero_mat = np.matrix([0.0]) # needed otherwise stuck in matrix format...
    for z, phi in zip(zenith,rotate):
        
        dev_ref = AM.delR(AM.h0,z*np.pi/180.,AM.a,AM.b,l0)
        # for the central wavelength
        vec,ray,xyoff = eval('{}(prism1,prism2,l0,phi*np.pi/180)'.format(model))
        offset = np.max([ray,zero_mat])
        offsets.append(offset) # offset in degrees
        rotxyvar.append(xyoff)
        x_vals.append(np.max([vec[0],zero_mat]))
        y_vals.append(np.min([vec[1],zero_mat]))
        

    # assuming ray is initally centered along axis, using central wavelenght of the band so will be.
    offsets = np.array(offsets)
    plt.subplot(1,2,1)
    plt.plot(zenith,offsets) 
    plt.xlabel('Zenith (deg)',size=16)
    plt.ylabel('Offset angle (deg)',size=16)
    plt.title('Offset for rotation angle needed at given zenith',size=18)
    ## x-y positioning
    plt.subplot(1,2,2)
    x_vals =np.array(x_vals)*180/np.pi
    y_vals =np.array(y_vals)*180/np.pi
    plt.plot(x_vals,y_vals,'.',lw=2)
    plt.plot([-np.max(abs(x_vals)),np.max(abs(x_vals))],[0,0],'--k')
    plt.plot([0,0],[-np.max(abs(x_vals)),np.max(abs(x_vals))],'--k')
    plt.xlabel('x co-ord (deg)',size=16)
    plt.ylabel('y co-ord (deg)',size=16)
    plt.title('The Drift in beam direction as the ADC rotates about x-y cordinate.',size=18)
    plt.axis([-np.max(abs(x_vals)),np.max(abs(x_vals)),-np.max(abs(x_vals)),np.max(abs(x_vals))])
    plt.grid()
    plt.show()
	
    plt.plot(zenith,(180.-rotate))
    plt.xlabel('Zenith (deg)')
    plt.ylabel('rotation angle required (deg)')
    plt.title('Optimal rotation angle needed for a given zenith (180 is zeroth position)')
    plt.show()
    
    prism_radius = prism1['Diameter']/2.
    beam_rad = 9./2.
    theta = np.linspace(0,2*np.pi)
    circ_x = prism_radius * np.cos(theta)
    circ_y = prism_radius * np.sin(theta)
    beam_x = beam_rad * np.cos(theta)
    beam_y = beam_rad * np.sin(theta) 	
    rotxyvar = np.array(rotxyvar)
    plt.plot(beam_x+np.max(rotxyvar[:,0]),beam_y+np.max(rotxyvar[:,1]),'-r',label='Beam exit')
    plt.plot(beam_x,beam_y,'-g',label='Beam enter')
    plt.plot(circ_x,circ_y,'-k',label='Prism Exit surface')
    plt.xlabel('x co-ord (mm)',size=16)
    plt.ylabel('y co-ord (mm)',size=16)
    plt.title('Lateral shift at exit pupil due to Dispersion',size=18)
    plt.legend(loc='best')
    plt.grid()
    plt.show()


def RotationModel(prisms,fname,polyOrd,model,Read=True,Plot=True):
	#loading the prism data.
    prism1 = prisms[0]
    prism2 = prisms[1]
    l0= 0.7
    bw =0.2
    l1 = l0 - bw/2.
    l2 = l0 + bw/2.
    if Read == True:
        data = fits.getdata(fname+'.fits')
    else:
        data = OptRotZen(prism1,prism2,fname,model)
    
    rotate = data[:,1]
    zenith = data[:,0]
    #---- Model section ----#
    # making my fit model for this particular combination
    index, = np.where(rotate != 180)
    ind = index[0] # so there are no zero values
    y = rotate[ind:]
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
        disp4phi = DispersionCalculation(zenith,rotate,prism1,prism2,model)
        plt.plot(zenith,disp4phi,label='Recorded values')

        #looking at the dispersion after the prisms correct for the atmosphere (model values)
        rotate =  p(zenith)
        disp4phi = DispersionCalculation(zenith,rotate,prism1,prism2,model)
        plt.plot(zenith,disp4phi,label='Predicted Values')

        plt.legend(loc='best')
        plt.xlabel('Zenith angle')
        plt.ylabel('Dispersion after correction (arcseconds)')
        plt.title('Dispersion after correction at given Zenith angles using model')
        plt.show()
    
    return co_effs


