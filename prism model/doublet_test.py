# this will be used as a test calculation for the doublet model using the prisms as well as functions from the PrismFuncitons

#contains funcitons
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.patches import Rectangle
path = '/home/martin/Desktop/LEXI ADC project/atmospheric model/'
sys.path.insert(1,path)
import AtmosModel as AM
import PrismFunctions as PF
import time

# general loading of data.
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


#loading the prism data.
from astropy.io import ascii
table = ascii.read('glass tables/Glass_worksheet1.txt')
glass_by_shape = table.group_by('Prism Type')
prisms_same = glass_by_shape.groups[0]
limit = 1./10.* Lambda*10**-6 / D *180/np.pi *3600 *462#arcsec 
print 'limit is {} arcseconds'.format(np.round(limit,4))
rotate = np.linspace(0,180,181)
zenith = np.linspace(0,60,7)

prism1 = prisms_same[38]
prism2 = prisms_same[100]
l0= 0.7
bw =0.2
l1 = l0 - bw/2.
l2 = l0 + bw/2.
print 'Combination of {} and {} glass with wedge angles of {} and {}'.format(prism1['Material'],prism2['Material'],prism1['Wedge angle'],prism2['Wedge angle'])

#new design of double doublet
PF.DoubleDoublet2(prism1,prism2,Lambda,offset=0.)

'''
# this is for looking up at zenith.
# Initial double singlets
Lambda = 0.7 #microns
ray_ref = PF.CorrectionSingle(prism1,prism2,Lambda,180,offset=0.)
Lambda = 0.8 #microns
ray_long = PF.CorrectionSingle(prism1,prism2,Lambda,180,offset=0.)
Lambda = 0.6 #microns
ray_small = PF.CorrectionSingle(prism1,prism2,Lambda,180,offset=0.)
print ray_ref,ray_long,ray_small
print np.max([ray_ref,ray_long,ray_small]) - np.min([ray_ref,ray_long,ray_small])

# Single compound prism
Lambda = 0.7 #microns
vec_ref,ray_ref = PF.CompoundDoublet(prism1,prism2,Lambda,offset=0.)
Lambda = 0.8 #microns
vec_long,ray_long = PF.CompoundDoublet(prism1,prism2,Lambda,offset=0.)
Lambda = 0.6 #microns
vec_small,ray_small = PF.CompoundDoublet(prism1,prism2,Lambda,offset=0.)
print ray_ref,ray_long,ray_small
print np.max([ray_ref,ray_long,ray_small]) - np.min([ray_ref,ray_long,ray_small])

# Double Doublet
Lambda = 0.7 #microns
vec_ref,ray_ref = PF.DoubleDoublet(prism1,prism2,Lambda,offset=0.)
Lambda = 0.8 #microns
vec_long,ray_long = PF.DoubleDoublet(prism1,prism2,Lambda,offset=0.)
Lambda = 0.6 #microns
vec_small,ray_small = PF.DoubleDoublet(prism1,prism2,Lambda,offset=0.)
print ray_ref,ray_long,ray_small
print np.max([ray_ref,ray_long,ray_small]) - np.min([ray_ref,ray_long,ray_small])
'''
