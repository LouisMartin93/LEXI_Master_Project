'''
Simple validation such that my prism ray tracing with rotations is valid (values are compared to zemax values)
'''


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

#validate the prism model with zemax first will check for N-BK7 and CaF2 at wedge angles of 1,5,10 and rotation angles of prisms of 0,180 


#first specify the prism characteristics
#NBK7
def NBK7(lam):
	B1,B2,B3 = (1.03961212,0.231792344,1.01046945)
	C1,C2,C3 = (0.00600069867,0.0200179144,103.560653)
	n = np.sqrt(1 + (B1*lam**2)/(lam**2 - C1) +(B2*lam**2)/(lam**2 - C2)+(B3*lam**2)/(lam**2 - C3))
	return np.round(n,4)

def UVFS(lam):
	B1,B2,B3 = (0.6961663,0.4079426,0.8974794)
	C1,C2,C3 = (0.0684043**2,0.1162414**2,9.896161**2)
	n = np.sqrt(1 + (B1*lam**2)/(lam**2 - C1) +(B2*lam**2)/(lam**2 - C2)+(B3*lam**2)/(lam**2 - C3))
	return np.round(n,4)

def MgF2(lam):
	B1,B2,B3 = (0.48755108,0.39875031,2.310353)
	C1,C2,C3 = (0.04338408,0.09461442,23.793604)
	n = np.sqrt(1 + (B1*lam**2)/(lam**2 - C1**2) +(B2*lam**2)/(lam**2 - C2**2)+(B3*lam**2)/(lam**2 - C3**2))
	return np.round(n,4)



mat1 = 'NBK7'
deg1 = 1.0
mat2 = 'MgF2'
deg2 = 3.0
prism1={}
prism1['Material'] = mat1
prism1['Wedge angle'] = deg1
prism1['Diameter'] = 25.
prism1['C'] = 10
prism1['B'] = 10
print prism1['Material']
prism2={}
prism2['Material'] = mat2
prism2['Wedge angle'] = deg2
prism2['Diameter'] = 25.
prism2['C'] = 10
prism2['B'] = 10
wlen1 = 0.6
wlen2 = 0.8


# this is for an zenith of 60. from the zenith... so of course it does not calculate the dispersion correctly
offsetb = 0.000748
offsetr = -0.000480
phi = np.linspace(0,180,)
phi = np.array([0,15,30,45,60,90,120,150,180])
mm =[]
for ang in phi:

	ray1 = PF.CorrectionSingle(prism1,prism2,wlen1,ang,offset=offsetb)
	ref = PF.CorrectionSingle(prism1,prism2,0.7,ang)
	ray2 = PF.CorrectionSingle(prism1,prism2,wlen2,ang,offset=offsetr)
	print ang
	print ray1,ref,ray2
	print np.max([ray1, ref,ray2])-np.min([ray1, ref,ray2])
	
	mm.append(ref)
zm=[10.548,10.4571,10.17956,9.7256,9.099,7.3971,5.2045,2.6888,0]
zm_rot = [0,15,30,45,60,90,120,150,180]
rotation = phi
#print rotation
zemax_model = np.array(zm)
my_model = np.array(mm)[:,0,0]
#print zemax_model
#print my_model
perc = (zemax_model-my_model)/zemax_model *100.
#print perc
fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax2 = ax1.twinx()
lns1 = ax1.plot(zm_rot,zemax_model,label='zemax')
lns2 = ax1.plot(rotation,my_model,label='my model')
#lns3 = ax2.plot(rotation,perc,'r-',label='% diff')
ax1.set_xlabel('Rotation of prisms wrt each other (deg)')
ax1.set_ylabel(r'Angle deviation at 0.7 $\mu m$ (deg)')
#ax2.set_ylabel('Percentage off set (%)')
lns =lns1+lns2#+lns3
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc='best')
plt.title('Validation for prisms {} and {} glass with wedge angles of {} and {}'.format(prism1['Material'],prism2['Material'],prism1['Wedge angle'],prism2['Wedge angle']))
plt.show()

