#contains funcitons
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.patches import Rectangle
path = '/home/martin/Desktop/LEXI ADC project/atmospheric model/'
sys.path.insert(1,path)
import AtmosModel as am

# calculating wedge angle given a deviation angle.
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


#----- formula------#
def WedgeAngle(n,D):
	D = D *np.pi/180.
	beta = 90. - np.arctan( (n - np.cos(D))/np.sin(D))*180./np.pi
	return beta

#-----working wavelngth----#
lambda0 = 0.633
T = 3.
Diameter = 25.4
H = 'n/a'
Manu = 'Sydor'
Prism = 'Wedge'
dev_angles = np.linspace(0.5,7.0,66)
#print dev_angles


from astropy.table import Table
from astropy.io import ascii

t = Table(names=('Material','nd', "abbe number", 'T(mm)', 'Wedge angle', 'Diameter', 'T_edge', 'wedge radians', 'B', 'C', 'H', 'Manufacturer,',"Prism Type"), dtype=('S5', 'f8', 'f8','f8','f8','f8','f8','f8','f8','f8','S5','S10','S10'))
print t

#--wedge angles for NBK7--#
print'Wedge angles for NBK7'
for angle in dev_angles:
	beta = angle#WedgeAngle(BaF2(lambda0),angle)
	material = 'NBK7'
	nd = NBK7(0.5876)
	abbe_num = 81.78
	T_edge = T + Diameter*np.tan(beta*np.pi/180)
	wedge_rad = beta*np.pi/180
	B = T_edge - T
	C = np.sqrt(Diameter**2. + B**2.)
	t.add_row([material, nd, abbe_num, T, np.round(beta,4), Diameter, T_edge, wedge_rad, B, C, H,Manu, Prism])

#--wedge angles for BaF2--#
print'Wedge angles for BAF2'
for angle in dev_angles:
	beta = angle#WedgeAngle(BaF2(lambda0),angle)
	material = 'BaF2'
	nd = BaF2(0.5876)
	abbe_num = 81.78
	T_edge = T + Diameter*np.tan(beta*np.pi/180)
	wedge_rad = beta*np.pi/180
	B = T_edge - T
	C = np.sqrt(Diameter**2. + B**2.)
	t.add_row([material, nd, abbe_num, T, np.round(beta,4), Diameter, T_edge, wedge_rad, B, C, H,Manu, Prism])


#--wedge angles for UVFS--#
print'Wedge angles for UVFS'
for angle in dev_angles:
	beta = angle#WedgeAngle(BaF2(lambda0),angle)
	material = 'UVFS'
	nd = UVFS(0.5876)
	abbe_num = 81.78
	T_edge = T + Diameter*np.tan(beta*np.pi/180)
	wedge_rad = beta*np.pi/180
	B = T_edge - T
	C = np.sqrt(Diameter**2. + B**2.)
	t.add_row([material, nd, abbe_num, T, np.round(beta,4), Diameter, T_edge, wedge_rad, B, C, H,Manu, Prism])


#--wedge angles for CaF2--#
print'Wedge angles for CAF2'

for angle in dev_angles:
	beta = angle #WedgeAngle(CaF2(lambda0),angle)
	material = 'CaF2'
	nd =CaF2(0.5876)
	abbe_num = 94.99
	T_edge = T + Diameter*np.tan(beta*np.pi/180)
	wedge_rad = beta*np.pi/180
	B = T_edge - T
	C = np.sqrt(Diameter**2. + B**2.)
	
	t.add_row([material, nd, abbe_num, T, np.round(beta,4), Diameter, T_edge, wedge_rad, B, C, H,Manu, Prism])


#--wedge angles for MgF2--#
print'Wedge angles for MGF2'
for angle in dev_angles:
	beta = angle#WedgeAngle(MgF2(lambda0),angle)
	material = 'MgF2'
	nd = MgF2(0.5876)
	abbe_num = 106.22
	T_edge = T + Diameter*np.tan(beta*np.pi/180)
	wedge_rad = beta*np.pi/180
	B = T_edge - T
	C = np.sqrt(Diameter**2. + B**2.)
	
	t.add_row([material, nd, abbe_num, T, np.round(beta,4), Diameter, T_edge, wedge_rad, B, C, H,Manu, Prism])




t = Table(names=('Material','nd', "abbe number", 'T(mm)', 'Wedge angle', 'Diameter', 'T_edge', 'wedge radians', 'B', 'C', 'H', 'Manufacturer,',"Prism Type"), dtype=('S5', 'f8', 'f8','f8','f8','f8','f8','f8','f8','f8','S5','S10','S10'))
filename = 'OptimalPrisms.fits'

#N-SF11 prism characteristics
T = 5.0
angle = 3.0
material = 'NSF11'
nd = NSF11(0.5876)
abbe_num = 26
T_edge = T + Diameter*np.tan(angle*np.pi/180)
wedge_rad = angle*np.pi/180
B = T_edge - T
C = np.sqrt(Diameter**2. + B**2.)
Manu = 'Perkins'
Prism ='Wedge'
t.add_row([material, nd, abbe_num, T, np.round(angle,4), Diameter, T_edge, wedge_rad, B, C, H,Manu, Prism])

#CaF2 wedge prism Characteristics
T = 5.0
angle = 5.37
material = 'CaF2'
nd = CaF2(0.5876)
abbe_num = 94.99
T_edge = T + Diameter*np.tan(angle*np.pi/180)
wedge_rad = angle*np.pi/180
B = T_edge - T
C = np.sqrt(Diameter**2. + B**2.)
Manu = 'Perkins'
Prism ='Wedge'
t.add_row([material, nd, abbe_num, T, np.round(angle,4), Diameter, T_edge, wedge_rad, B, C, H,Manu, Prism])

t.write(filename,format='fits',overwrite=True)
