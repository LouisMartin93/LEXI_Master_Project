import numpy as np
import matplotlib.pyplot as plt

z = np.linspace(0,np.pi/3.,100)

a_vals = [0.56,0.5,0.4,0.3,0.2,0.1]

for a in a_vals:
	func =  np.arccos(-a*np.tan(z))*180./np.pi
	plt.plot(z*180./np.pi, func-90.,label='a={}'.format(a))

plt.legend(loc='best')
plt.xlabel('zenith (deg)')
plt.ylabel('rotation (deg)')
plt.show()



