import numpy as np
import matplotlib.pyplot as plt

#load the inputs and store them in an array
inputs = np.loadtxt('../resources/sys_inputs.txt')

G = inputs[0]
r = inputs[11]
Torb = inputs[9]
thitaDot = 2.0*np.pi/Torb

a1 = inputs[3]
b1 = inputs[4]
c1 = inputs[5]
I1x = (b1**2 + c1**2)/5.0
I1y = (a1**2 + c1**2)/5.0
I1z = (a1**2 + b1**2)/5.0

#V_ellipsoid = 4*pi*a*b*c/3
V2 = 4.0*np.pi*inputs[6]*inputs[7]*inputs[8]/3.0

#construct a 2D grid with x-axis -> a2/b2 and y-axis b2/c2
x = np.linspace(1.0,1.5,50)
y = np.linspace(1.0,1.5,50)
a2b2, b2c2 = np.meshgrid(x,y)

b2 = (3*b2c2*V2/(4*np.pi*a2b2))**(1.0/3.0)
c2 = b2/b2c2
a2 = b2*a2b2
#NORMALIZED moments of inertia
I2x = (b2**2 + c2**2)/5.0
I2y = (a2**2 + c2**2)/5.0
I2z = (a2**2 + b2**2)/5.0

def inertia_multiplier(r):
    trace1 = I1x + I1y + I1z
    trace2 = I2x + I2y + I2z
    return 1.0 + (3.0/(2.0*r**2))*(trace1 + trace2 - \
                 (3.0/2.0)*( I1x + I1y - (I1y - I1x) + I2x + I2y - (I2y - I2x) ) )

#corrected mass estimation
M = (thitaDot**2)*(r**3)/(G*inertia_multiplier(r))
#export the data into a file
fpMass = open('../resources/mass_from_varying_secondary_axes.txt','w')
for i in range(len(x)):
    for j in range(len(y)):
        fpMass.write('%f %f %f\n'%(a2b2[i][j], b2c2[i][j], M[i][j]) )
fpMass.close()

#plot the mass of the binary M as a function of a2/b2, b2/c2 in a colormap
plt.figure(1)
graph = plt.scatter(a2b2,b2c2, c = M, marker = ',', cmap = plt.cm.coolwarm, s = 50)
cb = plt.colorbar(graph)
cb.set_label('M [kg]')
plt.xlabel('a2/b2')
plt.ylabel('b2/c2')
plt.title('M_DRA = 5.37e11 [kg]')
plt.savefig('../plots/mass_correction2.png')
