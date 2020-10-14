import numpy as np
import matplotlib.pyplot as plt

#load the inputs and store them in an array
inputs = np.loadtxt('../resources/sys_inputs.txt')

G = inputs[0]
r = inputs[11]
Torb = inputs[9]
thitaDot = 2.0*np.pi/Torb

a2 = inputs[6]
b2 = inputs[7]
c2 = inputs[8]
I2x = (b2**2 + c2**2)/5.0
I2y = (a2**2 + c2**2)/5.0
I2z = (a2**2 + b2**2)/5.0

V1 = 4.0*np.pi*inputs[3]*inputs[4]*inputs[5]/3.0 # V_ellipsoid = 4*pi*a*b*c/3

#construct a 2D grid with x-axis -> a1/b1 and y-axis b1/c1
x = np.linspace(1.0,1.5,50)
y = np.linspace(1.0,1.5,50)
a1b1, b1c1 = np.meshgrid(x,y)

b1 = (3*b1c1*V1/(4*np.pi*a1b1))**(1.0/3.0)
c1 = b1/b1c1
a1 = b1*a1b1
I1x = (b1**2 + c1**2)/5.0
I1y = (a1**2 + c1**2)/5.0
I1z = (a1**2 + b1**2)/5.0

def inertia_multiplier(r):
    trace1 = I1x + I1y + I1z
    trace2 = I2x + I2y + I2z
    return 1.0 + (3.0/(2.0*r**2))*(trace1 + trace2 - \
                 (3.0/2.0)*( I1x + I1y - (I1y - I1x) + I2x + I2y - (I2y - I2x) ) )

#corrected mass estimation
M = (thitaDot**2)*(r**3)/(G*inertia_multiplier(r))
#export the data into a file
fpMass = open('../resources/mass_from_varying_primary_axes.txt','w')
for i in range(len(x)):
    for j in range(len(y)):
        fpMass.write('%f %f %f\n'%(a1b1[i][j], b1c1[i][j], M[i][j]) )
fpMass.close()

#plot the mass of the binary M as a function of a1/b1, b1/c1 in a colormap
plt.figure(1)
graph = plt.scatter(a1b1,b1c1, c = M, marker = ',', cmap = plt.cm.coolwarm, s = 50)
cb = plt.colorbar(graph)
cb.set_label('M [kg]')
plt.xlabel('a1/b1')
plt.ylabel('b1/c1')
plt.title('M_DRA = 5.37e11 [kg]')
plt.savefig('../plots/mass_correction1.png')
