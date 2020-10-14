import numpy as np
import matplotlib.pyplot as plt

#load the inputs and store them in an array
inputs = np.loadtxt('../resources/sys_inputs.txt')

#gravitational constant
G = inputs[0]

#mass parameters
M1 = inputs[1]
M2 = inputs[2]
m = M1*M2/(M1 + M2)

#ellipsoid 1 semi axes and moments of inertia
a1 = inputs[3]
b1 = inputs[4]
c1 = inputs[5]
I1x = M1*(b1**2 + c1**2)/5.0
I1y = M1*(a1**2 + c1**2)/5.0
I1z = M1*(a1**2 + b1**2)/5.0

#ellipsoid 2 semi axes and moments of inertia
a2 = inputs[6]
b2 = inputs[7]
c2 = inputs[8]
I2x = M2*(b2**2 + c2**2)/5.0
I2y = M2*(a2**2 + c2**2)/5.0
I2z = M2*(a2**2 + b2**2)/5.0

#equilibrium relationship thita_dot = f(r)
def thita_dot(r):
    grav = G*(M1 + M2)/r**3
    trace1 = (I1x + I1y + I1z)/M1
    trace2 = (I2x + I2y + I2z)/M2
    return np.sqrt(grav*( 1.0 + (3.0/(2.0*r**2))*(trace1 + trace2 -
                  (3.0/2.0)*( (I1x + I1y - (I1y - I1x))/M1 + \
                              (I2x + I2y - (I2y - I2x))/M2 )   ) ) )

def main():
    #domain of calculation
    r0 = 1.0; rmax = 1.5; dr = 0.01;
    r = np.arange(r0,rmax,dr)
    thitaDot = thita_dot(r) #calculate the function thitaDot = f(r)

    #write the data (r,thitaDot(r)) into a file
    fpCurve = open('../resources/Scheeres_equilibrium_curve_data.txt','w')
    for i in range(len(thitaDot)):
        fpCurve.write('%lf %.15lf\n'%(r[i],thitaDot[i]))
    fpCurve.close()

    #plot the curve
    plt.figure(1)
    plt.title('Scheeres equilibrium solution r - thitaDot')
    plt.xlabel('r [km]')
    plt.ylabel('thitaDot [rad/sec]')
    plt.gcf().subplots_adjust(left = 0.15)
    plt.text(max(r)/1.15, max(thitaDot)/1.1,
            'a1 = ' + str(a1) + ' [km]\nb1 = ' + str(b1) + ' [km]\nc1 = ' + str(c1) + ' [km]')
    plt.text(max(r)/1.15, max(thitaDot)/1.25,
            'a2 = ' + str(a2) + ' [km]\nb2 = ' + str(b2) + ' [km]\nc2 = ' + str(c2) + ' [km]')
    plt.plot(r,thitaDot)
    plt.savefig('../plots/Shceeres_thitaDot_r.png')
    plt.show()

main()
