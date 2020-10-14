import numpy as np
from scipy.integrate import solve_ivp
import time

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

#boolean variable that states if the 2 bodies collided during the simulation
global collisionFlag
collisionFlag = False

################################################################################

#binary's kinetic energy
def kinetic_energy(vec):
    return 0.5*I1z*vec[6]**2 + 0.5*I2z*vec[7]**2 + 0.5*m*vec[4]**2 + \
           0.5*(I1z + I2z + m*vec[0]**2)*vec[5]**2 + (I1z*vec[6] + I2z*vec[7])*vec[5]

#binary's potential energy
def potential_energy(vec):
    trace1 = (I1x + I1y + I1z)/M1
    trace2 = (I2x + I2y + I2z)/M2
    return (-G*M1*M2/vec[0])*(1.0 + (1.0/(2.0*vec[0]**2))*( trace1 + trace2 - \
           (3.0/2.0)*( (I1x + I1y - np.cos(2.0*vec[2])*(I1y - I1x))/M1 + \
                       (I2x + I2y - np.cos(2.0*vec[3])*(I2y - I2x))/M2 ) ) )

#binary's angular momentum
def angular_momentum(vec):
    return (I1z + I2z + m*vec[0]**2)*vec[5] + I1z*vec[6] + I2z*vec[7]

#total energy and angular momentum of the binary
def constants_of_motion(vec):
    T = kinetic_energy(vec)
    V = potential_energy(vec)
    E = T + V
    L = angular_momentum(vec) #Lz
    return E,L

#check if collision happened during the simulation
def check_for_collsion(t, collisionFlag):
    if collisionFlag:
        print('Collision detection at t = %.4f [sec]'%(t))

#calculate energy and angular momentum from the resulted solution
#and prin them into files
def post_process_solution(sol):
    E,L = constants_of_motion(sol.y)
    fpElements = open('../resources/orbital_elements.txt','w')
    fpConstants = open('../resources/constants_of_motion.txt','w')
    for i in range(len(sol.t)):
        fpElements.write('%lf %lf %lf %lf %lf\n'%(sol.t[i],sol.y[0][i],sol.y[1][i],sol.y[2][i],sol.y[3][i]))
        fpConstants.write('%lf %.18lf %.18lf\n'%(sol.t[i],E[i],L[i]))
    fpElements.close()
    fpConstants.close()
    
################################################################################

#potential's r-derivative
def dV_dr(r,phi1,phi2):
    grav = G*M1*M2/r**2
    trace1 = (I1x + I1y + I1z)/M1
    trace2 = (I2x + I2y + I2z)/M2
    term1 = (I1x + I1y - np.cos(2.0*phi1)*(I1y - I1x))/M1
    term2 = (I2x + I2y - np.cos(2.0*phi2)*(I2y - I2x))/M2
    return grav*(1.0 + (3.0/(2.0*r**2))*( trace1 + trace2 - (3.0/2.0)*(term1 + term2) ) )

#potential's phi1-derivative
def dV_dphi1(r,phi1):
    return 3.0*G*M2*np.sin(2.0*phi1)*(I1y - I1x)/(2.0*r**3)

#potential's phi2-derivative
def dV_dphi2(r,phi2):
    return 3.0*G*M1*np.sin(2.0*phi2)*(I2y - I2x)/(2.0*r**3)

#8 odes in total
def odes(t,vec):
    #extrac the vec[] components into human readable symbols for simplicity
    r = vec[0]
    thita = vec[1]
    phi1 = vec[2]
    phi2 = vec[3]
    vr = vec[4]
    vthita = vec[5]
    vphi1 = vec[6]
    vphi2 = vec[7]

    #approximate collision detection between the 2 bodies
    #assuming a circumscribed sphere for each body, the radius of which is equal
    #to the longest semi axis a (a1 for body 1 and a2 for body 2)
    global collisionFlag
    if r <= a1 + a2:
        collisionFlag = True
        return

    #4 odes
    dr_dt = vr
    dthita_dt = vthita
    dphi1_dt = vphi1
    dphi2_dt = vphi2

    #and the rest 4 odes
    dvr_dt = r*vthita**2 - dV_dr(r,phi1,phi2)/m
    dvthita_dt = dV_dphi1(r,phi1)/(m*r**2) + dV_dphi2(r,phi2)/(m*r**2) - 2.0*vr*vthita/r
    dvphi1_dt = -(1.0 + m*r**2/I1z)*dV_dphi1(r,phi1)/(m*r**2) - dV_dphi2(r,phi2)/(m*r**2) + 2.0*vr*vthita/r
    dvphi2_dt = -(1.0 + m*r**2/I2z)*dV_dphi2(r,phi2)/(m*r**2) - dV_dphi1(r,phi1)/(m*r**2) + 2.0*vr*vthita/r

    return [dr_dt,dthita_dt,dphi1_dt,dphi2_dt, dvr_dt,dvthita_dt,dvphi1_dt,dvphi2_dt]

################################################################################

def main():
    #initial conditions
    r = inputs[11];       vr = inputs[15];
    thita = inputs[12];   vthita = inputs[16];
    phi1 = inputs[13];    vphi1 = inputs[17];
    phi2 = inputs[14];    vphi2 = inputs[18];

    #time parameters
    t0 = inputs[19]
    tmax = inputs[20]
    dt = inputs[21]

    print('Calculating orbits...')
    tStart = time.time() #start measuring cpu time
    
    initConditions = [r,thita,phi1,phi2, vr,vthita,vphi1,vphi2]
    tSpan = np.array([t0, tmax])
    tEval = np.arange(t0, tmax + dt, dt)
    sol = solve_ivp(odes, tSpan, initConditions, method = 'DOP853', t_eval = tEval, max_step = dt)
    print('Post-processing solution...')
    check_for_collsion(max(sol.t)+dt, collisionFlag)
    post_process_solution(sol)
    
    tEnd = time.time() #end measuring cpu time
    cpuTime = tEnd - tStart
    print('Done. Estimated completion time %f [sec]'%(cpuTime))

main()
