import numpy as np
import matplotlib.pyplot as plt

#load the orbital elements
orb = np.loadtxt('../resources/orbital_elements.txt')

print('Plotting r(t)')
plt.figure(1)
plt.xlabel('t [sec]')
plt.ylabel('r [km]')
plt.gcf().subplots_adjust(left = 0.15)
plt.plot(orb[:,0],orb[:,1])
plt.savefig('../plots/r(t).png')

print('Plotting thita(t)')
plt.figure(2)
plt.xlabel('t [sec]')
plt.ylabel('thita [rad]')
plt.plot(orb[:,0],orb[:,2])
plt.savefig('../plots/thita(t).png')

print('Plotting phi(t)')
plt.figure(3)
plt.xlabel('t [sec]')
plt.ylabel('phi1 [rad]')
plt.plot(orb[:,0],orb[:,3])
plt.savefig('../plots/phi1(t).png')

print('Plotting phi2(t)')
plt.figure(4)
plt.xlabel('t [sec]')
plt.ylabel('phi2 [rad]')
plt.gcf().subplots_adjust(left = 0.15)
plt.plot(orb[:,0],orb[:,4])
plt.savefig('../plots/phi2(t).png')

################################################################################

#load the energy and the angular momentum of the binary
consts = np.loadtxt('../resources/constants_of_motion.txt')

print('Plotting E(t)')
plt.figure(5)
plt.xlabel('t [sec]')
plt.ylabel('E [kg*km^2/sec^2]')
plt.plot(consts[:,0],consts[:,1])
plt.savefig('../plots/E(t).png')

print('Plotting Lz(t)')
plt.figure(6)
plt.xlabel('t [sec]')
plt.ylabel('Lz [kg*km^2/sec]')
plt.plot(consts[:,0],consts[:,2])
plt.savefig('../plots/Lz(t).png')

plt.show()
