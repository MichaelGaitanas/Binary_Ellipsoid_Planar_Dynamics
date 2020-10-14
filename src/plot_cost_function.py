import numpy as np
import matplotlib.pyplot as plt

#load the inputs and store them in an array
inputs = np.loadtxt('../resources/sys_inputs.txt')

#ellipsoid 1 axes
a1 = inputs[3]
b1 = inputs[4]
c1 = inputs[5]

#ellipsoid 2 axes
a2 = inputs[6]
b2 = inputs[7]
c2 = inputs[8]

costFunc = np.loadtxt('../resources/cost_function_data.txt') #(r,vthita,J) data set
equilibr = np.loadtxt('../resources/Scheeres_equilibrium_curve_data.txt') #analytical curve vthita = f(r)

################################################################################

r = costFunc[:,0]
vthita = costFunc[:,1]
J = costFunc[:,2]

"""plt.figure(1)
graph = plt.scatter(r,vthita, c = np.log10(J), marker = ',', cmap = plt.cm.coolwarm, s = 500)
cb = plt.colorbar(graph)
cb.set_label('J')
plt.xlabel('r')
plt.ylabel('vthita')
plt.savefig('cost_func_colormap.png')"""

################################################################################

#We choose the 'plotPoints' minimum values of J and save the corresponding pairs (r,vthita) into j[][]
j = [ [] , [] ]
plotPoints = 60
for i in range(plotPoints):
    minPos = np.where(J == np.amin(J))
    j[0].append(r[minPos[0][0]])
    j[1].append(vthita[minPos[0][0]])
    r = np.delete(r,minPos[0][0])
    vthita = np.delete(vthita,minPos[0][0])
    J = np.delete(J,minPos[0][0])

#plot the Scheeres equilibrium curve and the data j[][]
plt.figure(2)
plt.title('Cost function J(t,vthita) fitting')
plt.xlabel('r [km]')
plt.ylabel('thitaDot [rad/sec]')
plt.gcf().subplots_adjust(left = 0.15)
plt.plot(equilibr[:,0], equilibr[:,1])
plt.scatter(j[0],j[1], color = (1,0,0), alpha = 0.5, label = '| J (r,vthita) | < tolerance')
plt.text(max(r)/1.08, max(vthita)/0.95, 'a1 = ' + str(a1) + ' [km]\nb1 = ' + str(b1) + ' [km]\nc1 = ' + str(c1) + ' [km]')
plt.text(max(r)/1.08, max(vthita)/1.1, 'a2 = ' + str(a2) + ' [km]\nb2 = ' + str(b2) + ' [km]\nc2 = ' + str(c2) + ' [km]')
plt.legend( scatterpoints = 1, loc = 'upper right', fontsize = 8)
plt.savefig('../plots/fitting_J(r,vthita).png')
plt.show()
