# Binary_Ellipsoid_Planar_Dynamics


============================= Project description =============================

The planar ellipsoid - ellispoid full 2 body problem is being studied in this
directory. In a nutshell, you have 2 rigid triaxial ellipsoids in space, each having its
own physical charasteristics (semi-axes mass, and moments of inertia). You "throw" them in space with
a specific set of initial conditions (positions, linear velocities,
orientations and angular velocities) and let the binary evolve according to the
laws of Newtonian mechanics. Keep in mind that the differential equations are such
that the motion of the 2 ellipsoids is bounded to a constant plane.

The parent directory contains the following directories:


==================================== [src] ====================================

This is where all main calculations take place.

1) binary_evolution.py
This code solves the differential equations of motion of the binary for a specific set of initial conditions.
The program uses scipy.integrate.solve_ivp function to solve the differential equations, while the method of
numerical integration was chosen to be Dormand-Prince 853, although it can very easily be altered as a function
parameter. In order to check if the 2 ellipsoids came in contact during the simulation, an approximate collision
detection criterion was written.

3) plot_binary_evolution.py
This code just plots the results produced by code 1)

4) Sheeres_equilibrium_curve.py
This code calculates some theoretical equilibrium stuff.

5) masss_from_varying_primary_axes.py
This code does some data analysis stuff.

6) masss_from_varying_secondary_axes.py
This code does some data analysis stuff.

7) cost_function.cpp
This code does some data analysis stuff to verify some results from the previous
codes. In a nutshell, it runs code 1) many times in order to evaluate a cost function.
Due to the computational cost of the specific analysis, the code was written in C++.
In the beginning I wrote it in Python, but in order to get the same results I get from the
current C++ version, I had to wait for about 8 hours in my machine... The C++ version runs
in about 4 minutes in my machine... Since I had to solve the differential equations again for
different initial conditions, I used boost::odeint for the numerical integration and as far as
the calculation of the definite integral of the cost function, I write a Simpson integration algorithm.

8) plot_cost_function.py
This codes plots the data produced by the code 7)

=================================== [plots] ===================================
Here you will find all the plots generated by the codes in [src]

================================= [resources] =================================
All the data needed as inputs, as well as all the data produced from the codes in
[src] as outputs, are located in this directory. More specifically : "sys_inputs.txt" are the
necessary inputs to be read from the source codes. Each number of the "sys_inputs.txt"
file is explained in the "sys_inputs_explanation.txt". All the others .txt files in
this directory are produced from the codes in [src]. Basically they contain the
evolution of the mechanical charasteristics of the binary through time (e.g relative
distance, orientations on the plane of motion, energy and angular momentum), as well
as some "data analysis" stuff. Keep in mind that the numbers in "sys_inputs.txt" match
the 65803 Didymos binary asteroid, future target of the AIDA mission.

