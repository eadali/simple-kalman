#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  2 13:49:29 2022

@author: eadali
"""

from simple_kalman import Kalman
from advanced_pid.models import MassSpringDamper
from matplotlib import pyplot as plt
from numpy import diff

m, k, b = 1.0, 1.0, 0.2
system = MassSpringDamper(mass=m, spring_const=k, damping_const=b)
system.set_initial_value(initial_position=1.0, initial_velocity=0.0)

A = [[ 0.0,  1.0],
     [-k/m, -b/m]]
B = [[0.0],
     [1/m]]
C = [[1, 0]]
D = 0
Q = [[0.00001, 0],
     [0, 0.00001]]
R = 0.02

kalman = Kalman(A, B, C, D, Q, R)
# Control loop
time, meas, cont = [], [], []
for i in range(800):
    # Get current measurement from system
    timestamp, measurement = system.get_measurement()
    
    # Calculate control signal by using PID controller
    control = kalman(timestamp, 0, measurement)
    
    # Feed control signal to system
    system.set_input(0)
    
    # Record for plotting
    time.append(timestamp)
    meas.append(measurement)
    cont.append(control)

# Plot result
fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
fig.suptitle('Mass-Spring-Damper system')
ax1.set_ylabel('Measured Position [m]')
ax1.plot(time, meas, 'b')
ax1.grid()
ax2.set_ylabel('Force [N]')
ax2.plot(time, cont, 'g')
ax2.grid()
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('Derivative Term')
ax3.plot(time[1:], diff(meas)/diff(time), 'r')
ax3.grid()
plt.show()