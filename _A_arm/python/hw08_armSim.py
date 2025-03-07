import matplotlib.pyplot as plt
import numpy as np
import armParam as P
from signalGenerator import signalGenerator
from armAnimation import armAnimation
from dataPlotter import dataPlotter
from armDynamics import armDynamics
from ctrlPD import ctrlPD

# instantiate arm, controller, and reference classes
arm = armDynamics()
controller = ctrlPD()
reference = signalGenerator(amplitude=50*np.pi/180.0, 
                            frequency=0.05)
disturbance = signalGenerator(amplitude=0.0)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = armAnimation()

t = P.t_start  # time starts at t_start
y = arm.h()  # output of system at start of simulation

while t < P.t_end:  # main simulation loop
    # Get referenced inputs from signal generators
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

    # updates control and dynamics at faster simulation rate
    while t < t_next_plot: 
        r = reference.square(t)
        x = arm.state
        u = controller.update(r, x)  # update controller
        y = arm.update(u)  # propagate system
        t += P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(arm.state)
    dataPlot.update(t, arm.state, u, r)

    # the pause causes the figure to be displayed for simulation
    plt.pause(0.0001)  

# Keeps the program from closing until the user presses a button.
dataPlot.write_data_file()
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
