import matplotlib.pyplot as plt
import blockBeamParam as P
from blockBeamDynamics import blockBeamDynamics
from ctrlPID import ctrlPID
from signalGenerator import signalGenerator
from blockBeamAnimation import blockBeamAnimation
from dataPlotter import dataPlotter

# instantiate blockBeam, controller, and reference classes
blockBeam = blockBeamDynamics(alpha=0.2)
controller = ctrlPID()
reference = signalGenerator(amplitude=0.125, frequency=0.02, y_offset=0.25)
disturbance = signalGenerator(amplitude=0.25, frequency=0.0)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockBeamAnimation()
t = P.t_start  # time starts at t_start
y = blockBeam.h()  # output of system at start of simulation

while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

    while t < t_next_plot:  # updates control and dynamics at faster simulation rate
        r = reference.square(t)  # reference input
        d = 0.0 #disturbance.step(t)  # input disturbance
        n = 0.0  #noise.random(t)  # simulate sensor noise, will use in future assignments
        u = controller.update(r, y + n)  # update controller
        y = blockBeam.update(u + d)  # propagate system, "d" is a disturbance used in future assignments
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(blockBeam.state)
    dataPlot.update(t, r, blockBeam.state, u)
    plt.pause(0.0001)  # the pause causes the figure to be displayed during the simulation

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
