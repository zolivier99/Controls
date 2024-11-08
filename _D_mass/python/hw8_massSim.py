import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from massDynamics import massDynamics
from ctrlPD import ctrlPD
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter

mass = massDynamics()
controller = ctrlPD()
reference = signalGenerator(amplitutde=1.0, frequency=0.04)
disturbance = signalGenerator(amplitude=0.25)

dataPlot = dataPlotter()
animation = massAnimation()
t = P.t_start
y = mass.h()

while t < P.t_end:
    
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = reference.square(t)
        d = disturbance.step(t)
        n = 0.0
        x = mass.state
        u = controller.updtae(r, x)
        y = mass.update(u + d)
        t = t + P.Ts
        
    animation.update(mass.state)
    dataPlot.update(t, r, mass.state, u)
    plt.pause(0.01)
    
print('Press key to close')
plt.waitforbuttonpress()
plt.close()