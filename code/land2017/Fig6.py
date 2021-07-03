"""
This script approximately reproduces Figure 6 in Land et al (2017).
The Ca transient is taken from Rice et al (2008), but reparameterized
to be visually similar to the one used in Land et al (2017).
"""

import land2017 as land_model
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0,1000,101)

#Get indices representing active force and Cai, for plotting
force_ind = land_model.monitor_indices("Ta")
cai_ind = land_model.monitor_indices("cai")

#solve and plot for different values of lambda
for lmbda in np.linspace(0.9,1.1,5):
    #Ca parameters chosen to mimic left panel of Fig 6:
    p = land_model.init_parameter_values(lmbda=lmbda, Ca_amplitude=0.6, Ca_diastolic=0.17, tau1=75, tau2=175)
    #initialize and solve the model
    init = land_model.init_state_values()
    s = odeint(land_model.rhs,init,t,(p,))
    force = []
    cai = []

    #loop through the solution and compute outputs of interest:
    for tn,sn in zip(t,s):
        m = land_model.monitor(sn,tn,p)
        force.append(m[force_ind])
        cai.append(m[cai_ind])

    plt.plot(t,force,label=f'$\lambda$ = {lmbda:.2f}')

plt.title('Fig 6 right, isometric force')
plt.ylabel('Force (kPa)')
plt.xlabel('Time (ms)')
plt.legend()

plt.figure()
plt.plot(t,cai)
plt.title('Fig 6 left, Ca transient')
plt.ylabel('Intracellular calcium concentration ()$\mu$M)')
plt.xlabel('Time (ms)')

plt.show()
