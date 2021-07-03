import rice_model_2008 as rice
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


t = np.linspace(0,1000,1001)
sl_init = 2.2
KSE_vals = [50,10,5,3,2,1.4,1]

force_ind = rice.monitor_indices("active")
sl_ind = rice.state_indices('SL')
for kse in KSE_vals:
    init = rice.init_state_values(SL=sl_init)
    p = rice.init_parameter_values(SLset=2.2,KSE= kse)
    s = odeint(rice.rhs,init,t,(p,))
    force = []
    plt.figure(1)
    plt.plot(t,s[:,sl_ind])
   
    for tn,sn in zip(t,s):
        m = rice.monitor(sn,tn,p)
        force.append(m[force_ind])
    plt.figure(2)
    plt.plot(t,force)
    
plt.show()
