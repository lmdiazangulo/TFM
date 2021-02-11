import numpy as np
import scipy.constants as sp
from math import pi, sin, exp



class Mesh:
    def __init__(self, ncells, ddx, epsilon_r, sigma, start_m, end_m):    
        self.ncells=ncells
        self.ddx=ddx
        self.epsilon_r= epsilon_r
        self.sigma= sigma
        self.start_m= start_m
        self.end_m= end_m
    
    def material(self, k):
        self.k=k
        
        dt=self.ddx / sp.c

        ca = np.ones(self.ncells+1)
        cb = np.ones(self.ncells+1) * 0.5
        
        eaf = dt * self.sigma / (2 * sp.epsilon_0 * self.epsilon_r)
        ca[self.start_m : self.end_m] = (1 - eaf ) / (1 + eaf )
        cb[self.start_m : self.end_m] = 0.5 / (self.epsilon_r * (1 + eaf ))

        return ca[k], cb[k]   

    def boundarymur(self, ex, time_step): 
        self.ex=ex   
        self.time_step=time_step 

        if time_step >= 2:
                ex[time_step][0]=ex[time_step-2][1]
                ex[time_step][self.ncells]=ex[time_step-2][self.ncells-1]

        return ex        


class Source:
    def __init__(self, sourcetype, t_0, s_0):
        self.sourcetype=sourcetype
        self.t_0=t_0
        self.s_0=s_0

    def pulse(self, time_step):
        
        self.time_step=time_step
        
        if self.sourcetype == 'gauss':
            pulse = exp(-0.5*( (self.t_0 - time_step) / self.s_0 )**2)
        
        return pulse
