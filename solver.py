import numpy as np
from math import pi, sin, exp
import scipy.constants as sp

class FDTD:
    def __init__(self, mesh, pulse, time):
        self.mesh=mesh
        self.pulse=pulse
        self.time=time
      
    def FDTDLoop(self):
        
        dt=self.mesh.ddx / (2*sp.c)
        nsteps= int(self.time  / dt)

        # COMENTAR: Mejor quitar nsteps, no guardar siempre todo...
        ex=np.zeros((nsteps+1,self.mesh.ncells+1))
        hy=np.zeros((nsteps+1,self.mesh.ncells+1))
    
       
        for time_step in range(1, nsteps + 1):

            # Calculate the Ex field, for cycle using slice notation
            # COMENTAR: self.mesh.material() se evalúa en cada paso.
            ex[time_step][1:] = self.mesh.material()[0][1:] * ex[time_step-1][1:] + \
            self.mesh.material()[1][1:] * (hy[time_step-1][:-1] - hy[time_step-1][1:])

            
            ex[time_step][self.pulse.k_ini] +=  0.5*self.pulse.pulse(time_step) 
            
            #Condiciones de contorno
            self.mesh.boundarymur(ex,time_step)
            
            # Calculate the Hy field, slicing notation
            hy[time_step][:-1] = hy[time_step-1][:-1] + \
            0.5 * (ex[time_step][:-1] - ex[time_step][1:])   

            t= time_step+1/2
            hy[time_step][self.pulse.k_ini] += 0.25* self.pulse.pulse(t) 
            hy[time_step][self.pulse.k_ini-1] += 0.25* self.pulse.pulse(t)                

        return ex



#Clase para la Trasformada Rápida de Fourier
# COMENTAR: Esto es mas un namespace que una clase. 
# COMENTAR: Cuanto menos estado, mejor
class FFT:
    def __init__(self, e1t, e2t, k_1, k_2):
        self.e1t=e1t
        self.e2t=e2t
        self.k_1=k_1
        self.k_2=k_2


    def fftfunction(self):

        #Hay que cancelar la parte incidente
        self.e1t[:,self.k_1] += (-1.0)* self.e2t[:,self.k_1]
        
        e1wk_1=np.fft.fft(self.e1t[:,self.k_1])
        e2wk_1=np.fft.fft(self.e2t[:,self.k_1])

        e1wk_2=np.fft.fft(self.e1t[:,self.k_2])
        e2wk_2=np.fft.fft(self.e2t[:,self.k_2])
    

        return e1wk_1, e2wk_1, e1wk_2, e2wk_2
    
    def RyT(self):

        #R=np.abs(self.fftfunction()[1]-self.fftfunction()[0]) / np.abs(self.fftfunction()[1])
        R=np.abs(self.fftfunction()[0]) / np.abs(self.fftfunction()[1])
        T=np.abs(self.fftfunction()[2]) / np.abs(self.fftfunction()[3])
        
        return R, T

    def frequency(self,time):
        self.time=time

        N=len(self.e1t)

        freq= (1.0/time) * np.arange(N)         

        return freq