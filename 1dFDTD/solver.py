import numpy as np
from math import pi, sin, exp
import scipy.constants as sp


class FDTD:
    def __init__(self, mesh, pulse, time):
        self.mesh=mesh
        self.pulse=pulse
        self.time=time


    def FDTDLoop(self, k_ini):
        self.k_ini=k_ini
        
        dt=self.mesh.ddx / sp.c
        nsteps= int(self.time  / dt)

        ex=np.zeros((nsteps+1,self.mesh.ncells+1))
        hy=np.zeros((nsteps+1,self.mesh.ncells+1))
    

        for time_step in range(1, nsteps + 1):

            # Calculate the Ex field
            for k in range(1, self.mesh.ncells + 1):
                ex[time_step][k] = self.mesh.material(k)[0] * ex[time_step-1][k] + \
                self.mesh.material(k)[1] * (hy[time_step-1][k - 1] - hy[time_step-1][k])


            ex[time_step][self.k_ini] += 0.5 * self.pulse.pulse(time_step) 
            hy[time_step-1][self.k_ini-1] += 0.25 * self.pulse.pulse(time_step) 
            hy[time_step-1][self.k_ini] += 0.25 * self.pulse.pulse(time_step) 
            
            #Condiciones de contorno
            self.mesh.boundarymur(ex,time_step)
            
            # Calculate the Hy field
            for k in range(self.mesh.ncells):
                hy[time_step][k] = hy[time_step-1][k] + \
                0.5 * (ex[time_step][k] - ex[time_step][k + 1])    

        return ex



#Clase para la Trasformada Rápida de Fourier

class FFT:
    def __init__(self, e1t, e2t, k_1, k_2):
        self.e1t=e1t
        self.e2t=e2t
        self.k_1=k_1
        self.k_2=k_2


    def fft(self):
        
        #Hay que cancelar la parte incidente
        self.e2t[:,self.k_1] = self.e1t[:,self.k_1]- self.e2t[:,self.k_1]

        e1w=np.fft.fft(self.e1t[:,self.k_1])
        e1w2=np.fft.fft(self.e1t[:,self.k_2])

        e2w=np.fft.fft(self.e2t[:,self.k_1])
        e2w2=np.fft.fft(self.e2t[:,self.k_2])
    

        return e1w, e1w2, e2w, e2w2
    
    def RyT(self):

        R=np.abs(self.fft()[2]) / np.abs(self.fft()[0])
        T=np.abs(self.fft()[3]) / np.abs(self.fft()[1])

        return R, T     

    def frequency(self,time):
        self.time=time

        N=len(self.e1t)

        freq= (1.0/time) * np.arange(N)         

        return freq