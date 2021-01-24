import numpy as np
from math import pi, sin, exp
from matplotlib import pyplot as plt


class FDTD:
    def __init__(self, ncells, nsteps, ddx, t_0, s_0, epsilon_r, sigma, start_m, end_m):
        self.ncells=ncells
        self.nsteps=nsteps
        self.ddx=ddx
        self.t_0=t_0
        self.s_0=s_0
        self.epsilon_r= epsilon_r
        self.sigma= sigma
        self.start_m= start_m
        self.end_m= end_m

    def pulsogauss(self, time_step):
        
        self.time_step=time_step

        gausspulse = np.zeros(self.nsteps + 1)
        
        gausspulse[time_step] = exp(-0.5*( (self.t_0 - time_step) / self.s_0 )**2)
        
        return gausspulse[time_step]


    def material(self, k):
        self.k=k
        
        dt=self.ddx / 6e8

        epsz = 8.854e-12
        ca = np.ones(self.ncells+1)
        cb = np.ones(self.ncells+1) * 0.5
        
        eaf = dt * self.sigma / (2 * epsz * self.epsilon_r)
        ca[self.start_m : self.end_m] = (1 - eaf ) / (1 + eaf )
        cb[self.start_m : self.end_m] = 0.5 / (self.epsilon_r * (1 + eaf ))

        return ca[k], cb[k]
    
    def boundary(self, ex, b_l, b_h):
        self.ex=ex
        self.b_l=b_l
        self.b_h=b_h

        ex[0] = b_l.pop(0)
        b_l.append(ex[1])
        ex[self.ncells] = b_h.pop(0)
        b_h.append(ex[self.ncells -1])

        return ex, b_l, b_h


    def FDTDLoop(self, k_ini):
        self.k_ini=k_ini

        ex=np.zeros(self.ncells+1)
        hy=np.zeros(self.ncells+1)

        b_l = [0] * int(self.epsilon_r)
        b_h = [0] * int(self.epsilon_r)

        for time_step in range(1, self.nsteps + 1):

            # Calculate the Ex field
            for k in range(1, self.ncells + 1):
                ex[k] = self.material(k)[0] * ex[k] + self.material(k)[1] * (hy[k - 1] - hy[k])

            ex[self.k_ini] = 0.5 * self.pulsogauss(time_step) + ex[self.k_ini]
            hy[self.k_ini-1] = 0.25 * self.pulsogauss(time_step) + hy[self.k_ini-1]
            hy[self.k_ini+1] = 0.25 * self.pulsogauss(time_step) + hy[self.k_ini+1]
            
            #Sin el print statement no funciona correctamente
            """self.boundary(ex, b_l, b_h)
        
            print(self.boundary(ex, b_l, b_h))"""

            # Calculate the Hy field
            for k in range(self.ncells):
                hy[k] = hy[k] + 0.5 * (ex[k] - ex[k + 1])    

        return ex


class FastFourierTrasnform:
    def __init__(self, ex):
        self.ex=ex

    def fft(self):

        exw=np.fft.fft(self.ex)
        modexw=np.linalg.norm(np.absolute(exw))

        return modexw
    


print(FastFourierTrasnform(FDTD(200,480,0.001,40,12,4,0.04,110,140).FDTDLoop(20)).fft()) 





"""for i in range(1, 100):
    print(FDTD(200, 300).ex(i)+i)        

for time_step in range(1, 300):
    print(FDTD(200,300).pulsogauss(40,12,time_step)) """   



# Plot the outputs in Fig. 1.6
plt.rcParams['font.size'] = 12
plt.figure(figsize=(8, 2.25))


plt.plot(FDTD(200,482,0.001,40,12,4,0.04,110,140).FDTDLoop(20), color='k', linewidth=1)
plt.ylabel('E$_x$', fontsize='14')
plt.xticks(np.arange(0, 199, step=20))
plt.xlim(0, 199)
plt.yticks(np.arange(-1, 1.2, step=1))
plt.ylim(-1.2, 1.2)


plt.xlabel('FDTD cells')

plt.subplots_adjust(bottom=0.25, hspace=0.45)
plt.show()