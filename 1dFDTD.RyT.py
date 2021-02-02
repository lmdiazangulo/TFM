""" Cálculo de los coeficientes de reflexión y transmisión
    para un pulso Gaussiano que pasa a través de un material
    dieléctrico con una cierta conductividad.
"""
import numpy as np
from math import pi, sin, exp
from matplotlib import pyplot as plt

#Tamaño 
kmax = 200
ex = np.zeros(kmax+1)
hy = np.zeros(kmax+1)
ex2 = np.zeros(kmax+1)
hy2 = np.zeros(kmax+1)

# Tamaño de la celda
ddx = 0.001 

# Time step 
dt = ddx / 6e8 

boundary_low = [0,0,0,0]
boundary_high = [0,0,0,0] 

#Parametros del pulso
k_ini=20
k_1=100
k_2=130
t_0=40
s_0=12


# Create Dielectric Profile
epsz = 8.854e-12
epsilon = 4
sigma = 0.04
ca = np.ones(kmax+1)
cb = np.ones(kmax+1) * 0.5

cb_start = 110
cb_end=140
eaf = dt * sigma / (2 * epsz * epsilon)
ca[cb_start:cb_end] = (1 - eaf ) / (1 + eaf )
cb[cb_start:cb_end] = 0.5 / (epsilon * (1 + eaf ))

nsteps = 480


#Inicialización de algunos valores
ei=np.zeros(nsteps+1)
er=np.zeros(nsteps+1)
et=np.zeros(nsteps+1)

# FDTD Loop
for time_step in range(1, nsteps + 1):

    # Calculate the Ex field
    for k in range(1, kmax+1):
        ex[k] = ca[k] * ex[k] + cb[k] * (hy[k - 1] - hy[k])
        ex2[k]= ex2[k] + 0.5 * (hy2[k - 1] - hy2[k])
        if k==k_1:
            ei[time_step] = ex2[k]
            er[time_step] = ex[k]
        if k==k_2:
            et[time_step] = ex[k]    

    # Put a gaussian pulse at the low end (Con material y sin el)
    #Soft source
    pulse = exp(-0.5*((t_0-time_step)/s_0)**2)
    ex[k_ini] = 0.5* pulse + ex[k_ini]
    hy[k_ini-1] = 0.25* pulse + hy[k_ini-1]
    hy[k_ini] = 0.25*pulse + hy[k_ini]

    ex2[k_ini] = 0.5* pulse + ex2[k_ini]
    hy2[k_ini-1] = 0.25* pulse + hy2[k_ini-1]
    hy2[k_ini] = 0.25* pulse + hy2[k_ini]
    
    #Hard Source
    """pulse = exp(-0.5*((t_0-time_step)/s_0)**2)
    ex[k_ini] = (3.0/2.0)*pulse
    hy[k_ini-1] = (3.0/4.0)*pulse 
    hy[k_ini+1] =(3.0/4.0)* pulse 
    ex2[k_ini] = (3.0/2.0)*pulse 
    hy2[k_ini-1] = (3.0/4.0)*pulse 
    hy2[k_ini+1] = (3.0/4.0)*pulse"""

    # Absorbing Boundary Conditions
    ex[0] = boundary_low.pop(0)
    boundary_low.append(ex[1])
    ex[kmax] = boundary_high.pop(0)
    boundary_high.append(ex[kmax - 1])

    ex2[0] = boundary_low.pop(0)
    boundary_low.append(ex2[1])
    ex2[kmax] = boundary_high.pop(0)
    boundary_high.append(ex2[kmax - 1])

    # Calculate the Hy field
    for k in range(kmax):
        hy[k] = hy[k] + 0.5 * (ex[k] - ex[k + 1])
        hy2[k] = hy2[k] + 0.5 * (ex2[k] - ex2[k + 1])
    

er=ei-er

#Transformada de Fourier
eiw=np.fft.fft(ei)
erw=np.fft.fft(er)
etw=np.fft.fft(et)

#Se calcula el modulo de los numeros complejos
modeiw=np.linalg.norm(np.absolute(eiw))
moderw=np.linalg.norm(np.absolute(erw))
modetw=np.linalg.norm(np.absolute(etw))
 

R=moderw/modeiw
T=modetw/modeiw

print(moderw, modeiw)
print(R)
print(modetw, modeiw)
print(T)


# Plot the outputs in Fig. 1.6
plt.rcParams['font.size'] = 12
plt.figure(figsize=(8, 2.25))


plt.plot(ex, color='k', linewidth=1)
plt.ylabel('E$_x$', fontsize='14')
plt.xticks(np.arange(0, 199, step=20))
plt.xlim(0, 199)
plt.yticks(np.arange(-1, 1.2, step=1))
plt.ylim(-1.2, 1.2)
plt.text(50, 0.5, 'T = {}'.format(time_step),
         horizontalalignment='center')
plt.plot((0.5 / cb - 1) / 3, 'k--',
         linewidth=0.75) # The math on cb is just for scaling
plt.text(170, 0.5, 'Eps = {}'.format(epsilon),
         horizontalalignment='center')
plt.text(170, -0.5, 'Cond = {}'.format(sigma),
         horizontalalignment='center')
plt.xlabel('FDTD cells')

plt.subplots_adjust(bottom=0.25, hspace=0.45)
plt.show()