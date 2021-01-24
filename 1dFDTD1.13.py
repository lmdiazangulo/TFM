import math
import numpy as np
from math import exp
import matplotlib
from matplotlib import pyplot as plt


#Inicializacion de las cantidades

#Numero de celdas
kmax=200

#Inicializacion a cero de los campos
#Los vectores no incluyen el ultimo punto
#De esta forma va de 0 a 200 incluido
ex=np.zeros(kmax+1)
hy=np.zeros(kmax+1)

#Parametros del pulso
kc=int(kmax/2)
t_0=40
s_0=12

#Numero de iteraciones o pasos de tiempo
nsteps=500
#Ciclo para el tiempo en el metodo FDTD
for time_step in range(1,nsteps+1):

    #Calculo de ex
    #El rango debe empezar en 1 ya que aparece k-1 en hy
    for k in range(1, kmax+1):
        ex[k]=ex[k]+0.5*(hy[k-1]-hy[k])

    #Poner un pulso Gaussiano en el medio
    pulse=exp(-0.5*((t_0-time_step)/s_0)**2)
    hy[kc]=pulse
    hy[kc-1]=-hy[kc]
    #Calculo de hy
    #El rango no puede incluir kmax porque hay un k+1
    for k in range(kmax):
        hy[k]=hy[k]+0.5*(ex[k]-ex[k+1])


# Plot the outputs 
plt.rcParams['font.size'] = 12
plt.figure(figsize=(8, 3.5))

plt.subplot(211)
plt.plot(ex, color='k', linewidth=1)
plt.ylabel('E$_x$', fontsize='14')
plt.xticks(np.arange(0, 201, step=20))
plt.xlim(0, 200)
plt.yticks(np.arange(-1, 1.2, step=1))
plt.ylim(-1.2, 1.2)
plt.text(100, 0.5, 'T = {}'.format(time_step),
horizontalalignment='center')

plt.subplot(212)
plt.plot(hy, color='k', linewidth=1)
plt.ylabel('H$_y$', fontsize='14')
plt.xlabel('FDTD cells')
plt.xticks(np.arange(0, 201, step=20))
plt.xlim(0, 200)
plt.yticks(np.arange(-1, 1.2, step=1))
plt.ylim(-1.2, 1.2)
plt.subplots_adjust(bottom=0.2, hspace=0.45)
plt.show()