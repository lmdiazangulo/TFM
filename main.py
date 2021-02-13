from mesh import Mesh, Source
from solver import FDTD, FFT
from viewer import Animator




malla1=Mesh(200,0.001,4,0.04,110,140)
malla2=Mesh(200,0.001,1,0,110,140)
pulso=Source('gauss',40,12,20)

#cambiar k_ini al pulso
ex1= FDTD(malla1,pulso,5e-9).FDTDLoop()
ex2= FDTD(malla2,pulso,5e-9).FDTDLoop()


Animator().animationex(ex1,malla1)

r= FFT(ex1,ex2,40,160).RyT()[0]
t= FFT(ex1,ex2,40,160).RyT()[1]
freq=FFT(ex1,ex2,40,160).frequency(5e-9)




Animator().fftgraph(freq,r,t)


