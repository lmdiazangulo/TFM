from mesh import Mesh, Source
from solver import FDTD, FFT
from viewer import Animator




malla1=Mesh(200,0.001,4,0.04,110,140)
malla2=Mesh(200,0.001,1,0,110,140)
pulso=Source('gauss',40,12)

ex1= FDTD(malla1,pulso,3e-9).FDTDLoop(20)
ex2= FDTD(malla2,pulso,3e-9).FDTDLoop(20)



r= FFT(ex1,ex2,100,130).RyT()[0]
t= FFT(ex1,ex2,100,130).RyT()[1]
freq=FFT(ex1,ex2,100,130).frequency(3e-9)



Animator().animationex(ex1)

Animator().fftgraph(freq,r,t)