import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

N = 100
eP = np.zeros(N)
eN = np.zeros(N)
eP2 = np.zeros(N)
eN2 = np.zeros(N)

#define tight binding hamiltonian monolayer graphene
def Ep(kx,ky,a,t0,s0):
    f = np.exp(1J*kx*a/(np.sqrt(3)))+2*np.exp(-1J*kx*a/(2*np.sqrt(3)))*np.cos(ky*a/2)
    f_c = np.conjugate(f)
    sd = 1/(1-((s0**2)*(abs(f))**2))
    ss = s0*t0*(abs(f))**2
    a1 = np.array([[ss,-(t0)*f],[-(t0)*f_c,ss]])
    w, v = LA.eig(a1*sd)
    return np.real(w)
   
    #return z*t0*(abs(f))/(1-z*s0*(abs(f)))
a = 2.46
t = 3.033 
s0 = 0.129
k_x1 = 0
k_y1 = 0
k_x = (2*np.pi/(np.sqrt(3)*a))
k_y = (2*np.pi/(3*a))
for i in xrange(N):
    #calculate gamma-K path
    dx1 = (2*np.pi/(np.sqrt(3)*a))/N
    dy1 = (2*np.pi/(3*a))/N
    res = Ep(k_x1,k_y1,a,t,s0)
    eP2[i] = res[0]
    eN2[i] = res[1]
    k_x1 = k_x1+dx1
    k_y1 = k_y1+dy1
    
    #calculate K- M path
    dx = (2*np.pi/(np.sqrt(3)*a))
    dy = (2*np.pi/(3*a))/N
    res2 = Ep(k_x,k_y,a,t,s0)
    if res2[0]>res2[1]:
        eP[i] = res2[0]
        eN[i] = res2[1]
    else:
        eP[i] = res2[1]
        eN[i] = res2[0]
    k_x = dx
    k_y = k_y-dy
    
ET = np.concatenate((eP2,eP))
ET2 = np.concatenate((eN2,eN))

plt.plot(ET,'b')
plt.plot(ET2,'b')

plt.xticks( [0, 100, 200],['$\Gamma$', "K", 'M'])
plt.axvline(0, color='k', linestyle='dashed', linewidth=1)
plt.axvline(100, color='k', linestyle='dashed', linewidth=1)
plt.axvline(200, color='k', linestyle='dashed', linewidth=1)
plt.show()
