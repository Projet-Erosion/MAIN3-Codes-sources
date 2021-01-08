import numpy as np
from scipy.sparse import spdiags
from scipy.sparse import linalg as li
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# list abscices
nspacestep = 500  # nombre pas d'espace
ntimestep  = 500 # nombre de pas de temps
nsavestep  = 20  # time step qu'on dessine
time = 100000 # temps total de la simulation en annee
x = np.linspace(-10000.,10000.,nspacestep) # longeur du  domaine en m
dx = x[2] - x[1] # pas de d'espace
dt = time/(ntimestep-1) # pas de temps
periode = 20.e3; # periode de variation du niveau marin en annee
D = 0.5
# pente initiale
teta = 1.
# array ordonnees
z   = x*np.tan((teta*np.pi)/180) #topo initiale m
plt.plot(x,z)
plt.pause(2)
tec = 2.e-3; # vitesse de soulevement en m par an
h_star = 5. #m limite d'action des vagues/3 ou 4 m
C = 0.007 # m/an taux d'erosion cotiere maximum
K0 = 1. # m^2/an taux d'erosion par diffusion
A_sl = 40. # amplitude des variations de niveau marin.

z2D=z
tplot=0.
t=0.

for k in range(0,ntimestep):

    SL = A_sl*np.sin(t/(periode)*np.pi)
    h      = SL-z
    dz_ero = -dt*C*np.exp(-h/h_star)*(h>0)

    
    K = K0*(np.linspace(1, 1, x.size)*(h<0) + np.linspace(10, 10, x.size)*(h>=0))

    Alpha = list(K*dt/dx**2)


    diags = [[], [], []]
    diags[0] = [(Alpha[j] - (Alpha[j+1]/4 - Alpha[j-1]/4))*(D-1) for j in range(1, len(x)-1)]
    diags[0].append(0)
    diags[0].append(0)
    diags[1] = [2*Alpha[j]*(1-D) + 1 for j in range(0, len(x))]
    diags[2] = [(Alpha[j] + (Alpha[j+1]/4 - Alpha[j-1]/4))*(D-1) for j in range(1, len(x)-1)]
    diags[2].insert(0, 0)
    diags[2].insert(0, 0)

    sparse = spdiags(diags, [-1, 0, 1], len(x), len(x)).tolil()
    sparse[0,0] = 1
    sparse[len(x) - 1,len(x) - 1] = 1
    sparse[0,1] = 0


    z_temp = list(z)
    z_old = [Alpha[i+1]*D*z_temp[i+1] + (-2*Alpha[i+1]*D + 1)*z_temp[i] + Alpha[i+1]*D*z_temp[i-1] for i in range(1, len(x)-1)]
    z_old.insert(0, z_temp[0])
    z_old.append(z_temp[-1])
    z_old = np.array(z_old)

    z_old = z_old + dt*tec + dz_ero
    z = li.spsolve(sparse, z_old)
    if k%(nsavestep)==0:
       tplot=np.append(tplot,t)
       plt.plot(x,z)
       plt.pause(0.1)
    t=t+dt

plt.show()
fig = plt.figure()
ax = plt.axes(projection='3d')
for k in range(0,ntimestep):
	if k%(nsavestep)==0:
		z2D=np.vstack((z2D,z))
X , T = np.meshgrid(x, tplot)
ax.plot_surface(X,T,z2D,cmap='terrain', edgecolor='none')
ax.set_title('topographie au cours du temps')
plt.show()

