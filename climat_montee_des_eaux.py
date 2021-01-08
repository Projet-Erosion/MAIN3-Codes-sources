import numpy as np
from scipy.sparse import spdiags
from scipy.sparse import linalg as li
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pandas as pd

data = pd.read_csv("donnees.csv",sep=';') #lecture du fichier .csv
data_t = data.iloc[:,0]
data_SL = interpolate.interp1d(np.linspace(0, 492000, 3937), data.iloc[:,1])

# list abscices
nspacestep = 500  # nombre pas d'espace
ntimestep  = 500 # nombre de pas de temps
nsavestep  = 20  # time step qu'on dessine
time = 100 # temps total de la simulation en annee
x = np.linspace(-200.,200.,nspacestep) # longeur du  domaine en m

dx = x[2] - x[1] # pas de d'espace
dt = time/ntimestep # pas de temps

periode = 20.e3; # periode de variation du niveau marin en annee
D = 0.5
# pente initiale
teta = 1.
# array ordonnees
z   = x*np.tan((teta*np.pi)/180) #topo initiale m

tec = 2.e-3; # vitesse de soulevement en m par an
h_star = 5. #m limite d'action des vagues/3 ou 4 m
C = 0.007 # m/an taux d'erosion cotiere maximum
K0 = 1. # m^2/an taux d'erosion par diffusion

tplot=0.
t=0.

#_____________________________________________________________

# Nous sommes maintenant en 2000 et le niveau de la mer est à posé à 0

# D'après les prédictions, le niveau de la mer va évoluer exponetiellement d'ici 2100
# Nous allons étudier 2 cas :

# - Le niveau de la mer monte jusqu'à 2m d'ici 2100
# - Le niveau de la mer monte jusqu'à 1m d'ici 2100

plt.plot(x,z)

k = 0
print_SL_x = []
while list(x*(z <= 0))[k] != 0 :
    print_SL_x.append(list(x*(z <= 0))[k])
    k = k + 1

print_SL_x = np.array(print_SL_x)
print_SL_y = np.zeros(print_SL_x.size)

plt.plot(print_SL_x , print_SL_y, color="blue")

for k in range(0,ntimestep+1): # teamps géologiques et formation des térasses

    SL = 2*np.exp((t - 100)/100) # sea level pour la mer à 2m en 2100
    #SL = np.exp((t - 100)/100) # sea level pour la mer à 1m en 2100
    h = SL-z
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

    tplot=np.append(tplot,t)

    t=t+dt

plt.plot(x,z)

k = 0
print_SL_x = []
while list(x*(h >= 0))[k] != 0 :
    print_SL_x.append(list(x*(h >= 0))[k])
    k = k + 1

print_SL_x = np.array(print_SL_x)
print_SL_y = np.ones(print_SL_x.size)*SL

plt.plot(print_SL_x , print_SL_y, color="red")

plt.show()

z2D=np.ones(x.size)*z[0]
fig = plt.figure()
ax = plt.axes(projection='3d')
for k in range(0,ntimestep+1):
	z2D=np.vstack((z2D,z))
X , T = np.meshgrid(x, tplot)
ax . set_xlabel ( 'x' ) 
ax . set_ylabel ( 'y' ) 
ax . set_zlabel ( 'z' );
ax.plot_surface(X,T,z2D,cmap='terrain', edgecolor='none')
ax.set_title('topographie')
plt.show()