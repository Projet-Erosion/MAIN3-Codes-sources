import numpy as np
import scipy.sparse
from scipy.sparse import linalg as li
import matplotlib.pyplot as plt

# list abscices
x = [i for i in range(-200, 201, 10)]
dx = x[2] - x[1]

def gaussienne(sigma, mu, x) :
    z = [(100000/(sigma*np.sqrt(2*np.pi)))*np.exp(-((i - mu)**2)/(2*(sigma**2))) for i in x]
    return z

# list ordonnées gaussienne
z = gaussienne(100, 5, x)
dt = 2000
K = 1

Alpha = K*dt/dx**2
D = 0.5

plt.plot(x,z)
plt.pause(1)

diags = [[], [], []]
diags[0] = [Alpha*(D-1) for j in range(0, len(x))]
diags[1] = [2*Alpha*(1-D) + 1 for j in range(0, len(x))]
diags[2] = diags[0]
sparse = scipy.sparse.spdiags(diags, [-1, 0, 1], len(x), len(x)).tolil()
sparse[0,1] = 0
sparse[0,0] = 1
sparse[len(x) - 1,len(x) - 2] = 0
sparse[len(x) - 1,len(x) - 1] = 1

t = dt
while t <= 50000 :
    z_old = [Alpha*D*z[i+1] + (-2*Alpha*D + 1)*z[i] + Alpha*D*z[i-1] for i in range(1, len(x)-1)]
    z_old.insert(0, z[0])
    z_old.append(z[-1])
    z_old = np.array(z_old)

    z = list(li.spsolve(sparse, z_old))
    t = t + dt
    plt.plot(x,z)
    plt.pause(0.001)

#plt.title("Application du processus de diffusion à une gaussienne")
plt.xlabel("x (en m)")
plt.ylabel("f (en m)")
plt.legend()
plt.show()
