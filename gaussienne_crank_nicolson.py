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


diags2 = [[], [], []]
diags2[0] = [Alpha*D for j in range(0, len(x))]
diags2[1] = [-2*Alpha*D + 1 for j in range(0, len(x))]
diags2[2] = diags2[0]
sparse2 = scipy.sparse.spdiags(diags2, [-1, 0, 1], len(x), len(x)).tolil()
sparse2[0,1] = 0
sparse2[0,0] = 1
sparse2[len(x) - 1,len(x) - 2] = 0
sparse2[len(x) - 1,len(x) - 1] = 1

t = dt
while t <= 200000 :
	print(z)
	z_old = np.dot(sparse2.todense(), np.array(z))
	z_old = np.array(z_old)
	print(z_old)
	z = list(li.spsolve(sparse, z_old))
	t = t + dt
	plt.plot(x,z)
	plt.pause(0.001)
plt.show()
