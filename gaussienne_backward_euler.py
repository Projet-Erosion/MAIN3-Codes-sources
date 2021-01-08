import math
import matplotlib.pyplot as plt

# list abscices
x = [i for i in range(-200, 201, 2)]
dx = x[2] - x[1]

def gaussienne(sigma, mu, x) :
	z = [(100000/(sigma*math.sqrt(2*math.pi)))*math.exp(-((i - mu)**2)/(2*(sigma**2))) for i in x]
	return z

# list ordonnées gaussienne
z = gaussienne(100, 5, x)
dt = 1e-2
K = 1

# K*dt/dx² <= 1/2 : Stable
Alpha = K*dt/dx**2

plt.plot(x,z)
plt.pause(1)
t = dt
while t <= 200 :
	z_old = z
	for i in range(1, len(x)-1) :
		z[i] = Alpha*(z[i+1] - 2*z[i] + z[i-1]) + z[i]
	t = t + dt
	plt.plot(x,z)
	plt.pause(0.001)
plt.show()