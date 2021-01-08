import math
import matplotlib.pyplot as plt

# list abscices
x = [i for i in range(-10**3, 10**3+1, 50)]
dx = x[2] - x[1]

# pente initiale
teta = 1

# list ordonnees
z = [i*math.tan((teta*math.pi)/180) for i in x]
Xc = 0
dt = 10
tec = 10**-3
h_star = 1
C = 0.1
plt.plot(x,z)
plt.pause(2)
for t in range(0, 2010, dt) :
    SL=math.sin(t)
    h = []
    i = 0
    while z[i] <= SL :
        h.append(SL - z[i])
        i = i + 1
    k = 0
    for i in range(0, len(h)) :
        z[i] = z[i]+dt*(tec-C*math.exp(-h[i]/h_star))
    for i in range(len(h), len(x)) :
        z[i] = z[i] + dt*tec
    plt.plot(x,z)
    if t <= 40 :
        plt.pause(1)
    else :
        plt.pause(0.0001)
plt.show()
