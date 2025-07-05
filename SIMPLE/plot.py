import numpy as np
import matplotlib.pyplot as plt

p = np.loadtxt("p.txt", delimiter=',')
u = np.loadtxt("u.txt", delimiter=',')
v = np.loadtxt("v.txt", delimiter=',')

xfile = np.loadtxt("x_domain.txt", delimiter=",")
yfile = np.loadtxt("y_domain.txt", delimiter=",")

X = xfile.reshape(-1)
Y = yfile.reshape(-1)
LX = X[-1]
LY = Y[-1]

XX, YY = np.meshgrid(X, Y)
XG, YG = np.meshgrid(X, LY-Y)

plt.figure()
plt.contourf(XG, YG, p)
plt.colorbar()
plt.streamplot(XG[::-1, :], YG[::-1, :], u[::-1, :], v[::-1, :], color="black")
plt.quiver(XG, YG, u, v, scale=50)
#plt.title("Lid Driven Cavity")
plt.xlabel("x domain")
plt.ylabel("y domain")
plt.gca().set_aspect('equal')
plt.show()