import matplotlib.pyplot as plt
import numpy as np
import csv

# Set folder directory

folder_dir = "./data"

# Parameters

with open("{}/Parameters.csv".format(folder_dir)) as F:
    csvwriter = csv.reader(F, delimiter=' ')
    for i in csvwriter:
        Parameters = list(map(lambda x: float(x), i[:]))

size = int(Parameters[0])

# Getting rho

with open("{}/rho0.csv".format(folder_dir)) as F:
    csvwriter = csv.reader(F, delimiter=' ')
    for i in csvwriter:
        rho = list(map(lambda x: float(x), i[:-1]))

rho0 = np.array(rho).reshape((size, size))

with open("{}/rho1.csv".format(folder_dir)) as F:
    csvwriter = csv.reader(F, delimiter=' ')
    for i in csvwriter:
        rho = list(map(lambda x: float(x), i[:-1]))

rho1 = np.array(rho).reshape((size, size))

# Quiver Graph

X = np.linspace(0, 1, size)
Y = np.linspace(0, 1, size)
U, V = np.meshgrid(X, Y)

with open("{}/mx.csv".format(folder_dir)) as F:
    csvwriter = csv.reader(F,delimiter=' ')
    for i in csvwriter:
        U = list(map(lambda x: float(x), i[:-1]))

with open("{}/my.csv".format(folder_dir)) as F:
    csvwriter = csv.reader(F,delimiter=' ')
    for i in csvwriter:
        V = list(map(lambda x: float(x), i[:-1]))

U = np.array(U).reshape((size,size))
V = np.array(V).reshape((size,size))

fig, ax = plt.subplots(1,3,figsize=(10,3))

Z = rho0

cp = ax[0].contourf(X, Y, rho0, cmap='jet')
ax[0].set_axis_off()
ax[0].set_title("$\\rho_0$")

Z = rho1

cp = ax[1].contourf(X, Y, rho1, cmap='jet')
ax[1].set_axis_off()
ax[1].set_title("$\\rho_1$")

skip = 3

q = ax[2].quiver(X[0:size:skip], Y[0:size:skip], U[0:size:skip,0:size:skip], V[0:size:skip,0:size:skip])
ax[2].set_axis_off()
ax[2].set_title("Tropical L1")

plt.savefig("./images/size-{}.png".format(size), dpi=500, bbox_inches = 'tight',pad_inches = 0)

plt.show()
