# Ref: 1973 Electron modecule collision ionization in hydrogen and deuterium
import numpy as np
import math
import matplotlib.pyplot as plt


filename_old = "original/Ionization_H_to_H+1.dat"
filename_new = "data/Ionization_H_to_H+1.dat"
data = np.loadtxt(filename_old)
np.savetxt(filename_new, data, fmt='%1.5e')
plt.plot(data[:,0], data[:,1], label = "Ionization_H_to_H+1")

plt.legend()
plt.xlim((0.0, 500.0))
plt.savefig("fig/Ionization_H.png")
