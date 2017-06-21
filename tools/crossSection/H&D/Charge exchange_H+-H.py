# Ref:
# !!!The unit of H+ energy is eV/amu
import numpy as np
import math
import matplotlib.pyplot as plt


filename_old = "original/Charge_exchange_H+-H.txt"
filename_new = "data/Charge_exchange_H+-H.dat"
data = np.loadtxt(filename_old)
np.savetxt(filename_new, data, fmt='%1.5e')
plt.plot(data[:,0], data[:,1], label = "Ionization_H_to_H+1")

'''
filename_old = "original/Ionization_D_to_D+1.dat"
data = np.loadtxt(filename_old)
plt.plot(data[:,0], data[:,1], label = "Ionization_D_to_D+1")
'''

plt.legend()
plt.xlim((0.0, 100.0))
#plt.savefig("fig/Ionization_H&D.png")
plt.savefig("fig/Charge_exchange_H+-H.png")
