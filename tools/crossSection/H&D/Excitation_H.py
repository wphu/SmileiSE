# Ref: from NIST website
import numpy as np
import math
import matplotlib.pyplot as plt



for id in np.arange(2,11):
    filename_old = "original/Excitation_H-1s-" + str(id) + "p.txt"
    filename_new = "data/Excitation_H_1s-" + str(id) + "p.dat"
    data = np.loadtxt(filename_old)
    data[:,1] = data[:,1] * 1.0e-20
    #print data
    np.savetxt(filename_new, data, fmt='%1.5e')
    plt.plot(data[:,0], data[:,1], label = "Excitation_H_1s-" + str(id) + "p")

plt.legend()
plt.xlim((0.0, 1000.0))
#plt.ylim((0.0, 1.0e-21))
plt.savefig("fig/Excitation_H.png")
