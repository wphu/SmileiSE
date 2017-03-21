# Ref: Cross section database for carbon atoms and ions: Electron impact ionization
#      excitation, and charge exchange in collisions with hydrogen atoms
# Journal: Atomic data and nuclear data tables 92 (2006) 407-455
import numpy as np
import math
import matplotlib.pyplot as plt

# ============= Ionization: equation (17) ==========================

I = 0.0
A = np.zeros(5)

def ionization( E ):
    s = 0.0
    for i in range(2, 6):
        s = s + A[i-1] * math.pow( (1.0 - I / E), i-1 )
    return 1.0e-17 * ( A[0] * math.log( E / I ) + s ) / (I * E)



# C => C+1 =====================
I       = 10.6
A[0]    = 1.829
A[1]    = -1.975
A[2]    = 1.149
A[3]    = -3.583
A[4]    = 2.451

ionization_id = "data/Ionization_C_to_C+1"
data = [[I],[0.0]]
for energy in np.arange(12.0, 100.0, 1.0):
    CS = ionization(energy)
    data[0].append(energy)
    data[1].append(CS)

for energy in np.arange(100.0, 1000.0, 10.0):
    CS = ionization(energy)
    data[0].append(energy)
    data[1].append(CS)

data_np = np.asarray(data)
data_np = np.transpose(data_np)
np.savetxt(ionization_id+".dat", data_np, fmt='%1.4e')

plt.plot(data_np[:,0], data_np[:,1], label = ionization_id)
#plt.savefig(ionization_file+".png")


# C+1 => C+2 =====================
I       = 24.4
A[0]    = 8.390e-1
A[1]    = -7.950e-1
A[2]    = 3.263
A[3]    = -5.382
A[4]    = 3.476

ionization_id = "data/Ionization_C+1_to_C+2"
data = [[I],[0.0]]
for energy in np.arange(25.0, 100.0, 1.0):
    CS = ionization(energy)
    data[0].append(energy)
    data[1].append(CS)

for energy in np.arange(100.0, 1000.0, 10.0):
    CS = ionization(energy)
    data[0].append(energy)
    data[1].append(CS)

data_np = np.asarray(data)
data_np = np.transpose(data_np)
np.savetxt(ionization_id+".dat", data_np, fmt='%1.4e')

plt.plot(data_np[:,0], data_np[:,1], label = ionization_id)
#plt.savefig(ionization_file+".png")





# C+2 => C+3 =====================
I       = 41.4
A[0]    = 4.009e-1
A[1]    = -3.518e-1
A[2]    = 2.375
A[3]    = -3.992
A[4]    = 2.794

ionization_id = "data/Ionization_C+2_to_C+3"
data = [[I],[0.0]]
for energy in np.arange(42.0, 100.0, 1.0):
    CS = ionization(energy)
    data[0].append(energy)
    data[1].append(CS)

for energy in np.arange(100.0, 1000.0, 10.0):
    CS = ionization(energy)
    data[0].append(energy)
    data[1].append(CS)

data_np = np.asarray(data)
data_np = np.transpose(data_np)
np.savetxt(ionization_id+".dat", data_np, fmt='%1.4e')

plt.plot(data_np[:,0], data_np[:,1], label = ionization_id)

plt.legend()
plt.savefig("fig/Ionization_C.png")


# ======== Plot ===========================================
#plt.plot(data_np[:,0], data_np[:,1])
#plt.savefig("fig.png")
#plt.show()
