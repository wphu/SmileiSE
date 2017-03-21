# Ref: Cross section database for carbon atoms and ions: Electron impact ionization
#      excitation, and charge exchange in collisions with hydrogen atoms
# Journal: Atomic data and nuclear data tables 92 (2006) 407-455
import numpy as np
import math
import matplotlib.pyplot as plt

# ============= Excitation: equation (6/7) ==========================

Vif = 0.0
A = 0.0
B = 0.0
C = 0.0
D = 0.0
E = 0.0
F = 0.0

# Equation (6)
def excitation6( Te ):
    X = Te / Vif
    Qif =  A + B / X + C / (X*X) + D / (X*X*X) + E * math.log(X)
    return 1.1969e-19 * Qif / ( 1.0 * Te )
# Equation (7)
def excitation7( Te ):
    X = Te / Vif
    Qif =  A / (X*X) + B * math.exp( - F * X ) + C * math.exp(- 2.0 * F * X) + D * math.exp( - 3.0 * F * X ) + E * math.exp( - 4.0 * F * X )
    return 1.1969e-19 * Qif / ( 1.0 * Te )




# 1 Excitation_C+3_2s2S-2p2P =====================
Vif = 8.0
A   = 4.744
B   = 4.440
C   = 5.971e-1
D   = 0.0
E   = 4.519


excitation_id = "data/Excitation_C+3_2s2S-2p2P"
data = [[Vif],[0.0]]
for energy in np.arange(float(int(Vif+1)), 100.0, 1.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

for energy in np.arange(100.0, 1000.0, 10.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

data_np = np.asarray(data)
data_np = np.transpose(data_np)
np.savetxt(excitation_id+".dat", data_np, fmt='%1.4e')

plt.plot(data_np[:,0], data_np[:,1], label = excitation_id)


# 2 Excitation_C+3_2s2S-3s2S =====================
Vif = 37.6
A   = 4.846e-1
B   = -2.743e-1
C   = 4.768e-1
D   = -3.971e-1
E   = 0.0
X1  = 1.072
P   = 7.390e-1
Q   = 2.460e-1

excitation_id = "data/Excitation_C+3_2s2S-3s2S"
data = [[Vif],[0.0]]
for energy in np.arange(float(int(Vif+1)), 100.0, 1.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

for energy in np.arange(100.0, 1000.0, 10.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

data_np = np.asarray(data)
data_np = np.transpose(data_np)
np.savetxt(excitation_id+".dat", data_np, fmt='%1.4e')

plt.plot(data_np[:,0], data_np[:,1], label = excitation_id)



# 3 Excitation_C+3_2s2S-3p2P =====================
Vif = 40.0
A   = -5.731e-1
B   = 9.548e-1
C   = -8.931e-2
D   = 0.0
E   = 5.638e-1
X1  = 1.019
P   = -3.412
Q   = 9.141

excitation_id = "data/Excitation_C+3_2s2S-3p2P"
data = [[Vif],[0.0]]
for energy in np.arange(float(int(Vif+1)), 100.0, 1.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

for energy in np.arange(100.0, 1000.0, 10.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

data_np = np.asarray(data)
data_np = np.transpose(data_np)
np.savetxt(excitation_id+".dat", data_np, fmt='%1.4e')

plt.plot(data_np[:,0], data_np[:,1], label = excitation_id)



# 4 Excitation_C+3_2s2S-3d2D =====================
Vif = 40.0
A   = 1.252
B   = -1.243
C   = 4.93e-1
D   = 1.326e-1
E   = 0.0


excitation_id = "data/Excitation_C+3_2s2S-3d2D"
data = [[Vif],[0.0]]
for energy in np.arange(float(int(Vif+1)), 100.0, 1.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

for energy in np.arange(100.0, 1000.0, 10.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

data_np = np.asarray(data)
data_np = np.transpose(data_np)
np.savetxt(excitation_id+".dat", data_np, fmt='%1.4e')

plt.plot(data_np[:,0], data_np[:,1], label = excitation_id)



# 5 Excitation_C+3_2s2S-4p2P =====================
Vif = 50.6
A   = 2.775e-1
B   = -6.663e-1
C   = 4.959e-1
D   = 0.0
E   = 3.530e-6


excitation_id = "data/Excitation_C+3_2s2S-4p2P"
data = [[Vif],[0.0]]
for energy in np.arange(float(int(Vif+1)), 100.0, 1.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

for energy in np.arange(100.0, 1000.0, 10.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

data_np = np.asarray(data)
data_np = np.transpose(data_np)
np.savetxt(excitation_id+".dat", data_np, fmt='%1.4e')

plt.plot(data_np[:,0], data_np[:,1], label = excitation_id)



# 6 Excitation_C+3_2p2P-3s2S =====================
Vif = 29.4
A   = 1.273
B   = -5.038
C   = 5.374
D   = 0.0
E   = 2.270e-3


excitation_id = "data/Excitation_C+3_2p2P-3s2S"
data = [[Vif],[0.0]]
for energy in np.arange(float(int(Vif+1)), 100.0, 1.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

for energy in np.arange(100.0, 1000.0, 10.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

data_np = np.asarray(data)
data_np = np.transpose(data_np)
np.savetxt(excitation_id+".dat", data_np, fmt='%1.4e')

plt.plot(data_np[:,0], data_np[:,1], label = excitation_id)


# 7 Excitation_C+3_2p2P-3d2D =====================
Vif = 32.1
A   = 1.899e1
B   = -4.193e1
C   = 2.857e1
D   = 0.0
E   = 4.027e-3


excitation_id = "data/Excitation_C+3_2p2P-3d2D"
data = [[Vif],[0.0]]
for energy in np.arange(float(int(Vif+1)), 100.0, 1.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

for energy in np.arange(100.0, 1000.0, 10.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

data_np = np.asarray(data)
data_np = np.transpose(data_np)
np.savetxt(excitation_id+".dat", data_np, fmt='%1.4e')

plt.plot(data_np[:,0], data_np[:,1], label = excitation_id)



# 8 Excitation_C+3_2p2P-4s2S =====================
Vif = 41.8
A   = 1.964e-1
B   = -6.743e-1
C   = 7.116e-1
D   = 0.0
E   = 2.411e-4


excitation_id = "data/Excitation_C+3_2p2P-4s2S"
data = [[Vif],[0.0]]
for energy in np.arange(float(int(Vif+1)), 100.0, 1.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

for energy in np.arange(100.0, 1000.0, 10.0):
    CS = excitation6(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

data_np = np.asarray(data)
data_np = np.transpose(data_np)
np.savetxt(excitation_id+".dat", data_np, fmt='%1.4e')

plt.plot(data_np[:,0], data_np[:,1], label = excitation_id)



plt.xlim((0.0, 100))
#plt.legend()
plt.savefig("fig/Excitation_C+3.png")


# ======== Plot ===========================================
#plt.plot(data_np[:,0], data_np[:,1])
#plt.savefig("fig.png")
#plt.show()
