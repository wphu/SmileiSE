# Ref: Cross section database for carbon atoms and ions: Electron impact ionization
#      excitation, and charge exchange in collisions with hydrogen atoms
# Journal: Atomic data and nuclear data tables 92 (2006) 407-455
import numpy as np
import math
import matplotlib.pyplot as plt

# ============= Excitation: equation (6/7) ==========================
id = 0
maxCS = 0.0e-20

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



# 1 Excitation_C+1_2s22p2P-2s2p24P =====================
id = id + 1
Vif = 4.74
A   = 4.360e-1
B   = 1.627
C   = -2.923
D   = 1.094e1
E   = -6.580
F   = 5.769e-2

excitation_id = "data/Excitation_C+1_2s22p2P-2s2p24P"
data = [[Vif],[0.0]]
for energy in np.arange(float(int(Vif+1)), 100.0, 1.0):
    CS = excitation7(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

for energy in np.arange(100.0, 1000.0, 10.0):
    CS = excitation7(energy)
    if CS > 0.0 :
        data[0].append(energy)
        data[1].append(CS)

data_np = np.asarray(data)
data_np = np.transpose(data_np)
if data_np[:,1].max() >= maxCS:
    np.savetxt(excitation_id+".dat", data_np, fmt='%1.5e')
    plt.plot(data_np[:,0], data_np[:,1], label = str(id))





# 2 Excitation_C+1_2s22p2P-2s2p22D =====================
id = id + 1
Vif = 9.27
A   = -2.155
B   = 1.374e1
C   = -6.562
D   = 0.0
E   = 1.018


excitation_id = "data/Excitation_C+1_2s22p2P-2s2p22D"
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
if data_np[:,1].max() >= maxCS:
    np.savetxt(excitation_id+".dat", data_np, fmt='%1.5e')
    plt.plot(data_np[:,0], data_np[:,1], label = str(id))


# 3 Excitation_C+1_2s22p2P-2s2p22S =====================
id = id + 1
Vif = 12.0
A   = 8.572e-1
B   = 2.079e-1
C   = 7.058e-1
D   = 0.0
E   = 4.136


excitation_id = "data/Excitation_C+1_2s22p2P-2s2p22S"
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
if data_np[:,1].max() >= maxCS:
    np.savetxt(excitation_id+".dat", data_np, fmt='%1.5e')
    plt.plot(data_np[:,0], data_np[:,1], label = str(id))



# 4 Excitation_C+1_2s22p2P-2s2p22P =====================
id = id + 1
Vif = 13.7
A   = -6.883e-1
B   = 4.106
C   = 3.435e-1
D   = 0.0
E   = 1.631e1


excitation_id = "data/Excitation_C+1_2s22p2P-2s2p22P"
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
if data_np[:,1].max() >= maxCS:
    np.savetxt(excitation_id+".dat", data_np, fmt='%1.5e')
    plt.plot(data_np[:,0], data_np[:,1], label = str(id))


# 5 Excitation_C+1_2s22p2P-2s23s2S =====================
id = id + 1
Vif = 14.4
A   = -9.843e-1
B   = 2.983
C   = -9.975e-1
D   = 0.0
E   = 8.910


excitation_id = "data/Excitation_C+1_2s22p2P-2s23s2S"
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
if data_np[:,1].max() >= maxCS:
    np.savetxt(excitation_id+".dat", data_np, fmt='%1.5e')
    plt.plot(data_np[:,0], data_np[:,1], label = str(id))


# 6 Excitation_C+1_2s22p2P-2s23p2P =====================
id = id + 1
Vif = 16.3
A   = 1.678
B   = -2.238
C   = 4.829
D   = 0.0
E   = 0.0


excitation_id = "data/Excitation_C+1_2s22p2P-2s23p2P"
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
if data_np[:,1].max() >= maxCS:
    np.savetxt(excitation_id+".dat", data_np, fmt='%1.5e')
    plt.plot(data_np[:,0], data_np[:,1], label = str(id))



# 7 Excitation_C+1_2s22p2P-2s23d2D =====================
id = id + 1
Vif = 18.0
A   = -1.703
B   = 2.675
C   = 1.165
D   = 0.0
E   = 5.791


excitation_id = "data/Excitation_C+1_2s22p2P-2s23d2D"
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
if data_np[:,1].max() >= maxCS:
    np.savetxt(excitation_id+".dat", data_np, fmt='%1.5e')
    plt.plot(data_np[:,0], data_np[:,1], label = str(id))




plt.xlim((0.0, 100))
plt.legend()
plt.savefig("fig/Excitation_C+1.png")


# ======== Plot ===========================================
#plt.plot(data_np[:,0], data_np[:,1])
#plt.savefig("fig.png")
#plt.show()
