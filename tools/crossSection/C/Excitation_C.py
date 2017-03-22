# Ref: Cross section database for carbon atoms and ions: Electron impact ionization
#      excitation, and charge exchange in collisions with hydrogen atoms
# Journal: Atomic data and nuclear data tables 92 (2006) 407-455
import numpy as np
import math
import matplotlib.pyplot as plt

# ============= Excitation: equation (6/7) ==========================

id = 0
maxCS = 2.0e-20

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




# 1 Excitation_C_2s22p23P-2s22p21D =====================
id = id + 1
Vif = 1.26
A   = 1.277
B   = 1.304e1
C   = -3.314e1
D   = 6.357e1
E   = -4.621e1
F   = 4.0e-2

excitation_id = "data/Excitation_C_2s22p23P-2s22p21D"
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





# 2 Excitation_C_2s22p23P-2s22p21S =====================
id = id + 1
Vif = 2.68
A   = 1.598e-1
B   = 2.836
C   = -9.813
D   = -1.012e1
E   = 1.233e-1


excitation_id = "data/Excitation_C_2s22p23P-2s22p21S"
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



# 3 Excitation_C_2s22p23P-2s2p35S =====================
id = id + 1
Vif = 4.18
A   = -2.302
B   = 8.923
C   = -2.873e1
D   = 4.555e1
E   = -2.431e1
F   = 1.323e-1

excitation_id = "data/Excitation_C_2s22p23P-2s2p35S"
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




# 4 Excitation_C_2s22p23P-2s22p3s3P =====================
id = id + 1
Vif = 7.48
A   = -6.117
B   = 2.989
C   = 4.148
D   = 0.0
E   = 1.222e1

excitation_id = "data/Excitation_C_2s22p23P-2s22p3s3P"
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




# 5 Excitation_C_2s22p23P-2s22p3s1P =====================
id = id + 1
Vif = 7.68
A   = 2.423
B   = 7.590e-1
C   = -6.486
D   = 3.143
E   = -2.074e1
F   = 6.744e-1

excitation_id = "data/Excitation_C_2s22p23P-2s22p3s1P"
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


# 6 Excitation_C_2s22p23P-2s2p33D =====================
id = id + 1
Vif = 7.94
A   = 1.307e1
B   = -3.052e1
C   = 1.835e1
D   = 0.0
E   = 1.001e1


excitation_id = "data/Excitation_C_2s22p23P-2s2p33D"
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








# 7 Excitation_C_2s22p23P-2s22p3p1P =====================
id = id + 1
Vif = 8.53
A   = 1.909
B   = -6.101e-1
C   = -2.336
D   = -3.093e2
E   = 7.450e2
F   = 1.491

excitation_id = "data/Excitation_C_2s22p23P-2s22p3p1P"
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






# 8 Excitation_C_2s22p23P-2s22p3p3D =====================
id = id + 1
Vif = 8.64
A   = 9.753e-1
B   = 1.908
C   = 2.161
D   = -4.891
E   = 0.0

excitation_id = "data/Excitation_C_2s22p23P-2s22p3p3D"
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





# 9 Excitation_C_2s22p23P-2s22p3p3S =====================
id = id + 1
Vif = 8.77
A   = 1.886e-3
B   = 2.263e-1
C   = 8.836e-1
D   = -1.031
E   = 0.0


excitation_id = "data/Excitation_C_2s22p23P-2s22p3p3S"
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





# 10 Excitation_C_2s22p23P-2s22p3p3P =====================
id = id + 1
Vif = 8.85
A   = 7.728
B   = 8.661e-1
C   = -1.640e1
D   = 7.754
E   = 0.0


excitation_id = "data/Excitation_C_2s22p23P-2s22p3p3P"
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






# 11 Excitation_C_2s22p23P-2s22p3p1D =====================
id = id + 1
Vif = 9.00
A   = -1.389
B   = 1.062
C   = -3.742e-1
D   = 4.786
E   = -4.126
F   = 2.606e-1

excitation_id = "data/Excitation_C_2s22p23P-2s22p3p1D"
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




# 12 Excitation_C_2s22p23P-2s22p3p1S =====================
id = id + 1
Vif = 9.20
A   = 3.567e-2
B   = 5.081e-1
C   = -7.062e-1
D   = 1.463
E   = -2.569
F   = 4.101e-1

excitation_id = "data/Excitation_C_2s22p23P-2s22p3p1S"
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




# 13 Excitation_C_2s22p23P-2s2p33P =====================
id = id + 1
Vif = 9.33
A   = 1.146e1
B   = -2.477e1
C   = 1.281e1
D   = 0.0
E   = 6.328


excitation_id = "data/Excitation_C_2s22p23P-2s2p33P"
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




# 14 Excitation_C_2s22p23P-2s22p3d1D =====================
id = id + 1
Vif = 9.63
A   = 5.747e-2
B   = -1.444e-1
C   = 1.793
D   = -1.719
E   = 0.0


excitation_id = "data/Excitation_C_2s22p23P-2s22p3d1D"
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




# 15 Excitation_C_2s22p23P-2s22p3d3F =====================
id = id + 1
Vif = 9.69
A   = -3.067
B   = 2.352
C   = -3.127
D   = -9.123
E   = 1.656e1
F   = 1.476e-1

excitation_id = "data/Excitation_C_2s22p23P-2s22p3d3F"
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



# 16 Excitation_C_2s22p23P-2s22p3d3D =====================
id = id + 1
Vif = 9.71
A   = 4.081
B   = -9.274
C   = 5.423
D   = 0.0
E   = 4.970


excitation_id = "data/Excitation_C_2s22p23P-2s22p3d3D"
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



# 17 Excitation_C_2s22p23P-2s22p3d1F =====================
id = id + 1
Vif = 9.73
A   = 1.183
B   = -8.798e-1
C   = 7.839
D   = -1.718e1
E   = -8.149e-1
F   = 5.66e-1

excitation_id = "data/Excitation_C_2s22p23P-2s22p3d1F"
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



# 18 Excitation_C_2s22p23P-2s22p3d3P =====================
id = id + 1
Vif = 9.83
A   = 1.247
B   = -5.247
C   = 3.983
D   = 0.0
E   = 3.876


excitation_id = "data/Excitation_C_2s22p23P-2s22p3d3P"
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


# 19 Excitation_C_2s22p23P-2s2p31D =====================
id = id + 1
Vif = 12.2
A   = -3.820
B   = 2.082
C   = 2.801
D   = -9.506
E   = 1.552e1
F   = 2.937e-1

excitation_id = "data/Excitation_C_2s22p23P-2s2p31D"
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



# 20 Excitation_C_2s22p23P-2s2p33S =====================
id = id + 1
Vif = 13.1
A   = 8.801e-1
B   = -9.630
C   = 8.179
D   = 0.0
E   = 1.199e1


excitation_id = "data/Excitation_C_2s22p23P-2s2p33S"
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




# 21 Excitation_C_2s22p23P-2s2p31P =====================
id = id + 1
Vif = 14.9
A   = 2.975
B   = -4.065
C   = -2.943e-1
D   = 7.361e1
E   = -3.359e2
F   = 1.40

excitation_id = "data/Excitation_C_2s22p23P-2s2p31P"
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



# 22 Excitation_C_2s22p21D-2s22p21S =====================
id = id + 1
Vif = 1.42
A   = 2.029
B   = -1.093e1
C   = 2.049e1
D   = -1.221e1
E   = 0.0


excitation_id = "data/Excitation_C_2s22p21D-2s22p21S"
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


# 23 Excitation_C_2s22p21D-2s2p33D =====================
id = id + 1
Vif = 6.68
A   = 1.003e2
B   = -5.439e1
C   = -1.260e2
D   = 1.861e2
E   = -9.357e2
F   = 6.9e-1

excitation_id = "data/Excitation_C_2s22p21D-2s2p33D"
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



# 24 Excitation_C_2s22p21D-2s22p3s3P =====================
id = id + 1
Vif = 6.22
A   = 9.576
B   = -6.855e-1
C   = 2.121
D   = 3.983
E   = -7.383e1
F   = 4.750e-1

excitation_id = "data/Excitation_C_2s22p21D-2s22p3s3P"
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


# 25 Excitation_C_2s22p21D-2s22p3s1P =====================
id = id + 1
Vif = 6.42
A   = -1.795
B   = -1.134
C   = 2.995
D   = 0.0
E   = 4.794


excitation_id = "data/Excitation_C_2s22p21D-2s22p3s1P"
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
plt.savefig("fig/Excitation_C.png")


# ======== Plot ===========================================
#plt.plot(data_np[:,0], data_np[:,1])
#plt.savefig("fig.png")
#plt.show()
