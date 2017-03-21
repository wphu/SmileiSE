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




# 1 Excitation_C+2_2s21S-2s2p3P =====================
Vif = 6.5
A   = -4.024e-2
B   = 2.914
C   = -4.344
D   = 2.121
E   = 0.0
X1   = 2.0
P   = 5.434e-1
Q   = 3.453e-1

excitation_id = "data/Excitation_C+2_2s21S-2s2p3P"
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



# 2 Excitation_C+2_2s21S-2s2p1P =====================
Vif = 12.7
A   = 6.261e-1
B   = 3.044
C   = 5.384e-1
D   = 0.0
E   = 4.364


excitation_id = "data/Excitation_C+2_2s21S-2s2p1P"
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



# 3 Excitation_C+2_2s21S-2p23P =====================
Vif = 17.2
A   = -5.732e-3
B   = 8.467e-2
C   = -1.676e-1
D   = 1.307e-1
E   = 0.0
X1   = 2.0
P   = 5.389e-3
Q   = 7.167e-3

excitation_id = "data/Excitation_C+2_2s21S-2p23P"
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



# 4 Excitation_C+2_2s21S-2p21D =====================
Vif = 18.2
A   = 2.915e-1
B   = 1.655e-1
C   = -3.445e-1
D   = 3.720e-1
E   = 0.0


excitation_id = "data/Excitation_C+2_2s21S-2p21D"
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



# 5 Excitation_C+2_2s21S-2p21S =====================
Vif = 22.9
A   = 2.092e-2
B   = 2.289e-1
C   = -4.268e-1
D   = 2.137e-1
E   = 0.0


excitation_id = "data/Excitation_C+2_2s21S-2p21S"
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



# 6 Excitation_C+2_2s21S-2s3s3S =====================
Vif = 29.5
A   = 3.910e-3
B   = 2.740e-1
C   = 1.120
D   = 1.110
E   = 0.0


excitation_id = "data/Excitation_C+2_2s21S-2s3s3S"
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



# 7 Excitation_C+2_2s21S-2s3s1S =====================
Vif = 30.6
A   = 4.260e-1
B   = 4.350e-1
C   = -2.280
D   = 1.760
E   = 0.0


excitation_id = "data/Excitation_C+2_2s21S-2s3s1S"
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




# 8 Excitation_C+2_2s21S-2s3p1P =====================
Vif = 32.1
A   = -1.140e-1
B   = -5.230e-1
C   = 8.110e-1
D   = 0.0
E   = 3.730e-1


excitation_id = "data/Excitation_C+2_2s21S-2s3p1P"
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




# 9 Excitation_C+2_2s21S-2s3p3P =====================
Vif = 32.2
A   = -2.190e-2
B   = 3.660e-1
C   = -1.030
D   = 8.240e-1
E   = 0.0

excitation_id = "data/Excitation_C+2_2s21S-2s3p3P"
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




# 10 Excitation_C+2_2s21S-2s3d3D =====================
Vif = 33.5
A   = 4.190e-2
B   = 6.490e-2
C   = 1.100e-2
D   = 1.050e-1
E   = 0.0


excitation_id = "data/Excitation_C+2_2s21S-2s3d3D"
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



# 11 Excitation_C+2_2s21S-2s3d1D =====================
Vif = 34.3
A   = 7.330e-1
B   = 5.400e-2
C   = -2.110
D   = 1.490
E   = 0.0


excitation_id = "data/Excitation_C+2_2s21S-2s3d1D"
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




# 12 Excitation_C+2_2s2p3P-2p23P =====================
Vif = 10.6
A   = 9.842
B   = -8.660
C   = 2.184e1
D   = 0.0
E   = 1.631


excitation_id = "data/Excitation_C+2_2s2p3P-2p23P"
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



# 13 Excitation_C+2_2s2p3P-2p21D =====================
Vif = 11.7
A   = -1.822e-1
B   = 4.490
C   = -4.889
D   = 2.056
E   = 0.0


excitation_id = "data/Excitation_C+2_2s2p3P-2p21D"
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



# 14 Excitation_C+2_2s2p3P-2p21S =====================
Vif = 16.5
A   = 3.280e-4
B   = 1.676e-1
C   = 2.778e-1
D   = -4.524e-1
E   = 0.0


excitation_id = "data/Excitation_C+2_2s2p3P-2p21S"
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



# 15 Excitation_C+2_2s2p3P-2s3s1S =====================
Vif = 24.1
A   = -1.200e-2
B   = 1.310
C   = -4.030
D   = 3.590
E   = 0.0


excitation_id = "data/Excitation_C+2_2s2p3P-2s3s1S"
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



# 16 Excitation_C+2_2s2p3P-2s3p1P =====================
Vif = 25.6
A   = -1.310e-2
B   = 1.580
C   = -4.800
D   = 4.270
E   = 0.0


excitation_id = "data/Excitation_C+2_2s2p3P-2s3p1P"
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



# 17 Excitation_C+2_2s2p3P-2s3p3P =====================
Vif = 25.7
A   = 1.980
B   = 1.490e1
C   = -4.330e1
D   = 3.170e1
E   = 0.0


excitation_id = "data/Excitation_C+2_2s2p3P-2s3p3P"
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



# 18 Excitation_C+2_2s2p3P-2s3d3D  =====================
Vif = 27.0
A   = -2.900e-1
B   = 2.090
C   = 7.110
D   = 0.0
E   = 1.010e1


excitation_id = "data/Excitation_C+2_2s2p3P-2s3d3D"
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




# 19 Excitation_C+2_2s2p3P-2s3s1D =====================
Vif = 27.8
A   = 9.840e-3
B   = 5.080e-1
C   = -7.040e-1
D   = 9.550e-1
E   = 0.0


excitation_id = "data/Excitation_C+2_2s2p3P-2s3s1D"
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



# 20 Excitation_C+2_2s2p1P-2p23P =====================
Vif = 4.3
A   = -8.624e-3
B   = 4.141
C   = -7.089
D   = 3.824
E   = 0.0
X1  = 1.85
P   = 7.072e-1
Q   = 2.899e-1


excitation_id = "data/Excitation_C+2_2s2p1P-2p23P"
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



# 21 Excitation_C+2_2s2p1P-2p21D =====================
Vif = 4.3
A   = 3.762
B   = 9.351
C   = -3.004
D   = 0.0
E   = 7.320


excitation_id = "data/Excitation_C+2_2s2p1P-2p21D"
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




# 22 Excitation_C+2_2s2p1P-2p21S =====================
Vif = 10.1
A   = 3.032
B   = -2.357
C   = 2.675
D   = 0.0
E   = 2.991


excitation_id = "data/Excitation_C+2_2s2p1P-2p21S"
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




# 23 Excitation_C+2_2p23P-2p21D =====================
Vif = 1.01
A   = 1.925
B   = 1.440e1
C   = -3.560e1
D   = 2.727e1
E   = 0.0


excitation_id = "data/Excitation_C+2_2p23P-2p21D"
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



# 24 Excitation_C+2_2p23P-2p21S =====================
Vif = 5.8
A   = 4.319e-2
B   = 1.004
C   = -1.037
D   = 3.630e-1
E   = 0.0


excitation_id = "data/Excitation_C+2_2p23P-2p21S"
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



# 25 Excitation_C+2_2p21D-2p21S =====================
Vif = 4.79
A   = 8.948e-1
B   = -8.607e-1
C   = 1.007
D   = -4.136e-1
E   = 0.0


excitation_id = "data/Excitation_C+2_2p21D-2p21S"
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
plt.savefig("fig/Excitation_C+2.png")


# ======== Plot ===========================================
#plt.plot(data_np[:,0], data_np[:,1])
#plt.savefig("fig.png")
#plt.show()
