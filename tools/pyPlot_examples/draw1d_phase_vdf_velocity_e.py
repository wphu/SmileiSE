from template import *

f0 = 1.0e-7

##read data from file
f=h5.File("restore/Restore0_global.h5")


v_number = 150
# for electron
m_ov_2T = 9.109382616e-31 / (2.0 * 24.8 * 1.602e-19)

# for D+ ion
#m_ov_2T = 2.0 * 1.67262158e-27 / (2.0 * 60 * 1.602e-19)

##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

##============ velocity distribution in x direction ================================
ax0 = fig.add_subplot(2,1,1)

val = f["/e/momentum0"]
val = val[...]
particle_number = val.shape[0]
v = np.zeros(particle_number)



#======= vx from simulation
val = f["/e/momentum0"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v0[iv]

vmin = v.min() * 1.0
vmax = -vmin
v_step = (vmax - vmin) / v_number
v_step_half = 0.5 * v_step

v_x = np.zeros(v_number)
v_y = np.zeros(v_number)
v_y_maxwell = np.zeros(v_number)

v_x[0] = vmin
for i in np.arange(1, v_number):
	v_x[i] = v_x[i-1] + v_step

for v_i in v:
	i = int( (v_i -vmin + v_step_half) / v_step )
	if i < 0:
		v_y[0] = v_y[0] + 1
	elif i >= v_number:
		v_y[v_number-1] = v_y[v_number-1] + 1
	else:
		v_y[i] = v_y[i] + 1
v_y = 1.0 * v_y / (particle_number * v_step)
v_y = v_y / f0

line0=ax0.plot(v_x, v_y, label = r'PIC-$v_x$', linestyle = linestyles[0])


#======= vy from simulation
val = f["/e/momentum1"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v0[iv]



vmin = v.min() * 1.0
vmax = -vmin
v_step = (vmax - vmin) / v_number
v_step_half = 0.5 * v_step

v_x = np.zeros(v_number)
v_y = np.zeros(v_number)
v_y_maxwell = np.zeros(v_number)

v_x[0] = vmin
for i in np.arange(1, v_number):
	v_x[i] = v_x[i-1] + v_step

for v_i in v:
	i = int( (v_i -vmin + v_step_half) / v_step )
	if i < 0:
		v_y[0] = v_y[0] + 1
	elif i >= v_number:
		v_y[v_number-1] = v_y[v_number-1] + 1
	else:
		v_y[i] = v_y[i] + 1
v_y = 1.0 * v_y / (particle_number * v_step)
v_y = v_y / f0

line0=ax0.plot(v_x, v_y, label = r'PIC-$v_y$', linestyle = linestyles[1])


#======= vz from simulation
val = f["/e/momentum2"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v0[iv]



vmin = v.min() * 1.0
vmax = -vmin
v_step = (vmax - vmin) / v_number
v_step_half = 0.5 * v_step

v_x = np.zeros(v_number)
v_y = np.zeros(v_number)
v_y_maxwell = np.zeros(v_number)

v_x[0] = vmin
for i in np.arange(1, v_number):
	v_x[i] = v_x[i-1] + v_step

for v_i in v:
	i = int( (v_i -vmin + v_step_half) / v_step )
	if i < 0:
		v_y[0] = v_y[0] + 1
	elif i >= v_number:
		v_y[v_number-1] = v_y[v_number-1] + 1
	else:
		v_y[i] = v_y[i] + 1
v_y = 1.0 * v_y / (particle_number * v_step)
v_y = v_y / f0

line0=ax0.plot(v_x, v_y, label = r'PIC-$v_z$', linestyle = linestyles[2])



#======= from theory
for i in np.arange(0, v_number):
	v_y_maxwell[i] = math.sqrt(m_ov_2T/3.14) * math.exp( -m_ov_2T * v_x[i] * v_x[i] )
v_y_maxwell = v_y_maxwell / f0

line0=ax0.plot(v_x, v_y_maxwell, label = r'Theory', linestyle = linestyles[3])

ax0.set_xlim((vmin, vmax))
ax0.set_ylim((0, 2.3e-7))

ax0.grid()
ax0.legend(framealpha=1)

ax0.xaxis.set_major_formatter(yformatter)
ax0.yaxis.set_major_formatter(yformatter)

major_ticks = np.arange(0, 2.3, 1.0)
ax0.set_yticks(major_ticks)

ax0.set_xlabel(r'Velocity (m/s)')
ax0.set_ylabel(r'$F_M\ \mathrm{(10^{-7})}$')

ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)


v_number = 150
##==================================== Speed distribution =====================================
ax0 = fig.add_subplot(2,1,2)

val = f["/e/momentum0"]
val = val[...]
particle_number = val.shape[0]
v = np.zeros(particle_number)

val = f["/e/momentum0"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v[iv] + v0[iv] * v0[iv]

val = f["/e/momentum1"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v[iv] + v0[iv] * v0[iv]

val = f["/e/momentum2"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v[iv] + v0[iv] * v0[iv]


for iv in np.arange(0, particle_number):
	v[iv] = math.sqrt( v[iv] )



vmin = v.min() * 1.0
vmax = v.max() * 1.0
v_step = (vmax - vmin) / v_number
v_step_half = 0.5 * v_step

v_x = np.zeros(v_number)
v_y = np.zeros(v_number)
v_y_maxwell = np.zeros(v_number)

v_x[0] = vmin
for i in np.arange(1, v_number):
	v_x[i] = v_x[i-1] + v_step

#======= from theory
for i in np.arange(0, v_number):
	v_y_maxwell[i] = math.pow(m_ov_2T/3.14, 1.5) * math.exp( -m_ov_2T * v_x[i] * v_x[i] ) * 4.0 * 3.14 * v_x[i] * v_x[i]

#======= from simulation
for v_i in v:
	i = int( (v_i -vmin + v_step_half) / v_step )
	if i < 0:
		v_y[0] = v_y[0] + 1
	elif i >= v_number:
		v_y[v_number-1] = v_y[v_number-1] + 1
	else:
		v_y[i] = v_y[i] + 1
v_y = 1.0 * v_y / (particle_number * v_step)

v_y = v_y / f0
v_y_maxwell = v_y_maxwell / f0

line0=ax0.plot(v_x, v_y, label = r'PIC', linestyle = linestyles[0])
line0=ax0.plot(v_x, v_y_maxwell, label = r'Theory', linestyle = linestyles[1])

ax0.set_xlim((vmin, vmax))
ax0.set_ylim((0, 3.2e-7))

ax0.grid()
ax0.legend(framealpha=1)

ax0.xaxis.set_major_formatter(yformatter)
ax0.yaxis.set_major_formatter(yformatter)

major_ticks = np.arange(0, 3.2, 1.0)
ax0.set_yticks(major_ticks)

ax0.set_xlabel(r'Speed (m/s)')
ax0.set_ylabel(r'$f_M\ \mathrm{(10^{-7})}$')

ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)


#plt.legend()
fig.savefig("vdf_vx_e.pdf", dpi = 300)
#plt.axis('equal')
#plt.show()         #The command is OK
