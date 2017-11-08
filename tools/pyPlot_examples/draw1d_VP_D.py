from template import *

label_fontsize = 26
legend_fontsize = 23

t = 19


mass_e = 9.109382616e-31
mass_i = 2.0 * 1.67262158e-27
const_pi = 3.1415926
const_e = 1.602e-19


##read data from file
f=h5.File("ref/data_global.h5")
print f.keys()

group = f['/Fields']
dims = group.attrs['dims_global']
dims = dims[...]

print dims

nx = dims[3]


dx = 0.5e-2  # unit (mm)
x_axis = np.zeros(5)
x_sheath = np.zeros(5, dtype = np.int)
Te0 = np.zeros(5)
Ti0 = np.zeros(5)
P0 = np.zeros(5)
V0 = np.zeros(5)

x_axis[0] = 0.0
x_axis[1] = 0.2
x_axis[2] = 0.4
x_axis[3] = 0.6
x_axis[4] = 0.8

x_sheath[0] = 28
x_sheath[1] = 25
x_sheath[2] = 20
x_sheath[3] = 16
x_sheath[4] = 12



#xmin = x.min()
#xmax = x.max()
ymin = 0.0
ymax = 6.4e6

##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

sp_temp1=fig.add_subplot(1,1,1)
# ================================= Velocity an Potential=========================================
##============Ref =======
f=h5.File("ref/data_global.h5")

#Velocity
val0 = f["/Fields/Vparallel_global_D1_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
V0[0] = val0_1d[x_sheath[0]]

#Potential
val0 = f["/Fields/Phi_global_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
P0[0] = val0_1d[x_sheath[0]]

#Te
val0 = f["/Fields/T_global_e_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
Te0[0] = val0_1d[x_sheath[0]]

#Ti
val0 = f["/Fields/T_global_D1_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
Ti0[0] = val0_1d[x_sheath[0]]

P0[0] = -P0[0] / (0.5 * Te0[0] * math.log( (2.0*const_pi*mass_e/mass_i)*(1.0+Ti0[0]/Te0[0]) ))
V0[0] = -V0[0] / math.sqrt( const_e * ( Te0[0] + Ti0[0] ) / mass_i )




##============ Case1 =========
t = 13
f=h5.File("Re0.2/data_global.h5")

#Velocity
val0 = f["/Fields/Vparallel_global_D1_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
V0[1] = val0_1d[x_sheath[1]]

#Potential
val0 = f["/Fields/Phi_global_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
P0[1] = val0_1d[x_sheath[1]]

#Te
val0 = f["/Fields/T_global_e_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
Te0[1] = val0_1d[x_sheath[1]]

#Ti
val0 = f["/Fields/T_global_D1_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
Ti0[1] = val0_1d[x_sheath[1]]

P0[1] = -P0[1] / (0.5 * Te0[1] * math.log( (2.0*const_pi*mass_e/mass_i)*(1.0+Ti0[1]/Te0[1]) ))
V0[1] = -V0[1] / math.sqrt( const_e * ( Te0[1] + Ti0[1] ) / mass_i )


##============ Case2 =========
t = 13
f=h5.File("Re0.4/data_global.h5")

#Velocity
val0 = f["/Fields/Vparallel_global_D1_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
V0[2] = val0_1d[x_sheath[2]]

#Potential
val0 = f["/Fields/Phi_global_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
P0[2] = val0_1d[x_sheath[2]]

#Te
val0 = f["/Fields/T_global_e_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
Te0[2] = val0_1d[x_sheath[2]]

#Ti
val0 = f["/Fields/T_global_D1_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
Ti0[2] = val0_1d[x_sheath[2]]

P0[2] = -P0[2] / (0.5 * Te0[2] * math.log( (2.0*const_pi*mass_e/mass_i)*(1.0+Ti0[2]/Te0[2]) ))
V0[2] = -V0[2] / math.sqrt( const_e * ( Te0[2] + Ti0[2] ) / mass_i )




##============ Case3 ===========
t = 19
f=h5.File("Re0.6/data_global.h5")

#Velocity
val0 = f["/Fields/Vparallel_global_D1_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
V0[3] = val0_1d[x_sheath[3]]

#Potential
val0 = f["/Fields/Phi_global_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
P0[3] = val0_1d[x_sheath[3]]

#Te
val0 = f["/Fields/T_global_e_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
Te0[3] = val0_1d[x_sheath[3]]

#Ti
val0 = f["/Fields/T_global_D1_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
Ti0[3] = val0_1d[x_sheath[3]]

P0[3] = -P0[3] / (0.5 * Te0[3] * math.log( (2.0*const_pi*mass_e/mass_i)*(1.0+Ti0[3]/Te0[3]) ))
V0[3] = -V0[3] / math.sqrt( const_e * ( Te0[3] + Ti0[3] ) / mass_i )


##============ Case4 ===========
f=h5.File("Re0.8/data_global.h5")

#Velocity
val0 = f["/Fields/Vparallel_global_D1_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
V0[4] = val0_1d[x_sheath[4]]

#Potential
val0 = f["/Fields/Phi_global_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
P0[4] = val0_1d[x_sheath[4]]

#Te
val0 = f["/Fields/T_global_e_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
Te0[4] = val0_1d[x_sheath[4]]

#Ti
val0 = f["/Fields/T_global_D1_avg"]
val0 = val0[...]
val0_1d = np.transpose(val0[t, 0, 0, :])
Ti0[4] = val0_1d[x_sheath[4]]

P0[4] = -P0[4] / (0.5 * Te0[4] * math.log( (2.0*const_pi*mass_e/mass_i)*(1.0+Ti0[4]/Te0[4]) ))
V0[4] = -V0[4] / math.sqrt( const_e * ( Te0[4] + Ti0[4] ) / mass_i )



print P0
print V0


sp_temp1.yaxis.set_major_formatter(yformatter)
cf_temp1=sp_temp1.plot(x_axis, P0, marker = '8', label = r'$\mathrm{Normalized\ potential}$', linestyle = linestyles[0])
cf_temp1=sp_temp1.plot(x_axis, V0, marker = '^', label = r'$\mathrm{Normalized\ D^+\ ion\ velocity}$', linestyle = linestyles[1])

sp_temp1.set_ylim((0.6, 1.3))
major_ticks = np.arange(0.6, 1.3, 0.2)                                              
#minor_ticks = np.arange(0, 31, 5)                                          
sp_temp1.set_yticks(major_ticks)                                                       
#sp_temp1.set_yticks(minor_ticks, minor=True)  





#major_ticks = np.arange(0.2, 1.02, 0.2)                                              
#minor_ticks = np.arange(0, 31, 5)                                          
#sp_temp1.set_yticks(major_ticks)                                                       
#sp_temp1.set_yticks(minor_ticks, minor=True)  




sp_temp1.grid(True)
sp_temp1.legend(fontsize = legend_fontsize, loc = 7, framealpha=1)

sp_temp1.xaxis.set_ticks(x_axis)

#sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_ylim((ymin, ymax))

sp_temp1.set_xlabel(r'$c_r$', fontsize = label_fontsize)
#sp_temp1.set_ylabel(r'Normalized fluxes', fontsize = label_fontsize)





fig.savefig("all_VP_D.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK


##1d plot=============================
##ax1=fig.add_subplot(2,1,1)
#flux=f["/1d/pflux"]
#flux=flux[...]

#line1,=ax1.plot(flux[0,:],flux[1,:])
#line2,=ax1.plot(flux[0,:],flux[2,:])
