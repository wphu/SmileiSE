# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
#
# Remember: never override the following names:
#           SmileiComponent, Species, Laser, Collisions, DiagProbe, DiagParticles,
#           DiagScalar, DiagPhase or ExtField
#
import math

method = 'explicit'

l0 = 0.5e-5     # nu.norm_l is reference time, the value's unit before / is m (SI)
Lsim = [5000.*l0]	# length of the simulation

t0 = 0.5e-12
us = int( 1.0e-6 / t0 )
Tsim = 20 * us			# duration of the simulation
output_step = 20

# number of processes
n_procs = 24


B = 2.0
Bangle = 5.0

#source parameters
source_dn = 2.0e25
source_temp = 300.0
source_density = 1.0e19

# parameters of recycling
isRecycling = False
cr = 0.8
temp_recycling = 2.0 #eV

# parameters of Carbon flux
isIC = False
flux_C = 1.0e22
temp_C = 1.0  # eV



#> number of timestep of incrementing averaged electromagnetic fields
ntime_step_avg = 100000

#> Timestep to output some fields into hdf5 file
dump_step = int( Tsim / output_step )

timesteps_collision = 40


timesteps_restore = dump_step

collision_zoom_factor = 100.0

ion_step = 1


# dim: Geometry of the simulation
#      1d3v = cartesian grid with 1d in space + 3d in velocity
#      2d3v = cartesian grid with 2d in space + 3d in velocity
#      3d3v = cartesian grid with 3d in space + 3d in velocity
#      2drz = cylindrical (r,z) grid with 3d3v particles
#
dim = '1d3v'

number_of_procs = [n_procs]

#print sim_time / timestep
# ELECTROMAGNETIC BOUNDARY CONDITIONS
# bc_em_type_x/y/z : boundary conditions used for EM fields
#                    periodic = periodic BC (using MPI topology)
#                    silver-muller = injecting/absorbing BC
#                    reflective = consider the ghost-cells as a perfect conductor
#
bc_em_type_x = ['Dirichlet', 'Dirichlet']
#bc_em_type_x = ['Neumann', 'Dirichlet']

bc_em_value_x = [0.0, 0.0]

angle = Bangle * math.pi / 180.0
Bx = B * math.sin(angle)
By = B * math.cos(angle)
Bz = 0.0
externB = [Bx, By, Bz]

ion_sound_velocity = 0.0   #math.sqrt( (20.0 * 1.6021766208e-19) / (2.0 * 1.67262158e-27) )
vx = ion_sound_velocity * math.sin(angle)
vy = ion_sound_velocity * math.cos(angle)
vz = 0.0




# RANDOM seed
# this is used to randomize the random number generator
random_seed = 0


# order of interpolation
#
interpolation_order = 1
projection_order = 2


# SIMULATION BOX : for all space directions (use vector)
# cell_length: length of the cell
# sim_length: length of the simulation in units of the normalization wavelength
#
cell_length = [l0]
sim_length  = Lsim

# SIMULATION TIME
# timestep: duration of the timestep
# sim_time: duration of the simulation in units of the normalization period
#
timestep = t0
n_time = Tsim





# ================= DEFINE ALL SPECIES ===========================================
# species_type       = string, given name to the species (e.g. ion, electron, positron, test ...)
# initPosition_type  = string, "regular" or "random"
# initMomentum_type  = string "cold", "maxwell-juettner" or "rectangular"
# c_part_max         = float, factor on the memory reserved for the total number of particles
# mass               = float, particle mass in units of the electron mass
# dynamics_type      = string, type of species dynamics = "norm" or "rrLL"
# time_frozen        = float, time during which particles are frozen in units of the normalization time
# radiating          = boolean, if true, incoherent radiation calculated using the Larmor formula
# n_part_per_cell    = integer or function, number of particles/cell
# charge             = float or function, particle charge in units of the electron charge
# charge_density     = float or function, species charge density in units of the "critical" density
#     or nb_density for number density
# mean_velocity      = list of floats or functions, mean velocity in units of the speed of light
# temperature        = list of floats or functions, temperature in units of m_e c^2
# Predefined functions: constant, trapezoidal, gaussian, polygonal, cosine
#

Species(
	species_type = 'e',
	initPosition_type = 'random',
	initMomentum_type = 'maxwell',
	ionization_model = 'none',
	#Pusher_type = 'GC0',
	n_part_per_cell = 0,
	n_part_per_cell_for_weight = 160,
	c_part_max = 1.0,
	mass = 9.109382616e-31,
	charge = -1.6021766208e-19,
	nb_density = 2.0e19,
	temperature = [20.0],
	time_frozen = 0.,
	bc_part_type_west  = 'supp',
	bc_part_type_east  = 'supp',
)

Species(
	species_type = 'D1',
	initPosition_type = 'random',
	initMomentum_type = 'maxwell',
	ionization_model = 'none',
	timestep_zoom = ion_step,
	n_part_per_cell = 0,
	n_part_per_cell_for_weight = 160,
	c_part_max = 1.0,
	mass = 2.0 * 1.67262158e-27,
	charge = 1.6021766208e-19,
	nb_density = 2.0e19,
	temperature = [20.0],
	time_frozen = 0.0,
	bc_part_type_west  = 'supp',
	bc_part_type_east  = 'supp',

	diameter = 2.751E-10,
	ref_temperature = 273.,
	visc_temp_index = 0.75,
	vss_scat_inv = 1.
)

if isRecycling:
	Species(
		species_type = 'D',
		initPosition_type = 'random',
		initMomentum_type = 'maxwell',
		ionization_model = 'none',
		timestep_zoom = ion_step,
		n_part_per_cell = 0,
		n_part_per_cell_for_weight = 160,
		c_part_max = 1.0,
		mass = 2.0 * 1.67262158e-27,
		charge = 0.0,
		nb_density = 2.0e19,
		temperature = [20.0],
		time_frozen = 0.0,
		bc_part_type_west  = 'supp',
		bc_part_type_east  = 'supp',

		diameter = 2.751E-10,
		ref_temperature = 273.,
		visc_temp_index = 0.75,
		vss_scat_inv = 1.
	)

if isIC:
	Species(
		species_type = 'C',
		initPosition_type = 'random',
		initMomentum_type = 'rectangular',
		ionization_model = 'none',
		timestep_zoom = ion_step,
		n_part_per_cell = 0,
		n_part_per_cell_for_weight = 160,
		c_part_max = 1.0,
		mass = 1.993e-26,
		charge = 0.0,
		nb_density = 2.0e19,
		temperature = [20.0],
		mean_velocity = [vx, vy, vz],
		time_frozen = 0.0,
		bc_part_type_west  = 'supp',
		bc_part_type_east  = 'supp',

		diameter = 3.784E-10,
		ref_temperature = 273.,
		visc_temp_index = 0.75,
		vss_scat_inv = 1.
	)


	Species(
		species_type = 'C1',
		initPosition_type = 'random',
		initMomentum_type = 'rectangular',
		ionization_model = 'none',
		timestep_zoom = ion_step,
		n_part_per_cell = 0,
		n_part_per_cell_for_weight = 160,
		c_part_max = 1.0,
		mass = 1.993e-26,
		charge = 1.6021766208e-19,
		nb_density = 2.0e19,
		temperature = [20.0],
		mean_velocity = [vx, vy, vz],
		time_frozen = 0.0,
		bc_part_type_west  = 'supp',
		bc_part_type_east  = 'supp',

		diameter = 3.784E-10,
		ref_temperature = 273.,
		visc_temp_index = 0.75,
		vss_scat_inv = 1.
	)


	Species(
		species_type = 'C2',
		initPosition_type = 'random',
		initMomentum_type = 'rectangular',
		ionization_model = 'none',
		timestep_zoom = ion_step,
		n_part_per_cell = 0,
		n_part_per_cell_for_weight = 160,
		c_part_max = 1.0,
		mass = 1.993e-26,
		charge = 2.0 * 1.6021766208e-19,
		nb_density = 2.0e19,
		temperature = [20.0],
		mean_velocity = [vx, vy, vz],
		time_frozen = 0.0,
		bc_part_type_west  = 'supp',
		bc_part_type_east  = 'supp',

		diameter = 3.784E-10,
		ref_temperature = 273.,
		visc_temp_index = 0.75,
		vss_scat_inv = 1.
	)


	Species(
		species_type = 'C3',
		initPosition_type = 'random',
		initMomentum_type = 'rectangular',
		ionization_model = 'none',
		timestep_zoom = ion_step,
		n_part_per_cell = 0,
		n_part_per_cell_for_weight = 160,
		c_part_max = 1.0,
		mass = 1.993e-26,
		charge = 3.0 * 1.6021766208e-19,
		nb_density = 2.0e19,
		temperature = [20.0],
		mean_velocity = [vx, vy, vz],
		time_frozen = 0.0,
		bc_part_type_west  = 'supp',
		bc_part_type_east  = 'supp',

		diameter = 3.784E-10,
		ref_temperature = 273.,
		visc_temp_index = 0.75,
		vss_scat_inv = 1.
	)



# ================ PartSource ========================================
PartSource(
	species1 = ["e"],
	PartSource_type = "Load",
	loadKind = "nq",
	loadNumber = 4,
	loadDensity = source_density,
	loadTemperature = source_temp,
	loadTemperature_upLimit_factor = 2.1,
	loadDn = source_dn,
	loadPos_start 	= Lsim[0] / 2.0 - 500.0*l0,
	loadPos_end 	= Lsim[0] / 2.0 + 500.0*l0,
	species_dependent = ["e"],

)


PartSource(
	species1 = ["D1"],
	PartSource_type = "Load",
	loadKind = "nq",
	loadNumber = 4,
	loadDensity = source_density,
	loadTemperature = source_temp,
	loadTemperature_upLimit_factor = 2.1,
	loadDn = source_dn,
	loadPos_start = 	Lsim[0] / 2.0 - 500.0*l0,
	loadPos_end = 		Lsim[0] / 2.0 + 500.0*l0,
	species_dependent = ["e"],

)

if isIC:
	PartSource(
		PartSource_type = "Emit",
		species1 = ["C"],
		emitKind = "regular",
		emitPos = "left",
		emitNumber = 2,
		emitOffset = 0.2,
		emitTemp = temp_C,
		emitFlux = flux_C
	)

	PartSource(
		PartSource_type = "Emit",
		species1 = ["C"],
		emitKind = "regular",
		emitPos = "right",
		emitNumber = 2,
		emitOffset = 0.2,
		emitTemp = temp_C,
		emitFlux = flux_C
	)



#==================== PSI =================================
if isRecycling:
	PSI(
		species1 = ["D1"],
		species2 = ["D"],
		PSI_type = "Recycling",
		psiPos = "left",
		emitTemp = temp_recycling,
		recycling_factor = cr
	)

	PSI(
		species1 = ["D1"],
		species2 = ["D"],
		PSI_type = "Recycling",
		psiPos = "right",
		emitTemp = temp_recycling,
		recycling_factor = cr
	)



#==================== Collisions =================================

'''
# DSMC
Collisions(
	species1 = ["D", "D1"],
	species2 = ["C" ],
	collisions_type = "DSMC"
)
'''


# Coulomb collisions
if isIC:
	Collisions(
		species1 = ["e", "D1", "C1", "C2", "C3"],
		#coulomb_log = 1,
		collisions_type = "coulomb"
	)
else:
	Collisions(
		species1 = ["e", "D1"],
		#coulomb_log = 1,
		collisions_type = "coulomb"
	)



# For D ==================
if isRecycling:
	Collisions(
		collisions_type = "Ionization",
		species1 = ["e"],
		species2 = ["D"],
		species3 = ["D1"],
		crossSection_fileName = "crossSection/Ionization_D_to_D+1.dat"
	)

	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["D"],
		crossSection_fileName = "crossSection/Excitation_H_1s-2p.dat"
	)

	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["D"],
		crossSection_fileName = "crossSection/Excitation_H_1s-3p.dat"
	)

	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["D"],
		crossSection_fileName = "crossSection/Excitation_H_1s-4p.dat"
	)




# For C ==========
if isIC:
	Collisions(
		collisions_type = "Ionization",
		species1 = ["e"],
		species2 = ["C"],
		species3 = ["C1"],
		crossSection_fileName = "crossSection/Ionization_C_to_C+1.dat"
	)

	Collisions(
		collisions_type = "Ionization",
		species1 = ["e"],
		species2 = ["C1"],
		species3 = ["C2"],
		crossSection_fileName = "crossSection/Ionization_C+1_to_C+2.dat"
	)

	Collisions(
		collisions_type = "Ionization",
		species1 = ["e"],
		species2 = ["C2"],
		species3 = ["C3"],
		crossSection_fileName = "crossSection/Ionization_C+2_to_C+3.dat"
	)

	# For C Excitation
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C"],
		crossSection_fileName = "crossSection/Excitation_C_2s22p21D-2s2p33D.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C"],
		crossSection_fileName = "crossSection/Excitation_C_2s22p21D-2s22p3s1P.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C"],
		crossSection_fileName = "crossSection/Excitation_C_2s22p23P-2s2p33D.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C"],
		crossSection_fileName = "crossSection/Excitation_C_2s22p23P-2s2p33P.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C"],
		crossSection_fileName = "crossSection/Excitation_C_2s22p23P-2s2p33S.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C"],
		crossSection_fileName = "crossSection/Excitation_C_2s22p23P-2s2p35S.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C"],
		crossSection_fileName = "crossSection/Excitation_C_2s22p23P-2s22p3d3D.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C"],
		crossSection_fileName = "crossSection/Excitation_C_2s22p23P-2s22p3p3P.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C"],
		crossSection_fileName = "crossSection/Excitation_C_2s22p23P-2s22p3s3P.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C"],
		crossSection_fileName = "crossSection/Excitation_C_2s22p23P-2s22p21D.dat"
	)


	# For C+1 Excitation
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C1"],
		crossSection_fileName = "crossSection/Excitation_C+1_2s22p2P-2s2p22D.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C1"],
		crossSection_fileName = "crossSection/Excitation_C+1_2s22p2P-2s2p22P.dat"
	)

	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C1"],
		crossSection_fileName = "crossSection/Excitation_C+1_2s22p2P-2s2p22S.dat"
	)

	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C1"],
		crossSection_fileName = "crossSection/Excitation_C+1_2s22p2P-2s2p24P.dat"
	)

	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C1"],
		crossSection_fileName = "crossSection/Excitation_C+1_2s22p2P-2s23d2D.dat"
	)

	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C1"],
		crossSection_fileName = "crossSection/Excitation_C+1_2s22p2P-2s23p2P.dat"
	)

	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C1"],
		crossSection_fileName = "crossSection/Excitation_C+1_2s22p2P-2s23s2S.dat"
	)

	# For C+2 Excitation
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C2"],
		crossSection_fileName = "crossSection/Excitation_C+2_2p23P-2p21D.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C2"],
		crossSection_fileName = "crossSection/Excitation_C+2_2s2p1P-2p21D.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C2"],
		crossSection_fileName = "crossSection/Excitation_C+2_2s2p3P-2p23P.dat"
	)


	# For C+3 Excitation
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C3"],
		crossSection_fileName = "crossSection/Excitation_C+3_2p2P-3d2D.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C3"],
		crossSection_fileName = "crossSection/Excitation_C+3_2p2P-3s2S.dat"
	)
	Collisions(
		collisions_type = "Excitation",
		species1 = ["e"],
		species2 = ["C3"],
		crossSection_fileName = "crossSection/Excitation_C+3_2s2S-2p2P.dat"
	)
