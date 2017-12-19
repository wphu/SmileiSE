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
Lsim = [500.*l0]	# length of the simulation

t0 = 0.5e-12
ns = int( 1.0e-9 / t0 )
Tsim = 30 * ns			# duration of the simulation
output_step = 10

# number of processes
n_procs = 5



#> number of timestep of incrementing averaged electromagnetic fields
ntime_step_avg = ns

#> Timestep to output some fields into hdf5 file
dump_step = int( Tsim / output_step )


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

bc_em_value_x = [60.0, 0.0]



Bx = 0.0
By = 0.0
Bz = 0.0
externB = [Bx, By, Bz]

ion_sound_velocity = math.sqrt( (20.0 * 1.6021766208e-19) / (2.0 * 1.67262158e-27) )
vx = ion_sound_velocity
vy = 0.0
vz = 0.0




# RANDOM seed
# this is used to randomize the random number generator
random_seed = 0


# order of interpolation
#
interpolation_order = 1
projection_order = 1


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





# DEFINE ALL SPECIES
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
	n_part_per_cell = 100,
	n_part_per_cell_for_weight = 100,
	c_part_max = 1.0,
	mass = 9.109382616e-31,
	charge = -1.6021766208e-19,
	nb_density = 1.0e19,
	temperature = [20.0],
	mean_velocity = [vx, vy, vz],
	#mean_velocity = [0.0, 0.0, 0.0],
	time_frozen = 0.,
	bc_part_type_west  = 'supp',
	bc_part_type_east  = 'supp',
)

Species(
	species_type = 'D1',
	initPosition_type = 'random',
	initMomentum_type = 'maxwell',
	ionization_model = 'none',
	n_part_per_cell = 100,
	n_part_per_cell_for_weight = 100,
	c_part_max = 1.0,
	mass = 2.0 * 1.67262158e-27,
	charge = 1.6021766208e-19,
	nb_density = 1.0e19,
	temperature = [20.0],
	mean_velocity = [vx, vy, vz],
	time_frozen = 0.0,
	bc_part_type_west  = 'supp',
	bc_part_type_east  = 'supp',
)



# PartSource
# species1    = list of strings, the names of the first species that performPSI
# species2    = list of strings, the names of the second species that performPSI
#               (can be the same as species1)

nx_source_left = 200

PartSource(
	species1 = ["e"],
	PartSource_type = "Load",
	loadKind = "nT",
	loadStep = 100,
	loadDensity = 1.0e19,
	loadTemperature = 20.0,
	mean_velocity = [vx, vy ,vz],
	#loadDn = 2.0e25,
	loadPos_start 	= 0.0,
	loadPos_end 	= nx_source_left * l0,



)


PartSource(
	species1 = ["D1"],
	PartSource_type = "Load",
	loadKind = "nT",
	loadStep = 100,
	loadDensity = 1.0e19,
	loadTemperature = 20.0,
	mean_velocity = [vx, vy ,vz],
	#loadDn = 2.0e25,
	loadPos_start 	= 0.0,
	loadPos_end 	= nx_source_left * l0,

)
