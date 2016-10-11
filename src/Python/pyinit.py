"""@package pyinit
    Definition of Smilei components
"""


class SmileiComponentType(type):
    """Metaclass to all Smilei components"""

    # Constructor of classes
    def __init__(self, name, bases, attrs):
        self._list = []
        self.verify = True
        self.current = 0
        # Run standard metaclass init
        super(SmileiComponentType, self).__init__(name, bases, attrs)

    # Functions to define the iterator
    def __iter__(self):
        return self
    def next(self):
        if self.current >= len(self._list):
            raise StopIteration
        self.current += 1
        return self._list[self.current - 1]

    # Function to return one given instance, for example DiagParticles[0]
    # Special case: species can also be indexed by their name: Species["ion1"]
    def __getitem__(self, key):
        if self.__name__ == "Species" and type(key) is str:
            for obj in self._list:
                if obj.species_type == key:
                    return obj
        else:
            return self._list[key]

    # Function to return the number of instances, for example len(Species)
    def __len__(self):
        return len(self._list)

    # Function to display the content of the component
    def __repr__(self):
        if len(self._list)==0:
            return "<Empty list of "+self.__name__+">"
        else:
            l = []
            for obj in self._list: l.append(str(obj))
            return "["+", ".join(l)+"]"


class SmileiComponent(object):
    """Smilei component generic class"""
    __metaclass__ = SmileiComponentType

    # This constructor is used always for all child classes
    def __init__(self, **kwargs):
        if kwargs is not None: # add all kwargs as internal class variables
            for key, value in kwargs.iteritems():
                if key=="_list":
                    print "Python warning: in "+type(self).__name__+": cannot have argument named '_list'. Discarding."
                else:
                    setattr(self, key, value)
        type(self)._list.append(self) # add the current object to the static list "list"


class Species(SmileiComponent):
    """Species parameters"""
    species_type = None
    initPosition_type = None
    initMomentum_type = ""
    n_part_per_cell = None
    c_part_max = 1.0
    charge_density = None
    nb_density = None
    mean_velocity = [0.]
    temperature = [1e-10]
    dynamics_type = "norm"
    time_frozen = 0.0
    radiating = False
    bc_part_type_west = None
    bc_part_type_east = None
    bc_part_type_north = None
    bc_part_type_south = None
    ionization_model = "none"
    atomic_number = None
    atomic_mass = None
    surface_binding_energy = None
    isTest = False


class Laser(SmileiComponent):
    """Laser parameters"""
    boxSide = None
    a0 = None
    omega0 = 1.
    delta = 1.
    tchirp = 0.
    focus = []
    angle = 0.
    delay = 0.
    time_profile = None
    int_params = []
    double_params = []
    transv_profile = None
    int_params_transv = []
    double_params_transv = []

class Collisions(SmileiComponent):
    """Collisions parameters"""
    // species1 and species2 are list of unsigned integer
    species1 = None
    species2 = None
    coulomb_log = 0.
    debug_every = 0


class PSI(SmileiComponent):
    """PSI parameters"""
    # species1 and species2 are unsigned integer type
    species1 = None
    species2 = None
    PSI_type = None
    posEmit = None
    nPartEmit = 0
    tempEmit = 0.0
    jEmit = 0.0
    weight_const = 1.0
    emitOffset = 0.0
    a_FN = 1.5414e-6
    b_FN = 6.8309e9
    work_function = 4.5



#diagnostics
class DiagProbe(SmileiComponent):
    """Diagnostic probe"""
    every = 0
    time_range = [None, None]
    number = []
    pos = []
    pos_first = []
    pos_second = []
    pos_third = []
    fields = []

class DiagParticles(SmileiComponent):
    """Diagnostic particles"""
    output = None
    every = None
    time_average = 1
    species = None
    axes = []

class DiagPhase(SmileiComponent):
    """Diagnostic phase"""
    pass

class DiagScalar(SmileiComponent):
    """Diagnostic scalar"""
    every = None
    time_range = []
    precision = 10
    vars = []

# external fields
class ExtField(SmileiComponent):
    """External Field"""
    field = []
    profile = None

# default simulation values
output_script = "smilei.py"
dump_step = 0
dump_minutes = 0.0
exit_after_dump = True
restart = False
check_stop_file = False
dump_file_sequence = 2
sim_units = ""
wavelength_SI = 0.
dim = ""
interpolation_order = 2
res_time = None
res_space = []
timestep = None
cell_length = []
sim_time = None
sim_length = []
bc_em_type_x = []
bc_em_type_y = []
time_fields_frozen = 0.0
nspace_win_x = 0
t_move_win = 0.0
vx_win = 1.
clrw = 1
every = 0
number_of_procs = [None]
print_every = None
fieldDump_every = 0
fieldsToDump = []
avgfieldDump_every = None
ntime_step_avg = 0
particleDump_every = None # for backwards-compatibility
time_fields_frozen = 0.
