### Initialize ###

# Import modules
import datetime
import hoomd
import numpy as np
import math
import gsd.hoomd
import itertools

# Define simulation device
cpu = hoomd.device.CPU()

# Create simulation object
sim = hoomd.Simulation(device=cpu, seed=777)

################

### Define initial condition ###

# Parameters
particleNumber = 409
numberDensity = 0.03
temperature = 0.25
timeSteps = int(10e6)

# create a initial box size
box_size = math.pow(particleNumber/numberDensity, 1.0 / 3.0)

# Calculate the box size (different from the homework script)
volume = particleNumber / numberDensity

x = np.linspace(-box_size / 2, box_size / 2, 10, endpoint=False)
position = np.array(list(itertools.product(x, repeat=3)))
position = position[:particleNumber]

# Define snapshot
snapshot = gsd.hoomd.Snapshot()
snapshot.particles.N = len(position)
snapshot.particles.position = position
snapshot.particles.typeid = [0]*len(position)
snapshot.configuration.box = [box_size, box_size, box_size, 0, 0, 0]
snapshot.particles.types = ['A']

# Create snapshot
sim.create_state_from_snapshot(snapshot)

################################

## Define integration ###

# Parameters
kbt = 0.25
dt = 1.0E-3

# Set up integrator
integrator = hoomd.md.Integrator(dt=dt)

# Cell list
cell = hoomd.md.nlist.Cell(buffer=0.4)

potential_k = 8.00
potential_phi = 0.53
opp = hoomd.md.pair.OPP(cell)  # should we change the default_r_cut?
opp.params[('A', 'A')] = {
    'C1': 1., 'C2': 1., 'eta1': 15,
    'eta2': 3, 'k': potential_k, 'phi': potential_phi}
opp.r_cut[('A', 'A')] = 3.0
integrator.forces.append(opp)

# Define ensemble and add to integrator
full = hoomd.filter.All()
nvt = hoomd.md.methods.NVT(kT=kbt, filter=full, tau=1.0)
integrator.methods.append(nvt)

# Assign to simulation
sim.operations.integrator = integrator

# Thermalize
sim.state.thermalize_particle_momenta(filter=full, kT=kbt)

### Define logger for tps printing ###
logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(sim, quantities=['timestep', 'tps'])
class Status():
    def __init__(self, sim):
        self.sim = sim
    @property
    def seconds_remaining(self):
        try:
            return (self.sim.final_timestep - self.sim.timestep) / self.sim.tps
        except ZeroDivisionError:
            return 0
    @property
    def etr(self):
        return str(datetime.timedelta(seconds=self.seconds_remaining))
status = Status(sim)
logger[('Status', 'etr')] = (status, 'etr', 'string')
table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(period=5000),
                          logger=logger)
sim.operations.writers.append(table)
######################################


# Write to file
ndump = 1000
gsd_writer = hoomd.write.GSD(filename='sample.gsd',
                             trigger=hoomd.trigger.Periodic(ndump),
                             mode='wb')
sim.operations.writers.append(gsd_writer)

# Run
sim.run(1E6)

#########################

