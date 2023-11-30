import hoomd
import hoomd.md
import math
import numpy as np
import datetime
import gsd.hoomd
import itertools
impor numpy

# Define the device
cpu = hoomd.device.CPU()

# Initialize the simulation
sim = hoomd.Simulation(device=cpu, seed=777)

# Parameters
particleNumber = 409
numberDensity = 0.03
temperature = 0.25
timeSteps = int(10e6)


# create a initial box size
box_size = math.pow(particleNumber/numberDensity, 1.0 / 3.0)

# Calculate the box size (different from the homework script)
volume = particleNumber / numberDensity



# Initialize a system with particles placed at random
snapshot = gsd.hoomd.Snapshot()
snapshot.particles.N = particleNumber
snapshot.particles.types = ['A']
snapshot.particles.typeid = np.zeros(particleNumber, dtype=np.uint32)

# initial positions to uniformly distributed
# define the box
x = np.linspace(-box_size / 2, box_size / 2, 10, endpoint=False)
position = np.array(list(itertools.product(x, repeat=3)))
snapshot.particles.position = position[:particleNumber]
snapshot.particles.box = [box_size, box_size, box_size, 0, 0, 0]


# Create snapshot
sim.create_state_from_snapshot(snapshot)

##################
# define the integration

# Parameters (different from the original script)
dt = 1.0E-3
potential_k = 8.00
potential_phi = 0.53

# Set up integrator
integrator = hoomd.md.Integrator(dt=dt)

# Cell list
cell = hoomd.md.nlist.Cell(buffer=0.4)

# Define OPP potential and add to integrator
# default_r_cut need to be modified, or just leave it 3.0! should be an OK guess? (see the paper)
opp = hoomd.md.pair.OPP(cell, default_r_cut=2.5)  # should we change the default_r_cut?

opp.params[('A', 'A')] = {
    'C1': 1., 'C2': 1., 'eta1': 15,
    'eta2': 3, 'k': potential_k, 'phi': potential_phi}
opp.r_cut[('A', 'A')] = 3.0
integrator.forces.append(opp)

# Pair potential table?
# nl = hoomd.md.nlist.Cell(0.01)
# table = hoomd.md.pair.Table(nlist=cell)
# table.params[('A', 'A')] = dict(func=opp, rmin=0.5, rmax=2.5, coeff=dict(k=potential_k, phi=potential_phi))

# Integrate at constant temperature
full = hoomd.filter.All()
nvt = hoomd.md.methods.NVT(kT=temperature, filter=full, tau=1.0)
integrate = hoomd.md.Integrator(dt=0.001, methods=[nvt])
integrator.methods.append(nvt)

# assign to the simulation
sim.operations.integrator = integrator

##################

# thermalize i.e. setting randon velocities (different from the original script, should we modify it?)
sim.state.thermalize_particle_momenta(filter=full, kT=temperature)

# compute the properties of the system
thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=full)
sim.operations.computes.append(thermodynamic_properties)

# make all the properties available to the logger
sim.run(0)

# save the initial configuration to a gsd file
gsd_writer = hoomd.write.GSD(filename='initial.gsd',
                             trigger=hoomd.trigger.Periodic(1000),
                             mode='xb')

# define the logger for tps printing ### (what is tps printing?)
logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(sim, quantities=['timestep', 'tps'])

# not sure what is this for and commented out
# sim.operations.writers.append(hoomd.write.Table(trigger=hoomd.trigger.Periodic(period=100),
#                                                 logger=logger, mode='txt', file='log-output.log'))


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
gsd_writer = hoomd.write.GSD(filename='sample2.gsd',
                             trigger=hoomd.trigger.Periodic(ndump),
                             mode='wb')
sim.operations.writers.append(gsd_writer)

# Run the simulation
sim.run(timeSteps)
