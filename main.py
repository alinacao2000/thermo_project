import hoomd
import hoomd.md
import math
import numpy as np

# Parameters
particleNumber = 4096
numberDensity = 0.03
temperature = 0.25
potential_k = 8.00
potential_phi = 0.53
timeSteps = int(50e6)

# Calculate the box size
volume = particleNumber / numberDensity
box_size = volume ** (1. / 3.)

# Initialize a system with particles placed at random
snapshot = hoomd.Snapshot()
snapshot.particles.N = particleNumber
snapshot.particles.types = ['A']
snapshot.configuration.box = [box_size, box_size, box_size, 0, 0, 0]
snapshot.particles.position[:] = np.random.uniform(-box_size / 2., box_size / 2., size=(particleNumber, 3))

# Initialize the simulation
sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
sim.create_state_from_snapshot(snapshot)

# Pair potential table
nl = hoomd.md.nlist.Cell(0.01)
table = hoomd.md.pair.Table(nlist=nl)
table.params[('A', 'A')] = dict(func=OPP, rmin=0.5, rmax=2.5, coeff=dict(k=potential_k, phi=potential_phi))

# Integrate at constant temperature
nvt = hoomd.md.methods.NVT(kT=temperature, filter=hoomd.filter.All(), tau=1.0)
integrate = hoomd.md.Integrator(dt=0.005, methods=[nvt])
sim.operations.integrator = integrate

# Log the simulation data
logger = hoomd.logging.Logger()
logger.add(sim, quantities=['potential_energy', 'temperature'], list=True)
sim.operations.writers.append(hoomd.write.Table(trigger=hoomd.trigger.Periodic(period=100), logger=logger, mode='txt', file='log-output.log'))

# Run the simulation
sim.run(timeSteps)
