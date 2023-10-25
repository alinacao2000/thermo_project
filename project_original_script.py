# A simple molecular dynamics simulation script for the HOOMD-blue package # (available for download at: http://codeblue.umich.edu/hoomd-blue/)
# Purpose: Self-assembles an icosahedral quasicrystal.
# This script is part of the Supplementary Information of:
# M. Engel, P.F. Damasceno, C.L. Phillips, S.C. Glotzer
# "Computational self-assembly of a one-component icosahedral quasicrystal" # Nature Materials, doi: 10.1038/nmat4152
import hoomd
import hoomd.md
import math
import numpy as np

# Parameters
particleNumber = 4096
numberDensity = 0.03
temperature = 0.25

potential_k = 8.00
potential_phi = 0.53  # initial k and phi
timeSteps = int(50e6)
# Define the OPP
def OPP(r, rmin, rmax, k, phi):
    cos = math.cos(k * (r - 1.25) - phi)
    sin = math.sin(k * (r - 1.25) - phi)
    V = pow(r, -15) + cos * pow(r, -3)
    F = 15.0 * pow(r, -16) + 3.0 * cos * pow(r, -4) + k * sin * pow(r, -3) return (V, F)

# Determine the potential range by searching for extrema def determineRange(k, phi):
r = 0.5
extremaNum = 0
force1 = OPP(r, 0, 0, k, phi)[1] while (extremaNum < 6 and r < 5.0):
r += 1e-5
force2 = OPP(r, 0, 0, k, phi)[1] if (force1 * force2 < 0.0):
extremaNum += 1 force1 = force2
return r
# Initialize a system with particles placed at random with a given density system = init.create_random(N = particleNumber, phi_p = numberDensity)
# Generate the pair interaction table
range = determineRange(potential_k, potential_phi)
table = pair.table(width = 1000)
table.pair_coeff.set('A', 'A', func = OPP, rmin = 0.5, rmax = range,
coeff = dict(k = potential_k, phi = potential_phi))
# Start logging
filename = "quasicrystal_k" + str(potential_k) + "_phi" + str(potential_phi) filename += "_T" + str(temperature)
dump.xml(filename = filename, period = timeSteps * 1e-2)
logger = analyze.log(filename = filename + ".log", period = timeSteps * 1e-4,
quantities = ['time','potential_energy','pressure'])
# Integrate at constant temperature
integrate.nvt(group = group.all(), tau = 1.0, T = temperature) integrate.mode_standard(dt = 0.01)
run(timeSteps + 1)