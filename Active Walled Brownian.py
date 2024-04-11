##### System variable definitions #####
sigma = 1     # sigma is diameter of particles
eps = 1       # epsilon is a constant affecting Lennard Jones function steepness/energy
kT = 10        # kT is the system thermal energy
dt = 1e-5    # dt is time step size
N = 100       # N is number of particles

##### import python packages
import hoomd
from hoomd import md
import gsd.hoomd
import numpy as np
import datetime, time



##### Setting initial particle positions #####
particlePositions = []     # Initializing an empty array to store particle positions
for i in range(10):        # Initialize square packing of 100 particles [10X10]
  for j in range(10):
    particlePositions.append([1*i-12, 1*j-5, 0])

wall_particles = []
for i in range(16):
    if i < 6 or i > 9:
        wall_particles.append([0, i-15/2, 0])
#Check for overlapping particle positions here too
    wall_particles.append([-15 + i, 15/2, 0])
#Stops overlapping particle positions
    if i < 15:
        wall_particles.append([15, i - 15 / 2, 0])

# Now we make the system in hoomd
hoomd.context.initialize()

# A small shift to help with the periodic box
snap = hoomd.data.make_snapshot(N = N + len(wall_particles),
                                box = hoomd.data.boxdim(Lx=30,
                                                        Ly=15,
                                                        dimensions=2),
                                particle_types = ['A','B'])

# Set positions/types for all particles

snap.particles.position[:] = (particlePositions + wall_particles)[:]
snap.particles.typeid[:] = ([0] * N + [1] * len(wall_particles))[:]
snap.particles.types[:] = ['A','B'][:]

# Initialize the system
system = hoomd.init.read_snapshot(snap)
all = hoomd.group.all()
ljgroup = hoomd.group.type(name="particleA", type='A')
groups = []

#for i in unique_char_types:
#       groups.append(hoomd.group.type(type=i))


# Set particle potentials
groups.append(hoomd.group.type(name= "particleA", type='A'))
groups.append(hoomd.group.type(name= "particleB", type='B'))
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=0.75, nlist=nl)
####Check
lj.set_params(mode='shift')

#for i in range(0, len(unique_char_types)):
#    for j in range(i, len(unique_char_types)):
#        lj.pair_coeff.set(unique_char_types[i],
#                          unique_char_types[j],
#                          epsilon=self.eps, sigma=self.sigma)        # Brownian equilibration

lj.pair_coeff.set('A', 'A', epsilon=eps, sigma=sigma)
lj.pair_coeff.set('A', 'B', epsilon=eps, sigma=sigma)
lj.pair_coeff.set('B', 'B', epsilon=0, sigma=sigma)
brownEquil = 10000
hoomd.md.integrate.mode_standard(dt=dt)

f_lst = [(30, 0, 0)] * len(ljgroup)


#define the forces
force = hoomd.md.force.active(
    seed = 1,
    group=hoomd.group.type(name= "particleA", type='A'),
    f_lst = f_lst
)

#find out how to turn group to equal ljgroup for this integrator
bd = hoomd.md.integrate.brownian(group=hoomd.group.type(name= "particleA", type='A'), kT=kT, seed=1)        # Name the file from parameters        # Name the file from parameters

#name the file
out = "active_brownian_motion"
#out += "_pb" + str(1)
#out += "_phi" + str(1)
#out += "_eps" + str(self.eps)
#out += "_xa" + str(self.partFracA)
#out += "_pNum" + str(self.partNum)
#out += "_dtau" + "{:.1e}".format(dt)
out += ".gsd"

# Write dump
hoomd.dump.gsd(out,
               period=10000,
               group=all,
               overwrite=True,
               phase=-1,
               dynamic=['attribute', 'property', 'momentum'])        # Run
hoomd.run(1000000)