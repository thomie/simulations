#!/usr/bin/env python

from egfrd import *
from dumper import *
from model import *
import myrandom

'''Reaction network:

A     <-> B
A + B <-> C
C      -> 0
'''


'''Switch all toggles down.

'''
WORLD                        = None
MEMBRANE1                    = None
MEMBRANE2                    = None
DNA                          = None
BOX                          = None
REACTIONS                    = None
DECAY                        = None
SURFACE_BINDING_INTERACTIONS = None
SURFACE_UNBINDING_REACTIONS  = None
LOGGER                       = None
VTK_LOGGER                   = None


'''Toggles.

'''
#WORLD                        = True
#MEMBRANE1                    = True
#MEMBRANE2                    = True
DNA                          = True
#BOX                          = True
#REACTIONS                    = True
#DECAY                        = True
#SURFACE_BINDING_INTERACTIONS = True
#SURFACE_UNBINDING_REACTIONS  = True
#LOGGER                       = True
VTK_LOGGER                   = True


'''Settings

'''
# Make multis run 'BD_DT_FACTOR' times faster than normal.
BD_DT_FACTOR       = 1
RADIUS_FACTOR      = 1
SINGLE_RATE_FACTOR = 1
PAIR_RATE_FACTOR   = 1
N_PARTICLES_FACTOR = 1
WORLD_SIZE_FACTOR  = 1

N_STEPS = 100

MY_SEED = 0
myrandom.seed(MY_SEED)

L = WORLD_SIZE_FACTOR * 1e-6    # Size of simulation box.
D = 1e-12                       # Diffusion constant.
v = 0 # drift
radius = RADIUS_FACTOR * 2.5e-9 # Radius of particles.


'''Reaction rates

'''
sigma = 2 * radius
D_tot = D
tau = sigma * sigma / D_tot
# Bimolecular reaction type.
kf_2 = PAIR_RATE_FACTOR * 100 * sigma * D_tot
# Unimolecular reaction type.
kb_2 = PAIR_RATE_FACTOR * 0.1 / tau
kf_1 = SINGLE_RATE_FACTOR * 0.1 / tau
kb_1 = SINGLE_RATE_FACTOR * 0.1 / tau

# Surface binding interaction (molecule + surface).
kon = kf_2
# Surface unbinding reaction (unimolecular).
koff = kb_1




'''Create model

'''
from model import *
m = ParticleModel(L)


'''Add Surfaces and Regions. 

'''
speciesA = []
speciesB = []
speciesC = []
for i in range(4): #4
    for j in range(5): # 5
        x = L / 2 + (myrandom.random() - 0.5) * L / 3
        #x = 3 * L / 4 + i * L / 20.
        y = - L / 10.
        z =  L / 2 + (myrandom.random() - 0.5) * L / 4
        #z = 3 * L / 4 + j * L / 20.
        d = create_cylindrical_surface('d' + str(i) + str(j),
                                       [x, y, z],
                                       radius,
                                       [0, 1, 0],
                                       L*1.2)
        m.add_structure(d)


        '''Define and add SpeciesTypes.

        '''
        v = myrandom.random()
        Ad = Species('Ad' + str(i), D, radius, d.id, v)
        m.add_species_type(Ad)


'''Define simulator.

A cube of dimension (L x L x L) is created, with periodic boundary 
conditions.

'''
matrix_size = 10
w = create_world(m, matrix_size)
s = EGFRDSimulator(w)
s.bd_dt_factor = BD_DT_FACTOR


for st in m.species_types:
    print st
    throw_in_particles(s.world, st, N_PARTICLES_FACTOR * 1)
    '''Add particles.

    '''
    #throw_in_particles(s.world, Ad, N_PARTICLES_FACTOR * 1)
    #throw_in_particles(s.world, Bd, N_PARTICLES_FACTOR * 1)
    #throw_in_particles(s.world, Cd, N_PARTICLES_FACTOR * 4)


'''Define loggers.

'''
if VTK_LOGGER == True:
    from visualization import vtklogger
    # Write vtk files. See vtklogger.py. Use VTK or Paraview to 
    # visualize.  
    vtklogger = vtklogger.VTKLogger(s, 'data/run', extra_particle_step=False)

if LOGGER == True:
    from logger import Logger
    # Write log files. See logger.py. Use visualizer.py to visualize.
    l = Logger('example')
    l.start(s)
    #l.set_particle_out_interval(1e-4) #1e-3)


'''Run  simulation.

'''
def run():
    print dump_reaction_rules(m)

    s.initialize()
    for i in range(N_STEPS):
        s.step()
        if i % 10 == 0:
            vtklogger.log()

    if VTK_LOGGER == True:
        vtklogger.stop()

    s.print_report()

    s.stop(s.t)


if __name__ == '__main__':
    run()

