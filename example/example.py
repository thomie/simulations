#!/usr/bin/env python

'''
To run this script:
PYTHONPATH=../../ python features.py

or, without terminal output and faster:
PYTHONPATH=../../../ python -O features.py
'''

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
N_PARTICLES_FACTOR = 0.3
WORLD_SIZE_FACTOR  = 2

N_STEPS = 200

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
if MEMBRANE1 and MEMBRANE2 and BOX:
    b = create_cuboidal_region('b',
                               [0, 0, 2 * L / 10],
                               [L, L, 6 * L / 10])
    m.add_structure(b)

if MEMBRANE1:
    m1 = create_planar_surface('m1',
                               [0, 0, 2 * L / 10],
                               [1, 0, 0],
                               [0, 1, 0],
                               L,
                               L)
    m.add_structure(m1)

if MEMBRANE2:
    m2 = create_planar_surface('m2',
                               [0, 0, 8 * L / 10],
                               [1, 0, 0],
                               [0, 1, 0],
                               L,
                               L)
    m.add_structure(m2)

if DNA:
    d = create_cylindrical_surface('d',
                                   [0, 0, L / 2],
                                   radius,
                                   [0, 1, 0],
                                   L)
    m.add_structure(d)


'''Define and add SpeciesTypes.

'''
if WORLD:
    A = Species('A', D, radius)
    B = Species('B', D, radius)
    C = Species('C', D, radius)
    m.add_species_type(A)
    m.add_species_type(B)
    m.add_species_type(C)

if MEMBRANE1 and MEMBRANE2 and BOX:
    Ab = Species('Ab', D, radius, b.id)
    Bb = Species('Bb', D, radius, b.id)
    Cb = Species('Cb', D, radius, b.id)
    m.add_species_type(Ab)
    m.add_species_type(Bb)
    m.add_species_type(Cb)

if MEMBRANE1:
    Am1 = Species('Am1', D, radius, m1.id)
    Bm1 = Species('Bm1', D, radius, m1.id)
    Cm1 = Species('Cm1', D, radius, m1.id)
    m.add_species_type(Am1)
    m.add_species_type(Bm1)
    m.add_species_type(Cm1)

if MEMBRANE2:
    pass

if DNA:
    Ad = Species('Ad', D, radius, d.id, v)
    Bd = Species('Bd', D, radius, d.id, v)
    Cd = Species('Cd', D, radius, d.id, v)
    m.add_species_type(Ad)
    m.add_species_type(Bd)
    m.add_species_type(Cd)


'''Define simulator.

A cube of dimension (L x L x L) is created, with periodic boundary 
conditions.

'''
matrix_size = 3
w = create_world(m, matrix_size)
s = EGFRDSimulator(w)
s.bd_dt_factor = BD_DT_FACTOR


'''Stir (optional).

'''


'''Add reaction rules.

'''
if REACTIONS and WORLD:
    r1 = create_unimolecular_reaction_rule(A, B, kf_1)
    m.network_rules.add_reaction_rule(r1)

    r2 = create_unimolecular_reaction_rule(B, A, kb_1)
    m.network_rules.add_reaction_rule(r2)

    r3 = create_binding_reaction_rule(A, B, C, kf_2)
    m.network_rules.add_reaction_rule(r3)

    r4 = create_unbinding_reaction_rule(C, A, B, kb_2)
    m.network_rules.add_reaction_rule(r4)

    if DECAY:
        r5 = create_decay_reaction_rule(C, kf_1)
        m.network_rules.add_reaction_rule(r5)

if REACTIONS and MEMBRANE1 and MEMBRANE2 and BOX:
    r1 = create_unimolecular_reaction_rule(Ab, Bb, kf_1)
    m.network_rules.add_reaction_rule(r1)

    r2 = create_unimolecular_reaction_rule(Bb, Ab, kb_1)
    m.network_rules.add_reaction_rule(r2)

    r3 = create_binding_reaction_rule(Ab, Bb, Cb, kf_2)
    m.network_rules.add_reaction_rule(r3)

    r4 = create_unbinding_reaction_rule(Cb, Ab, Bb, kb_2)
    m.network_rules.add_reaction_rule(r4)

    if DECAY:
        r5 = create_decay_reaction_rule(Cb, kf_1)
        m.network_rules.add_reaction_rule(r5)


if REACTIONS and MEMBRANE1:
    r1 = create_unimolecular_reaction_rule(Am1, Bm1, kf_1)
    m.network_rules.add_reaction_rule(r1)

    r2 = create_unimolecular_reaction_rule(Bm1, Am1, kb_1)
    m.network_rules.add_reaction_rule(r2)

    r3 = create_binding_reaction_rule(Am1, Bm1, Cm1, kf_2)
    m.network_rules.add_reaction_rule(r3)

    r4 = create_unbinding_reaction_rule(Cm1, Am1, Bm1, kb_2)
    m.network_rules.add_reaction_rule(r4)

    if DECAY:
        r5 = create_decay_reaction_rule(Cm1, kf_1)
        m.network_rules.add_reaction_rule(r5)

if REACTIONS and MEMBRANE2:
    pass

if REACTIONS and DNA:
    r1 = create_unimolecular_reaction_rule(Ad, Bd, kf_1)
    m.network_rules.add_reaction_rule(r1)

    r2 = create_unimolecular_reaction_rule(Bd, Ad, kb_1)
    m.network_rules.add_reaction_rule(r2)

    r3 = create_binding_reaction_rule(Ad, Bd, Cd, kf_2)
    m.network_rules.add_reaction_rule(r3)

    r4 = create_unbinding_reaction_rule(Cd, Ad, Bd, kb_2)
    m.network_rules.add_reaction_rule(r4)

    if DECAY:
        r5 = create_decay_reaction_rule(Cd, kf_1)
        m.network_rules.add_reaction_rule(r5)


'''Add particles.

'''
if MEMBRANE1 and MEMBRANE2 and BOX:
    throw_in_particles(s.world, Ab, N_PARTICLES_FACTOR * 4)
    throw_in_particles(s.world, Bb, N_PARTICLES_FACTOR * 4)
    throw_in_particles(s.world, Cb, N_PARTICLES_FACTOR * 16)
elif WORLD:
    throw_in_particles(s.world, A, N_PARTICLES_FACTOR * 4)
    throw_in_particles(s.world, B, N_PARTICLES_FACTOR * 4)
    throw_in_particles(s.world, C, N_PARTICLES_FACTOR * 16)

if MEMBRANE1:
    throw_in_particles(s.world, Am1, N_PARTICLES_FACTOR * 2)
    throw_in_particles(s.world, Bm1, N_PARTICLES_FACTOR * 2)
    throw_in_particles(s.world, Cm1, N_PARTICLES_FACTOR * 8)

if MEMBRANE2:
    pass

if DNA:
    throw_in_particles(s.world, Ad, N_PARTICLES_FACTOR * 1)
    throw_in_particles(s.world, Bd, N_PARTICLES_FACTOR * 1)
    throw_in_particles(s.world, Cd, N_PARTICLES_FACTOR * 4)


'''Define loggers.

'''
if VTK_LOGGER == True:
    from visualization import vtklogger
    # Write vtk files. See vtklogger.py. Use VTK or Paraview to 
    # visualize.  
    vtklogger = vtklogger.VTKLogger(s, 'data/run', extra_particle_step=True)

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
        try:
            if VTK_LOGGER == True:
                vtklogger.log()
            if LOGGER == True and s.last_reaction:
                l.log(s, s.t)
            s.step()
        except Exception, message:
            print 'Exception at step = ', i
            print message
            if VTK_LOGGER == True:
                vtklogger.stop()
            s.print_report()
            raise

    if VTK_LOGGER == True:
        vtklogger.stop()

    s.print_report()


    #Todo: add test for this stuff, because it keeps breaking.
    '''
    for x in dump_domains(s):
        print x
    for x in dump_species(s):
        print x
        for y in dump_particles_by_sid(s, x.id):
            print y
    for x in dump_particles(s):
        print x
    '''

    s.stop(s.t)


if __name__ == '__main__':
    run()

