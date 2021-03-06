#!/usr/bin/env python

'''
To run this script:
PYTHONPATH=../../ python features.py

or, without terminal output and faster:
PYTHONPATH=../../../ python -O features.py
'''

from egfrd import *
from dumper import *
import myrandom


'''Reaction network:

A     <-> B
A + B <-> C
C      -> 0
'''


'''Switch all toggles down.

'''
EGFRD_MODEL                  = None
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
EGFRD_MODEL                  = True
WORLD                        = True
#MEMBRANE1                    = True
#MEMBRANE2                    = True
#DNA                          = True
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
RADIUS_FACTOR      = 2
SINGLE_RATE_FACTOR = 1
PAIR_RATE_FACTOR   = 1
N_PARTICLES_FACTOR = 4
WORLD_SIZE_FACTOR  = 1

N_STEPS = 1500

MY_SEED = 0
myrandom.seed(MY_SEED)

L = WORLD_SIZE_FACTOR * 1e-6    # Size of simulation box.
D = 1e-12                       # Diffusion constant.
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
if EGFRD_MODEL == True:
    from egfrd_model import *
    m = EGFRDModel(L)
else:
    from model import *
    m = ParticleModel(L)


if EGFRD_MODEL == True:
    '''Add Regions and Surfaces. 

    '''
    if BOX:
        b = m.add_cuboidal_region(id='b',
                                  corner=[0, 0, 2 * L / 10],
                                  diagonal=[L, L, 6 * L / 10])

    if MEMBRANE1:
        m1 = m.add_planar_surface(id='m1',
                                  corner=[0, 0, 2 * L / 10],
                                  unit_x=[1, 0, 0],
                                  unit_y=[0, 1, 0])

    if MEMBRANE2:
        m2 = m.add_planar_surface(id='m2',
                                  corner=[0, 0, 8 * L / 10],
                                  unit_x=[1, 0, 0],
                                  unit_y=[0, 1, 0])

    if DNA:
        d = m.add_cylindrical_surface(id='d',
                                      corner=[0, 0, L / 2],
                                      radius=radius,
                                      orientation=[0, 1, 0])


    '''Define ParticleTypes.

    '''
    A = ParticleType('A', D, radius)
    A2 = ParticleType('A2', D, radius)
    B = ParticleType('B', D, radius)
    B2 = ParticleType('B2', D, radius)
    C = ParticleType('C', D, radius)
    C2 = ParticleType('C2', D, radius)
    D2 = ParticleType('D2', D, radius)


    '''Add ParticleTypes (to specific Regions or Surfaces).

    '''
    if WORLD:
        m.add_particle_type(A)
        m.add_particle_type(A2)
        m.add_particle_type(B)
        m.add_particle_type(B2)
        m.add_particle_type(C)
        m.add_particle_type(C2)
        m.add_particle_type(D2)

    if BOX:
        m.add_particle_type(A, b)
        m.add_particle_type(B, b)
        m.add_particle_type(C, b)

    if MEMBRANE1:
        m.add_particle_type(A, m1)
        m.add_particle_type(B, m1)
        m.add_particle_type(C, m1)

    if MEMBRANE2:
        pass

    if DNA:
        m.add_particle_type(A, d, v=10)
        m.add_particle_type(B, d, v=-10)
        m.add_particle_type(C, d)


    '''Add reaction rules.

    '''
    if REACTIONS and WORLD:
        m.add_reaction_rule(A,      B,      kf_1)
        m.add_reaction_rule(B,      A,      kb_1)
        m.add_reaction_rule([A, B], C,      kf_2)
        m.add_reaction_rule(C,      [A, B], kb_2)
        if DECAY:
            m.add_reaction_rule(C, [], kf_1)

    if REACTIONS and BOX:
        m.add_reaction_rule((A, b),           (B, b),           kf_1)
        m.add_reaction_rule((B, b),           (A, b),           kb_1)
        m.add_reaction_rule([(A, b), (B, b)], (C, b),           kf_2)
        m.add_reaction_rule((C, b),           [(A, b), (B, b)], kb_2)
        if DECAY:
            m.add_reaction_rule([(C, b)], [], kf_1)

    if REACTIONS and MEMBRANE1:
        m.add_reaction_rule((A, m1),            (B, m1),            kf_1)
        m.add_reaction_rule((B, m1),            (A, m1),            kb_1)
        m.add_reaction_rule([(A, m1), (B, m1)], (C, m1),            kf_2)
        m.add_reaction_rule((C, m1),            [(A, m1), (B, m1)], kb_2)
        if DECAY:
            m.add_reaction_rule((C, m1), [], kf_1)

    if REACTIONS and MEMBRANE2:
        pass

    if REACTIONS and DNA:
        m.add_reaction_rule((A, d),           (B, d),           kf_1)
        m.add_reaction_rule((B, d),           (A, d),           kb_1)
        m.add_reaction_rule([(A, d), (B, d)], (C, d),           kf_2)
        m.add_reaction_rule((C, d),           [(A, d), (B, d)], kb_2)
        if DECAY:
            m.add_reaction_rule((C, d), [], kf_1)


    '''Add surface binding interaction rules.

    '''
    if SURFACE_BINDING_INTERACTIONS and MEMBRANE1 and BOX:
        # ParticleType C can bind from the box to membrane1. The 
        # membrane is reflective, by default, to ParticleTypes A and B.
        m.add_reaction_rule((C, b), (C, m1), kon)

    if SURFACE_BINDING_INTERACTIONS and MEMBRANE2 and BOX:
        # Membrane 2 absorbs all particles.
        m.add_reaction_rule((A, b), (0, m2), kon)
        m.add_reaction_rule((B, b), (0, m2), kon)
        m.add_reaction_rule((C, b), (0, m2), kon)

    if SURFACE_BINDING_INTERACTIONS and DNA and BOX:
        # ParticleType C can bind from the box to the dna. The dna is 
        # reflective, by default, to ParticleTypes A and B.
        m.add_reaction_rule((C, b), (C, d), kon)


    ''' Add surface unbinding reaction rules.

    '''
    if SURFACE_UNBINDING_REACTIONS and MEMBRANE1 and BOX:
        # Species C can unbind from membrane1 to the box.
        m.add_reaction_rule((C, m1), (C, b), koff)

    if SURFACE_UNBINDING_REACTIONS and MEMBRANE2 and BOX:
        pass

    if SURFACE_UNBINDING_REACTIONS and DNA and BOX:
        # Species C can unbind from the dna to the box.
        m.add_reaction_rule((C, d), (C, b), koff)


else:
    '''Add Surfaces and Regions. 

    '''
    if BOX:
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
        A = Species('A', D, radius, m.get_structure("world").id)
        B = Species('B', D, radius, m.get_structure("world").id)
        C = Species('C', D, radius, m.get_structure("world").id)
        m.add_species_type(A)
        m.add_species_type(B)
        m.add_species_type(C)

    if BOX:
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
        Ad = Species('Ad', D, radius, d.id)
        Bd = Species('Bd', D, radius, d.id)
        Cd = Species('Cd', D, radius, d.id)
        m.add_species_type(Ad)
        m.add_species_type(Bd)
        m.add_species_type(Cd)


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

    if REACTIONS and BOX:
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


    '''Add surface binding interaction rules.

    TODO
    '''
    pass


    ''' Add surface unbinding reaction rules.

    TODO
    '''
    pass


'''Define simulator.

A cube of dimension (L x L x L) is created, with periodic boundary 
conditions.

'''
matrix_size = 3
w = create_world(m, matrix_size)
nrw = NetworkRulesWrapper(m.network_rules)
s = EGFRDSimulator(w, myrandom.rng, nrw)
s.bd_dt_factor = BD_DT_FACTOR


if EGFRD_MODEL == True:
    '''Add particles.

    '''
    if MEMBRANE1 and MEMBRANE2 and BOX:
        # Add world particles inside box.
        throw_in(s, (A, b), N_PARTICLES_FACTOR * 4)
        throw_in(s, (B, b), N_PARTICLES_FACTOR * 4)
        throw_in(s, (C, b), N_PARTICLES_FACTOR * 16)
    elif WORLD:
        throw_in(s, A, N_PARTICLES_FACTOR * 8)
        throw_in(s, A2, N_PARTICLES_FACTOR * 2)
        throw_in(s, B, N_PARTICLES_FACTOR * 4)
        throw_in(s, B2, N_PARTICLES_FACTOR * 8)
        throw_in(s, C, N_PARTICLES_FACTOR * 2)
        throw_in(s, C2, N_PARTICLES_FACTOR * 4)
        throw_in(s, D2, N_PARTICLES_FACTOR * 8)

    if MEMBRANE1:
        throw_in(s, (A, m1), N_PARTICLES_FACTOR * 2)
        throw_in(s, (B, m1), N_PARTICLES_FACTOR * 2)
        throw_in(s, (C, m1), N_PARTICLES_FACTOR * 8)

    if MEMBRANE2:
        pass

    if DNA:
        throw_in(s, (A, d), N_PARTICLES_FACTOR * 1)
        throw_in(s, (B, d), N_PARTICLES_FACTOR * 1)
        throw_in(s, (C, d), N_PARTICLES_FACTOR * 4)


else:
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
    from vtklogger import VTKLogger
    # Write vtk files. See vtklogger.py. Use VTK or Paraview to 
    # visualize.  
    vtklogger = VTKLogger(s, 'data/run', extra_particle_step=True, 
                          buffer_size=20)

if LOGGER == True:
    from logger import Logger
    # Write log files. See logger.py. Use visualizer.py to visualize.
    l = Logger('example')
    l.start(s)
    #l.set_particle_out_interval(1e-4) #1e-3)


'''Run  simulation.

'''
def run():
    print s.dump_reaction_rules()

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

