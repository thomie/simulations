#!/usr/bin/env python

'''
To run this script:
PYTHONPATH=../../egfrd python example.py

or, without terminal output and faster:
PYTHONPATH=../../egfrd python -O example.py
'''

from egfrd import *
from dumper import *
import myrandom


'''Reaction network:

A     <-> B
A + B <-> C
C      -> 0
'''


WORLD                      = None
MEMBRANE1                  = None
MEMBRANE2                  = None
DNA                        = None
DECAY                      = None
VTK_LOGGER                 = None
LOGGER                     = None
ALTERNATIVE_USER_INTERFACE = None


'''Settings.

'''
WORLD                      = True
MEMBRANE1                  = True
#MEMBRANE2                  = True # Todo. But first fix interactions.
DNA                        = True
DECAY                      = True
VTK_LOGGER                 = True
#LOGGER                     = True # Todo.
ALTERNATIVE_USER_INTERFACE = True

if ALTERNATIVE_USER_INTERFACE == True:
    from model2 import *
else:
    from model import *

# Make multis run 'BD_DT_FACTOR' times faster than normal.
BD_DT_FACTOR = 1
RADIUS_FACTOR = 1
SINGLE_RATE_FACTOR = 1
PAIR_RATE_FACTOR = 1
N_PARTICLES_FACTOR = 0.6

MY_SEED = 0
myrandom.seed(MY_SEED)


def run():
    '''Dimensions.

    '''
    L = 1e-6                        # Size of simulation box.
    D = 1e-12                       # Diffusion constant.
    radius = RADIUS_FACTOR * 2.5e-9 # Radius of particles and surfaces.
    thickness = radius


    '''Create model

    '''
    m = ParticleModel(L)


    '''Add surfaces. 

    Define the surfaces to add to the model.

    Alternative user interface:
        m.addPlanarSurface(name, origin, vectorX, vectorY, Lx, Ly, Lz)
    , where origin is the middle of the plane, and Lz is the thickness.
        m.addCylindricalSurface(name, origin, radius, orientation, size)
    Also see the docstrings. 

    Note that surface are not allowed to touch or overlap. Surfaces should 
    extend over the whole simulation box, so finite surfaces are not 
    supported.

    Both methods return the surface they create. Assign it 
    to some variable, so it can be used when adding species.

    '''
    if MEMBRANE1:
        m1 = m.add_planar_surface(id='m1',
                                  origin=[L/2, L/2, 2*L/10],
                                  unit_x=[1, 0, 0],
                                  unit_y=[0, 1, 0],
                                  Lx=L/2,
                                  Ly=L/2)
                                  #Lz=thickness)

    if MEMBRANE2:
        m2 = m.add_planar_surface(id='m2',
                                  origin=[L/2, L/2, 8*L/10],
                                  unit_x=[1, 0, 0],
                                  unit_y=[0, 1, 0],
                                  Lx=L/2,
                                  Ly=L/2)
                                  #Lz=thickness)

    if DNA:
        d = m.add_cylindrical_surface(id='d',
                                      origin=[L/2, L/2, L/2],
                                      radius=radius,
                                      orientation=[0, 1, 0],
                                      size=L/2)


    '''Define species.

    Define the species that can exist.

    Alternative user interface:
        species = Species('name', D, radius)

    If no D or radius is specified, it has to be specified when adding it to a 
    specific surface.

    '''
    if ALTERNATIVE_USER_INTERFACE == True:
        A = Species('A', D, radius)
        B = Species('B', D, radius)
        C = Species('C', D, radius)


    '''Add species to 3D.

    Alternative user interface:
        m.add_species(species)

    When adding a species to the model without explicitly assigning the 
    surface it can live on, it is assumed to be in the 'world'. The 'world' is 
    a 3D space, also referred to as the 'bulk' or 'cytoplasm'.

    '''
    if WORLD:
        if ALTERNATIVE_USER_INTERFACE == True:
            m.add_species(A)
            m.add_species(B)
            m.add_species(C)
        else:
            A = Species('A', D, radius, m.get_structure("world").id)
            B = Species('B', D, radius, m.get_structure("world").id)
            C = Species('C', D, radius, m.get_structure("world").id)
            m.add_species_type(A)
            m.add_species_type(B)
            m.add_species_type(C)


    '''Add species to surfaces.

    Specify which species can exist on which surface.

    Alternative user interface:
        m.add_species(species, surface, D, radius)
    See the docstring.

    Per surface a different diffusion constant D and radius can be specified. 
    By default the ones for the 'world' are used. 

    Note: if particles of species 'A' should be able to exist in the world as 
    well as on one of the previously added surfaces, then it should be added 
    twice. Ones with and ones without an explicit surface argument.

    '''
    if MEMBRANE1:
        if ALTERNATIVE_USER_INTERFACE == True:
            m.add_species(A, m1)
            m.add_species(B, m1)
            m.add_species(C, m1)
        else:
            Am1 = Species('Am1', D, radius, m1.id)
            Bm1 = Species('Bm1', D, radius, m1.id)
            Cm1 = Species('Cm1', D, radius, m1.id)
            m.add_species_type(Am1)
            m.add_species_type(Bm1)
            m.add_species_type(Cm1)

    if MEMBRANE2:
        # No species can live on membrane 2.
        pass

    if DNA:
        if ALTERNATIVE_USER_INTERFACE == True:
            m.add_species(A, d)
            m.add_species(B, d)
            m.add_species(C, d)
        else:
            Ad = Species('Ad', D, radius, d.id)
            Bd = Species('Bd', D, radius, d.id)
            Cd = Species('Cd', D, radius, d.id)
            m.add_species_type(Ad)
            m.add_species_type(Bd)
            m.add_species_type(Cd)


    '''Add reactions in 3D.

    Alternative user interface:
        m.add_reaction([reactants], [products], rate)

    For now: a bimolecular reaction can only have 1 product species.

    '''
    sigma = 2 * radius
    D_tot = D
    tau = sigma * sigma / D_tot
    # Bimolecular. 
    kf_2 = PAIR_RATE_FACTOR * 100 * sigma * D_tot
    # Unimolecular.
    kb_2 = PAIR_RATE_FACTOR * 0.1 / tau
    kf_1 = SINGLE_RATE_FACTOR * 0.1 / tau
    kb_1 = SINGLE_RATE_FACTOR * 0.1 / tau

    if WORLD:
        if ALTERNATIVE_USER_INTERFACE == True:
            m.add_reaction([A],    [B],    kf_1)
            m.add_reaction([B],    [A],    kb_1)
            m.add_reaction([A, B], [C],    kf_2)
            m.add_reaction([C],    [A, B], kb_2)
            if DECAY:
                m.add_reaction([C], [], kf_1)
        else:
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


    '''Add reactions on surfaces.

    Alternative user interface:
        m.add_reaction([reactants], [products], rate)
    , where one reactant or product is a tuple: (species, surface).

    Combinations of species which do not appear together as reactants in any 
    reaction are made repulsive by default on every surface.

    '''
    if MEMBRANE1:
        if ALTERNATIVE_USER_INTERFACE == True:
            m.add_reaction([(A, m1)],          [(B, m1)],          kf_1)
            m.add_reaction([(B, m1)],          [(A, m1)],          kb_1)
            m.add_reaction([(A, m1), (B, m1)], [(C, m1)],          kf_2)
            m.add_reaction([(C, m1)],          [(A, m1), (B, m1)], kb_2)
            if DECAY:
                m.add_reaction([(C, m1)], [], kf_1)
        else:
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

    if MEMBRANE2:
        # No particles can live on membrane2.
        pass

    if DNA:
        if ALTERNATIVE_USER_INTERFACE == True:
            m.add_reaction([(A, d)],         [(B, d)],         kf_1)
            m.add_reaction([(B, d)],         [(A, d)],         kb_1)
            m.add_reaction([(A, d), (B, d)], [(C, d)],         kf_2)
            m.add_reaction([(C, d)],         [(A, d), (B, d)], kb_2)
            if DECAY:
                m.add_reaction([(C, d)], [], kf_1)
        else:
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


    ''' Surface binding reactions.

    Alternative user interface:
        m.add_reaction([reactant], [product], rate))
    , where product is a tuple: (species, surface).

    The reactant species for every surface binding reaction is always a 
    species that can only exist in the 'world', so no surface has to be 
    specified there.

    If a surface should be absorbive, specify (0, surface) as the product.

    When no surface binding reaction is defined for a combination of a species 
    and a surface, they are made repulsive by default.

    '''
    # Molecule + surface.
    kon = kf_2

    if MEMBRANE1 and WORLD:
        if ALTERNATIVE_USER_INTERFACE == True:
            # Species C can bind to the membrane. The membrane is reflective, 
            # by default, to species A and B.
            m.add_reaction([C], [(C, m1)], kon)
        else:
            pass

    if MEMBRANE2 and WORLD:
        if ALTERNATIVE_USER_INTERFACE == True:
            # Membrane 2 absorbs all particles.
            m.add_reaction([A], [(0, m2)], kon)
            m.add_reaction([B], [(0, m2)], kon)
            m.add_reaction([C], [(0, m2)], kon)
            pass
        else:
            pass

    if DNA and WORLD:
        if ALTERNATIVE_USER_INTERFACE == True:
            # Species C can bind to the dna. The dna is reflective, by 
            # default, to species A and B.
            m.add_reaction([C], [(C, d)], kon)
        else:
            pass


    ''' Surface unbinding reactions.

    Alternative user interface:
        m.add_reaction([reactant], [product], k))
    , where reactant is a tuple: (species, surface).

    Unbinding is a Poissonian reaction, from a surface to the 'world'.

    The product species for every surface binding reaction is always a 
    species that can only exist in the 'world', so no surface has to be 
    specified there.

    '''
    # Unimolecular.
    koff = kb_1

    if MEMBRANE1 and WORLD:
        if ALTERNATIVE_USER_INTERFACE == True:
            # Species C can unbind to the 'world'.
            m.add_reaction([(C, m1)], [C], koff)
        else:
            pass

    if MEMBRANE2 and WORLD:
        if ALTERNATIVE_USER_INTERFACE == True:
            # No particles can live on membrane2.
            pass
        else:
            pass

    if DNA and WORLD:
        if ALTERNATIVE_USER_INTERFACE == True:
            # Species C can unbind to the 'world'.
            m.add_reaction([(C, d)], [C], koff)
        else:
            pass



    '''Simulator.

    A cube of dimension (L x L x L) is created, with periodic boundary 
    conditions.

    '''
    matrix_size = 3
    w = create_world(m, matrix_size)
    nrw = NetworkRulesWrapper(m.network_rules)
    s = EGFRDSimulator(w, myrandom.rng, nrw)
    s.bd_dt_factor = BD_DT_FACTOR


    '''Add particles.

    '''
    if WORLD:
        if ALTERNATIVE_USER_INTERFACE == True:
            throw_in(s.world, A, N_PARTICLES_FACTOR * 4)
            throw_in(s.world, B, N_PARTICLES_FACTOR * 4)
            throw_in(s.world, C, N_PARTICLES_FACTOR * 16)
        else:
            # Particles are added to world (3D) by default.
            throw_in_particles(s.world, A, N_PARTICLES_FACTOR * 4)
            throw_in_particles(s.world, B, N_PARTICLES_FACTOR * 4)
            throw_in_particles(s.world, C, N_PARTICLES_FACTOR * 16)
        '''
        if MEMBRANE1 and MEMBRANE2:
            # Add world particles inside the two planes.
            # Note that a CuboidalRegion is defined by 2 corners.
            box1 = CuboidalRegion([0, 0, 2 * L / 10], [L, L, 8 * L / 10])
            s.throw_in_particles(C, N_PARTICLES_FACTOR * 4, box1)
            s.throw_in_particles(C, N_PARTICLES_FACTOR * 4, box1)
            s.throw_in_particles(C, N_PARTICLES_FACTOR * 16, box1)
        '''

    if MEMBRANE1:
        if ALTERNATIVE_USER_INTERFACE == True:
            throw_in(s.world, A, N_PARTICLES_FACTOR * 2, m1)
            throw_in(s.world, B, N_PARTICLES_FACTOR * 2, m1)
            throw_in(s.world, C, N_PARTICLES_FACTOR * 8, m1)
        else:
            throw_in_particles(s.world, Am1, N_PARTICLES_FACTOR * 2)
            throw_in_particles(s.world, Bm1, N_PARTICLES_FACTOR * 2)
            throw_in_particles(s.world, Cm1, N_PARTICLES_FACTOR * 8)

    if MEMBRANE2: pass

    if DNA:
        if ALTERNATIVE_USER_INTERFACE == True:
            throw_in(w, A, N_PARTICLES_FACTOR * 1, d)
            throw_in(w, B, N_PARTICLES_FACTOR * 1, d)
            throw_in(w, C, N_PARTICLES_FACTOR * 4, d)
        else:
            throw_in_particles(s.world, Ad, N_PARTICLES_FACTOR * 1)
            throw_in_particles(s.world, Bd, N_PARTICLES_FACTOR * 1)
            throw_in_particles(s.world, Cd, N_PARTICLES_FACTOR * 4)


    '''Define loggers.

    '''
    if VTK_LOGGER == True:
        from vtklogger import VTKLogger
        # Write vtk files. See vtklogger.py. Use VTK or Paraview to visualize.  
        vtklogger = VTKLogger(s, 'data/run', 100)

    if LOGGER == True:
        from logger import Logger
        # Write log files. See logger.py. Use visualizer.py to visualize.
        l = Logger(s, 'example')
        #l.set_particle_out_interval(1e-4) #1e-3)

    print s.dump_reaction_rules()


    '''Run  simulation.

    '''
    s.initialize()
    for i in range(100):
        try:
            if VTK_LOGGER == True:
                vtklogger.log()
            if LOGGER == True:
                l.log()
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

