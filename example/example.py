#!/usr/bin/env python

'''
Running this script to test if everything works:
PYTHONPATH=../../ python alternative-user-interface.py

or, without terminal output and faster:
PYTHONPATH=../../ python -O alternative-user-interface.py
'''

from egfrd import *
from dumper import *
import myrandom


'''Reaction network:

A     <-> B
A + B <-> C
C      -> 0
'''


'''Settings.

'''
ALTERNATIVE_USER_INTERFACE = True

# Make multis run 'BD_DT_FACTOR' times faster than normal.

BD_DT_FACTOR = 1
RADIUS_FACTOR = 1
SINGLE_RATE_FACTOR = 1
PAIR_RATE_FACTOR = 1

WORLD = True
MEMBRANE1 = True
MEMBRANE2 = False
DNA = False

DECAY = True

VTK_LOGGER = True
LOGGER = False

MY_SEED = 0
myrandom.seed(MY_SEED)


def run():
    '''Dimensions.

    '''
    L = 1e-6            # Size of simulation box.
    D = 1e-12           # Diffusion constant.
    radius = RADIUS_FACTOR * 2.5e-9 # Radius of particles and surfaces.
    thickness = radius


    '''Simulator.

    A cube of dimension (L x L x L) is created, with periodic boundary 
    conditions.

    '''
    matrixSize = 3
    w = World(L, matrixSize)
    s = EGFRDSimulator(w)
    s.bd_dt_factor = BD_DT_FACTOR


    '''Create model

    '''
    model = ParticleModel()


    '''Add surfaces. 

    Define the surfaces to add to the model.

    Alternative user interface:
        model.addPlanarSurface(name, origin, vectorX, vectorY, Lx, Ly, Lz)
    , where origin is the middle of the plane, and Lz is the thickness.
        model.addCylindricalSurface(name, origin, radius, orientation, size)
    Also see the docstrings. 

    Note that surface are not allowed to touch or overlap. Surfaces should 
    extend over the whole simulation box, so finite surfaces are not 
    supported.

    Both methods return the surface they create. Assign it 
    to some variable, so it can be used when adding species.

    '''
    if MEMBRANE1:
        m1 = model.addPlanarSurface(origin=[L/2, L/2, 2*L/10],
                                    vectorX=[1, 0, 0],
                                    vectorY=[0, 1, 0],
                                    Lx=L/2,
                                    Ly=L/2,
                                    Lz=thickness,
                                    name='m1')

    if MEMBRANE2:
        m2 = model.addPlanarSurface(origin=[L/2, L/2, 8*L/10],
                                    vectorX=[1, 0, 0], 
                                    vectorY=[0, 1, 0],
                                    Lx=L/2,
                                    Ly=L/2,
                                    Lz=thickness,
                                    name='m2')

    if DNA:
        d = model.addCylindricalSurface(origin=[L/2, L/2, L/2],
                                        radius=radius,
                                        orientation=[0, 1, 0],
                                        size=L/2,
                                        name='d')


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
        model.addSpecies(species)

    When adding a species to the model without explicitly assigning the 
    surface it can live on, it is assumed to be in the 'world'. The 'world' is 
    a 3D space, also referred to as the 'bulk' or 'cytoplasm'.

    '''
    if WORLD:
        if ALTERNATIVE_USER_INTERFACE == True:
            model.addSpecies(A)
            model.addSpecies(B)
            model.addSpecies(C)
        else:
            A = model.new_species_type('A', D, radius)
            A['surface'] = model.defaultSurface.name
            B = model.new_species_type('B', D, radius)
            B['surface'] = model.defaultSurface.name
            C = model.new_species_type('C', D, radius)
            C['surface'] = model.defaultSurface.name


    '''Add species to surfaces.

    Specify which species can exist on which surface.

    Alternative user interface:
        model.addSpecies(species, surface, D, radius)
    See the docstring.

    Per surface a different diffusion constant D and radius can be specified. 
    By default the ones for the 'world' are used. 

    Note: if particles of species 'A' should be able to exist in the world as 
    well as on one of the previously added surfaces, then it should be added 
    twice. Ones with and ones without an explicit surface argument.

    '''
    if MEMBRANE1:
        if ALTERNATIVE_USER_INTERFACE == True:
            model.addSpecies(A, m1)
            model.addSpecies(B, m1)
            model.addSpecies(C, m1)
        else:
            Am1 = model.new_species_type('Am1', D, radius)
            Am1['surface'] = m1.name
            Bm1 = model.new_species_type('Bm1', D, radius)
            Bm1['surface'] = m1.name
            Cm1 = model.new_species_type('Cm1', D, radius)
            Cm1['surface'] = m1.name

    if MEMBRANE2:
        # No species can live on membrane 2.
        pass

    if DNA:
        if ALTERNATIVE_USER_INTERFACE == True:
            model.addSpecies(A, d)
            model.addSpecies(B, d)
            model.addSpecies(C, d)
        else:
            Ad = model.new_species_type('Ad', D, radius)
            Ad['surface'] = d.name
            Bd = model.new_species_type('Bd', D, radius)
            Bd['surface'] = d.name
            Cd = model.new_species_type('Cd', D, radius)
            Cd['surface'] = d.name


    '''Add reactions in 3D.

    Alternative user interface:
        model.addReaction([reactants], [products], rate)

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
            model.addReaction([A],    [B],    kf_1)
            model.addReaction([B],    [A],    kb_1)
            model.addReaction([A, B], [C],    kf_2)
            model.addReaction([C],    [A, B], kb_2)
            if DECAY:
                model.addReaction([C], [], kf_1)
        else:
            r1 = createUnimolecularReactionRule(A, B, kf_1)
            model.network_rules.add_reaction_rule(r1)

            r2 = createUnimolecularReactionRule(B, A, kb_1)
            model.network_rules.add_reaction_rule(r2)

            r3 = createBindingReactionRule(A, B, C, kf_2)
            model.network_rules.add_reaction_rule(r3)

            r4 = createUnbindingReactionRule(C, A, B, kb_2)
            model.network_rules.add_reaction_rule(r4)

            if DECAY:
                r5 = createDecayReactionRule(C, kf_1)
                model.network_rules.add_reaction_rule(r5)


    '''Add reactions on surfaces.

    Alternative user interface:
        model.addReaction([reactants], [products], rate)
    , where one reactant or product is a tuple: (species, surface).

    Combinations of species which do not appear together as reactants in any 
    reaction are made repulsive by default on every surface.

    '''
    if MEMBRANE1:
        if ALTERNATIVE_USER_INTERFACE == True:
            model.addReaction([(A, m1)],          [(B, m1)],          kf_1)
            model.addReaction([(B, m1)],          [(A, m1)],          kb_1)
            model.addReaction([(A, m1), (B, m1)], [(C, m1)],          kf_2)
            model.addReaction([(C, m1)],          [(A, m1), (B, m1)], kb_2)
            if DECAY:
                model.addReaction([(C, m1)], [], kf_1)
        else:
            r1 = createUnimolecularReactionRule(Am1, Bm1, kf_1)
            model.network_rules.add_reaction_rule(r1)

            r2 = createUnimolecularReactionRule(Bm1, Am1, kb_1)
            model.network_rules.add_reaction_rule(r2)

            r3 = createBindingReactionRule(Am1, Bm1, Cm1, kf_2)
            model.network_rules.add_reaction_rule(r3)

            r4 = createUnbindingReactionRule(Cm1, Am1, Bm1, kb_2)
            model.network_rules.add_reaction_rule(r4)

            if DECAY:
                r5 = createDecayReactionRule(Cm1, kf_1)
                model.network_rules.add_reaction_rule(r5)

    if MEMBRANE2:
        # No particles can live on membrane2.
        pass

    if DNA:
        if ALTERNATIVE_USER_INTERFACE == True:
            model.addReaction([(A, d)],         [(B, d)],         kf_1)
            model.addReaction([(B, d)],         [(A, d)],         kb_1)
            model.addReaction([(A, d), (B, d)], [(C, d)],         kf_2)
            model.addReaction([(C, d)],         [(A, d), (B, d)], kb_2)
            if DECAY:
                model.addReaction([(C, d)], [], kf_1)
        else:
            r1 = createUnimolecularReactionRule(Ad, Bd, kf_1)
            model.network_rules.add_reaction_rule(r1)

            r2 = createUnimolecularReactionRule(Bd, Ad, kb_1)
            model.network_rules.add_reaction_rule(r2)

            r3 = createBindingReactionRule(Ad, Bd, Cd, kf_2)
            model.network_rules.add_reaction_rule(r3)

            r4 = createUnbindingReactionRule(Cd, Ad, Bd, kb_2)
            model.network_rules.add_reaction_rule(r4)

            if DECAY:
                r5 = createDecayReactionRule(Cd, kf_1)
                model.network_rules.add_reaction_rule(r5)


    ''' Surface binding reactions.

    Alternative user interface:
        model.addReaction([reactant], [product], rate))
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
            model.addReaction([C], [(C, m1)], kon)
        else:
            pass

    if MEMBRANE2 and WORLD:
        if ALTERNATIVE_USER_INTERFACE == True:
            # Membrane 2 absorbs all particles.
            model.addReaction([A], [(0, m2)], kon)
            model.addReaction([B], [(0, m2)], kon)
            model.addReaction([C], [(0, m2)], kon)
            pass
        else:
            pass

    if DNA and WORLD:
        if ALTERNATIVE_USER_INTERFACE == True:
            # Species C can bind to the dna. The dna is reflective, by 
            # default, to species A and B.
            model.addReaction([C], [(C, d)], kon)
        else:
            pass


    ''' Surface unbinding reactions.

    Alternative user interface:
        model.addReaction([reactant], [product], k))
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
            model.addReaction([(C, m1)], [C], koff)
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
            model.addReaction([(C, d)], [C], koff)
        else:
            pass



    '''Set model.

    '''
    model.set_all_repulsive()
    s.setModel(model)


    '''Add particles.

    '''
    if WORLD:
        if MEMBRANE1 and MEMBRANE2:
            # Add world particles inside the two planes.
            # Note that a CuboidalRegion is defined by 2 corners.
            box1 = CuboidalRegion([0, 0, 2 * L / 10], [L, L, 8 * L / 10])
            s.throwInParticles(C, 4, box1)
            s.throwInParticles(C, 4, box1)
            s.throwInParticles(C, 16, box1)
        else:
            # Particles are added to world (3D) by default.
            s.throwInParticles(A, 4)
            s.throwInParticles(B, 4)
            s.throwInParticles(C, 16)

    if MEMBRANE1:
        if ALTERNATIVE_USER_INTERFACE == True:
            s.throwInParticles(A, 2, m1)
            s.throwInParticles(B, 2, m1)
            s.throwInParticles(C, 8, m1)
        else:
            s.throwInParticles(Am1, 2)
            s.throwInParticles(Bm1, 2)
            s.throwInParticles(Cm1, 8)

    if MEMBRANE2: pass

    if DNA:
        if ALTERNATIVE_USER_INTERFACE == True:
            s.throwInParticles(A, 1, d)
            s.throwInParticles(B, 1, d)
            s.throwInParticles(C, 4, d)
        else:
            s.throwInParticles(Ad, 1)
            s.throwInParticles(Bd, 1)
            s.throwInParticles(Cd, 4)


    '''Define loggers.

    '''
    if VTK_LOGGER == True:
        from vtklogger import VTKLogger
        # Write vtk files. See vtklogger.py. Use VTK or Paraview to visualize.  
        vtklogger = VTKLogger(s, 'run', 100)

    if LOGGER == True:
        from logger import Logger
        # Write log files. See logger.py. Use visualizer.py to visualize.
        l = Logger(s, 'example')
        l.setParticleOutInterval(1e-4) #1e-3)

    print s.dump_reaction_rules()


    '''Run  simulation.

    '''
    s.initialize()
    for i in range(100):
	print i
        try:
            if VTK_LOGGER == True:
                vtklogger.log()
            if LOGGER == True:
                l.log()
            s.step()
        except Exception, message:
            print 'Exception at step = ', i
            if VTK_LOGGER == True:
                vtklogger.stop()
            s.print_report()
            raise

    if VTK_LOGGER == True:
        vtklogger.stop()

    '''Todo: add test for this stuff.
    s.print_report()
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

