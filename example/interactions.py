
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
    B = ParticleType('B', D, radius)
    C = ParticleType('C', D, 2 * radius)


    '''Add ParticleTypes (to specific Regions or Surfaces).

    '''
    if WORLD:
        m.add_particle_type(A)
        m.add_particle_type(B)
        m.add_particle_type(C)

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
if EGFRD_MODEL == True:
    '''Add particles.

    '''
    if MEMBRANE1 and MEMBRANE2 and BOX:
        # Add world particles inside box.
        throw_in(s, (A, b), N_PARTICLES_FACTOR * 4)
        throw_in(s, (B, b), N_PARTICLES_FACTOR * 4)
        throw_in(s, (C, b), N_PARTICLES_FACTOR * 16)
    elif WORLD:
        throw_in(s, A, N_PARTICLES_FACTOR * 4)
        throw_in(s, B, N_PARTICLES_FACTOR * 4)
        throw_in(s, C, N_PARTICLES_FACTOR * 16)

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
