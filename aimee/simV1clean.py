#!/usr/bin/env python


from egfrd import *
from model2 import *

from logger import *
#from vtklogger import VTKLogger


#import sys
#import math

'''

See: "rasSimsNotesV1.pdf" for more details. 

The reactants/products (all in the membrane):

RasD   Ras-GDP
RasT   Ras-GTP
Sos    Sos, no allosteric binding
SosD   Sos, Ras-GDP bound to the allosteric site (Sos-Ras-GDP)
SosT   Sos, Ras-GTP bound to the allosteric site (Sos-Ras-GTP)
SosDD  Sos-Ras-GDP, with Ras-GDP bound
SosTD  Sos-Ras-GTP, with Ras-GTP bound


The reaction network:

1. Membrane binding/unbinding: 

None for now.

2. Membrane reactions:

1. Sos + RasD    <->  SosD          (k1/km1)
2. Sos + RasT    <->  SosT          (k2/km2)
3. SosD + RasD   <->  SosDD         (k3/km3)
4. SosDD          ->  SosD + RasT   (k3cat)
5. SosT + RasD   <->  SosTD         (k4/km4)
6. SosTD          ->  SosT + RasT   (k4cat)
7. RasT + GAP    <->  RasTGAP       (k5)
8. RasTGAP        ->  RasD + GAP    (k5cat)
 

3. Purely cytosolic reactions: 

None for now. 


'''

def run( ):

    L  = 2e-6            # X-size of simulation box.
    D  = 1e-12           # Diffusion constant. (this was measured in experiment)  kon must not be larger than kD for the k_a() and k_d() functions. 
    radius = 1.7e-9      # Radius of particles.


    m = ParticleModel(L)

    WORLD = True
    MEMBRANE = True
    MEMBRANE2 = False

    if MEMBRANE:

        m1 = m.add_planar_surface(id='m1',
                                  origin=[L/2, L/2, L/2 ],
                                  unit_x=[1, 0, 0],
                                  unit_y=[0, 1, 0],
                                  Lx=(L / 2),
                                  Ly=(L / 2))
#                                 Lz=(L / 2),
#                                 Lz=(radius),


    RasD    = Species( 'RasD',    D, radius )
    RasT    = Species( 'RasT',    D, radius )
    Sos     = Species( 'Sos',     D, radius )
    SosD    = Species( 'SosD',    D, radius )
    SosT    = Species( 'SosT',    D, radius )
    SosDD   = Species( 'SosDD',   D, radius )
    SosDT   = Species( 'SosDT',   D, radius )
    SosTD   = Species( 'SosTD',   D, radius )
    SosTT   = Species( 'SosTT',   D, radius )
    GAP     = Species( 'GAP',     D, radius )
    RasTGAP = Species( 'RasTGAP', D, radius )


       
    if MEMBRANE:
        m.add_species(RasD,    m1)
        m.add_species(RasT,    m1)
        m.add_species(Sos,     m1)
        m.add_species(SosD,    m1)
        m.add_species(SosT,    m1)
        m.add_species(SosDD,   m1)
        m.add_species(SosDT,   m1)
        m.add_species(SosTD,   m1)
        m.add_species(SosTT,   m1)
        m.add_species(GAP,     m1)
        m.add_species(RasTGAP, m1)
        


    
    sigma = 2.0*radius
#    kD    = k_D( 2.0*D, sigma) # This is for 3D
    kD    = k_D( D, sigma) # This is for 2D

    k1    = 2.0e-22
    km1   = 3.0
    k2    = 1.8333e-22
    km2   = 0.4
    k3    = 8.3333e-23
    km3   = 0.1
    k3cat = 0.038  
    k4    = 1.1667e-22
    km4   = 1.0
    k4cat = 0.003
    k5    = 2.9000e-21
    km5   = 0.2    
    k5cat = 0.1

#    k1    = 2.0e-10
#    k2    = 1.833e-10
#    k3    = 8.33e-11
#    k4    = 1.166e-10
#    k5    = 1.74e-10   I don't know how I ended up calculating these.  This number appears to differ in the value, not only 12 orders of magnitude. 


#    '''
    k1a   = k_a( k1, kD )
    k1d   = k_d( km1, k1, kD )
    k2a   = k_a( k2, kD )
    k2d   = k_d( km2, k2, kD )
    k3a   = k_a( k3, kD )
    k3d   = k_d( km3, k3, kD ) 
    k4a   = k_a( k4, kD ) 
    k4d   = k_d( km4, k4, kD )  
    k5a   = k_a( k5, kD )  
    k5d   = k_d( km5, k5, kD )   

    if MEMBRANE:
        
        m.add_reaction([(Sos, m1), (RasD, m1) ], [(SosD, m1)            ], k1a  )
#       m.add_reaction([(SosD, m1)            ], [(Sos, m1), (RasT, m1)  ], 
#       k1d          )
        m.add_reaction([(SosD, m1)            ], [(Sos, m1), (RasD, m1) ], k1d  )
        m.add_reaction([(Sos, m1), (RasT, m1) ], [(SosT, m1)            ], k2a  )
        m.add_reaction([(SosT, m1)            ], [(Sos, m1), (RasT, m1) ], k2d  )
        m.add_reaction([(SosT, m1), (RasD, m1)], [(SosTD, m1)           ], k3a  )
        m.add_reaction([(SosTD, m1)           ], [(SosT, m1), (RasD, m1)], k3d  )
        m.add_reaction([(SosTD, m1)           ], [(SosT, m1), (RasT, m1)], k3cat)
        m.add_reaction([(SosD, m1), (RasD, m1)], [(SosDD, m1)           ], k4a  )
        m.add_reaction([(SosDD, m1)           ], [(SosD, m1), (RasD, m1)], k4d  )
        m.add_reaction([(SosDD, m1)           ], [(SosD, m1), (RasT, m1)], k4cat)
        m.add_reaction([(RasT, m1), (GAP, m1) ], [(RasTGAP, m1)         ], k5a  )
        m.add_reaction([(RasTGAP, m1)         ], [(RasT, m1), (GAP, m1) ], k5d  )
        m.add_reaction([(RasTGAP, m1)         ], [(RasD, m1), (GAP, m1) ], k5cat)




    '''        
    if MEMBRANE:
        s.addReaction( [ (Sos, m), (RasD, m)   ], [ (SosD, m)            ], k1    )
        s.addReaction( [ (SosD, m)             ], [ (Sos, m), (RasT, m)  ], km1   )
        s.addReaction( [ (Sos, m), (RasT, m)   ], [ (SosT, m)            ], k2    )
        s.addReaction( [ (SosT, m)             ], [ (Sos, m), (RasT, m)  ], km2   )
        s.addReaction( [ (SosT, m), (RasD, m)  ], [ (SosTD, m)           ], k3    )
        s.addReaction( [ (SosTD, m)            ], [ (SosT, m), (RasD, m) ], km3   )
        s.addReaction( [ (SosTD, m)            ], [ (SosT, m), (RasT, m) ], k3cat )
        s.addReaction( [ (SosD, m), (RasD, m)  ], [ (SosDD, m)           ], k4    )
        s.addReaction( [ (SosDD, m)            ], [ (SosD, m), (RasD, m) ], km4   )
        s.addReaction( [ (SosDD, m)            ], [ (SosD, m), (RasT, m) ], k4cat )
        s.addReaction( [ (RasT, m), (GAP, m)   ], [ (RasTGAP, m)         ], k5    )
        s.addReaction( [ (RasTGAP, m)          ], [ (RasT, m), (GAP, m)  ], km5   )
        s.addReaction( [ (RasTGAP, m)          ], [ (RasD, m), (GAP, m)  ], k5cat )



    k1 = 1e-18
    k2 = 10
    if MEMBRANE:
        s.addReaction( [ (Sos, m), (RasD, m)   ], [ (SosD, m)            ], k2   )
        s.addReaction( [ (SosD, m)             ], [ (Sos, m), (RasT, m)  ], k1   )
        s.addReaction( [ (Sos, m), (RasT, m)   ], [ (SosT, m)            ], k2   )
        s.addReaction( [ (SosT, m)             ], [ (Sos, m), (RasT, m)  ], k2   )
        s.addReaction( [ (SosT, m), (RasD, m)  ], [ (SosTD, m)           ], k2   )
        s.addReaction( [ (SosTD, m)            ], [ (SosT, m), (RasD, m) ], k1   )
        s.addReaction( [ (SosTD, m)            ], [ (SosT, m), (RasT, m) ], k1   )
        s.addReaction( [ (SosD, m), (RasD, m)  ], [ (SosDD, m)           ], k2   )
        s.addReaction( [ (SosDD, m)            ], [ (SosD, m), (RasD, m) ], k1   )
        s.addReaction( [ (SosDD, m)            ], [ (SosD, m), (RasT, m) ], k1   )
        s.addReaction( [ (RasT, m), (GAP, m)   ], [ (RasTGAP, m)         ], k2   )
        s.addReaction( [ (RasTGAP, m)          ], [ (RasT, m), (GAP, m)  ], k1   )
        s.addReaction( [ (RasTGAP, m)          ], [ (RasD, m), (GAP, m)  ], k1   )

    '''


    matrix_size = 3
    w = create_world(m, matrix_size)
    nrw = NetworkRulesWrapper(m.network_rules)
    s = EGFRDSimulator(w, myrandom.rng, nrw)


    if MEMBRANE: 
        throw_in(s.world, Sos,  100, m1)  
        throw_in(s.world, RasD,  50, m1)
        throw_in(s.world, RasT,  50, m1)
        #throw_in(s.world, GAP,  100, m1) 




    s.initialize()
    print s.dump_reaction_rules()






    '''
    vtklogger = VTKLogger( s, 'run', 100 )
    for i in range( 1000 ):
        try:
            vtklogger.log()
            s.step()
        except RuntimeError, message:
            print message
            break

    vtklogger.stop()
    s.stop( s.t )
    '''



    
    endTime = 2.0
  
    l = Logger( s ) 
#    l.setInterval( 1e-6 )
    rfile = open( 'data/reactions.dat', 'w')

    l.log() # this writes to data/log_tc.dat

    
    while s.t < endTime: 
        s.step()

    
        if s.lastReaction: 
            r = s.lastReaction
            line = '(%18.18g, \t%s, \t%s ) \n' % (s.t, r.reactants, r.products ) 
            rfile.write( line )
            rfile.flush()

            l.log()
    
    


# End "def run( ):"



if __name__ == '__main__':
    run( )
    
    


