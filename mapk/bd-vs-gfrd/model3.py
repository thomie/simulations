#!/usr/bin/env python

# python -O model3.py gfrd
# python -O model3.py bd

from logger import *
from vtklogger import *
import sys
import time

import math



LOGGER = False
INTERVAL = 1 #6
CPU_SECONDS = 10 #1800

V = 1e-15
D_ratio = 1
D_mode = 'normal'
ti = 0 # Lot's of rebindings.
ti = 1e-6 # Bifocation point.
KK_VS_P = 1.0

simulator = sys.argv[1]



# V in liter, L in meter
L = math.pow( V * 1e-3, 1.0 / 3.0 )
N = 180
matrixSize = min( max( 3, int( (3 * N) ** (1.0/3.0) ) ), 60 )
#print 'matrixSize=', matrixSize


if simulator == 'gfrd':
    from egfrd import *
    w = World(L, matrixSize)
    s = EGFRDSimulator(w)
    dataDirectory = 'gfrd-'
    no_shells = False
elif simulator == 'bd':
    from bd import *
    w = World(L, matrixSize)
    s = BDSimulator(w)
    #s.core.dtFactor = DEFAULT_DT_FACTOR
    dataDirectory = 'bd-'
    no_shells = True

dataDirectory += str(CPU_SECONDS) + '_cpu_s-' + str(KK_VS_P) + '_kk_vs_p'



model='mapk3'

if ti == 0:
    ki = float( 'inf' )
else:
    ki = math.log( 2 ) / ti

D_ref = 1e-12

D_move = D_ref * D_ratio

if D_mode == 'normal':
    D_react = D_move
elif D_mode == 'fixed':
    D_react = D_ref

box1 = CuboidalRegion( [0,0,0],[L,L,L] )

radius = 2.5e-9

m = ParticleModel()

K = m.new_species_type( 'K', D_move, radius )
KK = m.new_species_type( 'KK', D_move, radius )
P = m.new_species_type( 'P', D_move, radius )
Kp = m.new_species_type( 'Kp', D_move, radius )
Kpp = m.new_species_type( 'Kpp', D_move, radius )
K_KK = m.new_species_type( 'K_KK', D_move, radius )
Kp_KK = m.new_species_type( 'Kp_KK', D_move, radius )
Kpp_P = m.new_species_type( 'Kpp_P', D_move, radius )
Kp_P = m.new_species_type( 'Kp_P', D_move, radius )

# inactive forms
KKi = m.new_species_type( 'KKi', D_move, radius )
Pi = m.new_species_type( 'Pi', D_move, radius )

color_dict = {
              2:1, # enzyme KK
              3:1, # enzyme P
              10:1, # enzyme KKi
              11:1, # enzyme KKi

              1:4, # K
              6:4, # intermediate1 K_KK
              9:4, # intermediate1 Kp_P
              4:4, # Kp
              7:4, # intermediate2 Kp_KK
              8:4, # intermediate2 Kpp_P
              5:7, # Kpp
             }


s.setModel(m)

#  1 2   K + KK   <-> K_KK
#  3     K_KK       -> Kp + KKi
#  4 5   Kp + KK  <-> Kp_KK
#  6     Kp_KK      -> Kpp + KKi 
#  7 8   Kpp + P <-> Kpp_P
#  9     Kpp_P     -> Kp + Pi
# 10 11  Kp + P  <-> Kp_P
# 12     Kp_P      -> K + Pi
# 13     KKi     -> KK
# 14     Pi      -> P


sigma = radius * 2
kD = k_D( D_react * 2, sigma )

N_K = C2N( 200e-9, V ) 
N_KK = KK_VS_P * 2 * C2N( 50e-9, V )
N_P = (1 - KK_VS_P) * 2 * C2N( 50e-9, V )


s.throwInParticles( K, N_K, box1 )
print 'Number of K =', N_K
s.throwInParticles( KK, N_KK, box1 )
print 'Number of KK =', N_KK
s.throwInParticles( P, N_P, box1 )
print 'Number of P =', N_P

if simulator == 'gfrd':
    # Initialize scheduler, cut gfrd some slack.
    endTime = 0
    while 1:
        s.step()
        nextTime = s.scheduler.getTopTime()
        if nextTime > endTime:
            s.stop( endTime )
            break
s.reset()

k1 = k_a( Mtom3( 0.02e9 ), kD )
k2 = k_d( 1.0, Mtom3( 0.02e9 ), kD )
k3 = 1.5
k4 = k_a( Mtom3( 0.032e9 ), kD )
k5 = k_d( 1.0, Mtom3( 0.032e9 ), kD )
k6 = 15.0

r1 = createBindingReactionRule( K, KK, K_KK, k1 )
m.network_rules.add_reaction_rule( r1 )
r2 = createUnbindingReactionRule( K_KK, K, KK, k2 )
m.network_rules.add_reaction_rule( r2 )
r3 = createUnbindingReactionRule( K_KK, Kp, KKi, k3 )
m.network_rules.add_reaction_rule( r3 )

r4 = createBindingReactionRule( Kp, KK, Kp_KK, k4 )
m.network_rules.add_reaction_rule( r4 )
r5 = createUnbindingReactionRule( Kp_KK, Kp, KK, k5 )
m.network_rules.add_reaction_rule( r5 )
r6 = createUnbindingReactionRule( Kp_KK, Kpp, KKi, k6 )
m.network_rules.add_reaction_rule( r6 )


r7 = createBindingReactionRule( Kpp, P, Kpp_P, k1 )
m.network_rules.add_reaction_rule( r7 )
r8 = createUnbindingReactionRule( Kpp_P, Kpp, P, k2 )
m.network_rules.add_reaction_rule( r8 )
r9 = createUnbindingReactionRule( Kpp_P, Kp, Pi, k3 )
m.network_rules.add_reaction_rule( r9 )

r10 = createBindingReactionRule( Kp, P, Kp_P, k4 )
m.network_rules.add_reaction_rule( r10 )
r11 = createUnbindingReactionRule( Kp_P, Kp, P, k5 )
m.network_rules.add_reaction_rule( r11 )
r12 = createUnbindingReactionRule( Kp_P, K, Pi, k6 )
m.network_rules.add_reaction_rule( r12 )


r13 = createUnimolecularReactionRule( KKi, KK, ki )
m.network_rules.add_reaction_rule( r13 )
r14 = createUnimolecularReactionRule( Pi, P, ki )
m.network_rules.add_reaction_rule( r14 )

s.setModel(m)


vtklogger = VTKLogger(s, dataDirectory, bufferSize=None, no_shells=no_shells, 
                      extraParticleStep=False, color_dict=color_dict)


if LOGGER == True:
    logname = model
    l = Logger( s, logname = logname )

    rfile = open( 'data/' + logname + '_reactions.dat', 'w' )


    #l.setParticleOutput( ('Ea','X','EaX','Xp','Xpp','EaI') )
    l.setParticleOutInterval( 1e-0 )
    l.log()


#print s.dump_reaction_rules()



counter = 0

currentTime = time.clock()
nextTime = currentTime
endTime = currentTime + CPU_SECONDS

while currentTime < endTime:
    if feq(nextTime, currentTime) == True:
        counter += 1
        nextTime = nextTime + INTERVAL
        vtklogger.log()
    s.step()
    currentTime = time.clock()

    if LOGGER == True:
        if s.lastReaction:
            r = s.lastReaction
            line = '( %18.18g,\t%s,\t%s )\n' % ( s.t, r.reactants, r.products )
            rfile.write( line )
            rfile.flush()

            l.log()

print 'counter = ', counter
vtklogger.stop()

if simulator == 'bd':
    print 'simulation time:', s.t
    print 'number of steps:', s.stepCounter
    print 'bd.dt:', s.dt
elif simulator == 'gfrd':
    s.print_report()
