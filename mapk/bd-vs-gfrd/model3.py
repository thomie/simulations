#!/usr/bin/env python

# python -O model3.py gfrd
# python -O model3.py bd

from logger import *
from visualization import vtklogger
from model import *
import sys
import time

import math



LOGGER = False

# First 1000 frames (50s) of slow motion (zooming out).
INTERVAL_SIMULATION_TIME = 0.000001
TOTAL_SIMULATION_TIME = 0.001
INITIAL_SIMULATION_TIME = 0
INITIAL_VTK_STEP = 0

V = 1e-15
D_ratio = 1
D_mode = 'normal'
ti = 0 # Lot's of rebindings.
ti = 1e-6 # Bifocation point.
KK_VS_P = 1.0

simulator = sys.argv[1]



# V in liter, L in meter
L = math.pow(V * 1e-3, 1.0 / 3.0)
N = 180
matrixSize = min(max(3, int((3 * N) ** (1.0/3.0))), 60)
#print 'matrixSize=', matrixSize


if simulator == 'gfrd':
    from egfrd import *
    dataDirectory = 'gfrd-'
    show_shells = True
elif simulator == 'bd':
    from bd import *
    w = World(L, matrixSize)
    s = BDSimulator(w)
    #s.core.dtFactor = DEFAULT_DT_FACTOR
    dataDirectory = 'bd-'
    show_shells = False

dataDirectory += str(TOTAL_SIMULATION_TIME) + '_s-' + str(INTERVAL_SIMULATION_TIME) + '_s_intervals-' + str(KK_VS_P) + '_kk_vs_p-zoomout'
#str(CPU_SECONDS) + '_cpu_s-' 


model='mapk3'

if ti == 0:
    ki = float('inf')
else:
    ki = math.log(2) / ti

D_ref = 1e-12

D_move = D_ref * D_ratio

if D_mode == 'normal':
    D_react = D_move
elif D_mode == 'fixed':
    D_react = D_ref

radius = 2.5e-9

m = ParticleModel(L)

K = Species('K', D_move, radius)
KK = Species('KK', D_move, radius)
P = Species('P', D_move, radius)
Kp = Species('Kp', D_move, radius)
Kpp = Species('Kpp', D_move, radius)
K_KK = Species('K_KK', D_move, radius)
Kp_KK = Species('Kp_KK', D_move, radius)
Kpp_P = Species('Kpp_P', D_move, radius)
Kp_P = Species('Kp_P', D_move, radius)

# inactive forms
KKi = Species('KKi', D_move, radius)
Pi = Species('Pi', D_move, radius)

m.add_species_type(K)
m.add_species_type(KK)
m.add_species_type(P)
m.add_species_type(Kp)
m.add_species_type(Kpp)
m.add_species_type(K_KK)
m.add_species_type(Kp_KK)
m.add_species_type(Kpp_P)
m.add_species_type(Kp_P)
m.add_species_type(KKi)
m.add_species_type(Pi)

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


if simulator == 'gfrd':
    w = create_world(m, matrixSize)
    s = EGFRDSimulator(w)
elif simulator == 'bd':
    w = create_world(m, matrixSize)
    s = BDSimulator(w)

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
kD = k_D(D_react * 2, sigma)

N_K = C2N(200e-9, V) 
N_Kpp = 0

#N_K = 0
#N_Kpp = C2N(200e-9, V) 

N_KK = KK_VS_P * 2 * C2N(50e-9, V)
N_P = (1 - KK_VS_P) * 2 * C2N(50e-9, V)


throw_in_particles(s.world, K, N_K)
print 'Number of K =', N_K
throw_in_particles(s.world, Kpp, N_Kpp)
print 'Number of Kpp =', N_Kpp
throw_in_particles(s.world, KK, N_KK)
print 'Number of KK =', N_KK
throw_in_particles(s.world, P, N_P)
print 'Number of P =', N_P

if simulator == 'gfrd':
    # Initialize scheduler, cut gfrd some slack.
    endTime = 0
    while 1:
        s.step()
        nextTime = s.scheduler.top[1].time
        if nextTime > endTime:
            s.stop(endTime)
            break
s.reset()

k1 = k_a(per_M_to_m3(0.02e9), kD)
k2 = k_d(1.0, per_M_to_m3(0.02e9), kD)
k3 = 1.5
k4 = k_a(per_M_to_m3(0.032e9), kD)
k5 = k_d(1.0, per_M_to_m3(0.032e9), kD)
k6 = 15.0

r1 = create_binding_reaction_rule(K, KK, K_KK, k1)
m.network_rules.add_reaction_rule(r1)
r2 = create_unbinding_reaction_rule(K_KK, K, KK, k2)
m.network_rules.add_reaction_rule(r2)
r3 = create_unbinding_reaction_rule(K_KK, Kp, KKi, k3)
m.network_rules.add_reaction_rule(r3)

r4 = create_binding_reaction_rule(Kp, KK, Kp_KK, k4)
m.network_rules.add_reaction_rule(r4)
r5 = create_unbinding_reaction_rule(Kp_KK, Kp, KK, k5)
m.network_rules.add_reaction_rule(r5)
r6 = create_unbinding_reaction_rule(Kp_KK, Kpp, KKi, k6)
m.network_rules.add_reaction_rule(r6)


r7 = create_binding_reaction_rule(Kpp, P, Kpp_P, k1)
m.network_rules.add_reaction_rule(r7)
r8 = create_unbinding_reaction_rule(Kpp_P, Kpp, P, k2)
m.network_rules.add_reaction_rule(r8)
r9 = create_unbinding_reaction_rule(Kpp_P, Kp, Pi, k3)
m.network_rules.add_reaction_rule(r9)

r10 = create_binding_reaction_rule(Kp, P, Kp_P, k4)
m.network_rules.add_reaction_rule(r10)
r11 = create_unbinding_reaction_rule(Kp_P, Kp, P, k5)
m.network_rules.add_reaction_rule(r11)
r12 = create_unbinding_reaction_rule(Kp_P, K, Pi, k6)
m.network_rules.add_reaction_rule(r12)


r13 = create_unimolecular_reaction_rule(KKi, KK, ki)
m.network_rules.add_reaction_rule(r13)
r14 = create_unimolecular_reaction_rule(Pi, P, ki)
m.network_rules.add_reaction_rule(r14)


vtklogger = vtklogger.VTKLogger(s, dataDirectory, buffer_size=None,     
                                show_shells=show_shells, 
                                extra_particle_step=False,
                                color_dict=color_dict)


if LOGGER == True:
    logname = model
    l = Logger(s, logname = logname)

    rfile = open('data/' + logname + '_reactions.dat', 'w')


    #l.setParticleOutput(('Ea','X','EaX','Xp','Xpp','EaI'))
    l.setParticleOutInterval(1e-0)
    l.log()


#print s.dump_reaction_rules()



counter = 0
next_simulation_time = INITIAL_SIMULATION_TIME

s.t = INITIAL_SIMULATION_TIME
vtklogger.i = INITIAL_VTK_STEP

while s.t <= INITIAL_SIMULATION_TIME + TOTAL_SIMULATION_TIME +  INTERVAL_SIMULATION_TIME:
    if s.t > next_simulation_time:
        counter += 1
        interval= INTERVAL_SIMULATION_TIME

        '''
        if s.t > TOTAL_SIMULATION_TIME1:
            interval = INTERVAL_SIMULATION_TIME1 * pow(INTERVAL_SIMULATION_FACTOR, (counter - 600))
        else:
            interval= INTERVAL_SIMULATION_TIME1 
        '''

        print '%d %f %f' % (counter, interval, next_simulation_time)

        next_simulation_time += interval
        vtklogger.log()
    s.step()

    if LOGGER == True:
        if s.lastReaction:
            r = s.lastReaction
            line = '(%18.18g,\t%s,\t%s)\n' % (s.t, r.reactants, r.products)
            rfile.write(line)
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
