#!/usr/bin/python
# Requirements on Ubuntu: python-sympy or python-mpmath.

# test script to sample a distribution of first-passage times
#
# the setup is: one absorbing, one reflecting plane which are parallel. Particles are started
# close to the absorbing plane, all at a given distance to it.
# Multiple particles are started at the same time to incorporate also pair domains.
# The survival times are collected and histogrammed; they are then compared to the known
# analytical result
#

# general imports 
import os
import time
import shutil

# gfrd
from egfrd import *
from logger import *

# for calculating and sampling
import numpy as np
import matplotlib.pyplot as plt
try:
    import mpmath as m
except ImportError:
    import sympy.mpmath as m


# settings and parameters for the sampling run
#
# files
# --------------ADAPT THIS!-----------
DATADIR='/home/thomas/code/runs/fp_sampling'
DATAFN='sample'
#
# outer loop count. each time there are 9 particles running in parallel.
# therefore, the total number of counts in the histogram is about 9*NUMEROFRUNS
NUMBEROFRUNS = int(1e3)
#
# upper cutoff time for sampling: if a particle stays alive longer than this time, 
# it is discarded
CUTOFFTIME = .2
# the geometry and kinetic constants, collected here for easy reference
# these values may be adapted when using the calculated theoretical curve.
#
# box size
L = 5e-6    
# the run length of the effective 1-d diffusion in z direction is 
# slabThickness - 2*stopThickness/2 -2*rA (see below)
slabThickness = L / 2
# attention: it seems that Lz is really only half the thickness! centered!
stopThickness = 2e-8
# sort of standard values for the particles: 1mu^2/s, 10nm size.
DA = 1e-12
rA = 1e-8#
# association, irreversible. This is supposed to be essentially INF.
kAss = 1e5
# whether to include a fake species on the plane that decays later. Obsolete since
# direct absorption is implemented.
DUMMYDECAY = False  
# decay on membrane
if DUMMYDECAY: kDec = 1e3
# starting distance
distToAbsorbing = 1e-7
#
# settings and parameters for the histogram and plot
#
# bins
HISTBINS = 25
# plot?
MAKEFIGURE = True


def setupSim(l):        
    # setting up the simulator; has to be done in each loop to get proper reset.
    s = EGFRDSimulator(worldSize=L)
    # logging
    l.sim = s
    # geometry
    #
    # For PlanarSurface, the first argument is the center of mass
    # this is a list of the two surfaces.   
    stopSurf = dict([(i[1], s.addPlanarSurface(origin=[L / 2, L / 2, L / 2 + i[0] * slabThickness],
                                vectorX=[1, 0, 0],
                                vectorY=[0, 1, 0],
                                Lx=L / 2,
                                Ly=L / 2,
                                Lz=stopThickness / 2,
                                name=i[1])) for i in [(0.5, 'up'), (-0.5, 'down')] ])
    # Species
    # create particles in the bulk
    A = Species('A', DA, rA)
    s.addSpecies(A)           
    if DUMMYDECAY:
        for su in stopSurf.itervalues():
            s.addSpecies(A, su, 0.1 * DA, rA) # the adsorbed species: smaller diffusion constant on membrane
    # Reactions
    # absorption only for the last surface, 'down'
    for su in stopSurf['down'], : # ugly but it works
        if DUMMYDECAY:
            s.addReaction([A], [(A, su)], kAss)
            s.addReaction([(A, su)], [(0, su)], kDec)
        else:
            s.addReaction([A], [(0, su)], kAss)
    #
    # a grid of 9 particles. by default, they are far enough apart to not generate
    # multis. it may be nice to test those too but that just seems to take too long.
    for i in range(3): 
        for j in range(3):
            s.placeParticle(A, [L / 3 * j, L / 3 * i, stopSurf['down'].origin[2] + distToAbsorbing ])
    return s  


def run():
    '''sample the arrival times using the global settings.
    
    '''
    # initialization of the data directory.
    shutil.rmtree(DATADIR, ignore_errors=True)

    '''
    inner loop
    '''
    def onerun(s, l):
        old = s.t
        new = old
        startTime = old
        # this is the upper cutoff time.
        while (new - startTime) < CUTOFFTIME  :
            try:
                s.step()
            # when there are no more particles:
            except RuntimeError, message:
                #print message        
                break
            # in any case:
            finally:
                 l.logTimeCourse()
            new = s.t
            old = new
    '''
    Simulation.
    '''
    # outer loop
    ttt = time.time()
    dummysim = EGFRDSimulator(L)
    log = Logger(dummysim , DATAFN, DATADIR, 'first passage sampling')
    for i in range(NUMBEROFRUNS):
        if i % 10 == 0:
            print i
        try:
            sim = setupSim(log)
            onerun(sim, log)
            del sim
        except RuntimeError as e:
            print e
    log.timecourseFile.close()
    print time.time() - ttt



'''
making the histograms
'''


class aNamespace():
    pass

hi = aNamespace()

# utility
def logbins(c, n):
    '''n bins that are uniform in the log domain and span the range of c.
    '''
    mi, ma = np.log((min(c), max(c)))
    dm = (ma - mi) / float(n)
    return np.exp(np.arange(mi, ma, dm)) 


# take the histogram   
def makeHist():
    '''take the histogram from the sampled file
    '''
    hi.timetable = np.loadtxt(DATADIR + os.sep + DATAFN + '_tc.dat', skiprows=3)[:, 0]
    hi.timetable.sort()
    hi.cts, hi.bins = np.histogram(hi.timetable, bins=logbins(hi.timetable, HISTBINS))
    # dividing by the bin widths to get the right pdf
    hi.cts = hi.cts / np.diff(hi.bins)
    # bin midpoints and normalization to 1
    hi.binmids = np.convolve(hi.bins, (.5, .5), mode='valid')
    hi.rcts = hi.cts / np.trapz(hi.cts, hi.binmids)
    return (hi.binmids, hi.rcts)

# plot it


pn = aNamespace()

def makeFig():
    pn.f3 = plt.figure()
    pn.ax3 = pn.f3.add_subplot(111)
    pn.ax3.set_xscale('log')
    pn.ax3.set_yscale('log')
    pn.ax3.set_xlabel('hitting time (s)')
    pn.ax3.set_ylabel('PDF (1/s)')

    
def plotHist(hist):
    histline, = pn.ax3.plot(*hist)
    histline.set_marker('o')


    

'''
calculating the theoretical histogram
'''


# this function does in principle the right thing.
# however, for small arguments q to jtheta, we 
# need a lot of extra precision for evaluation.
# this is probably the same numerical effect as seen for 
# the inverse laplace transformed fp densities.
def calcHist():
    # ridiculous amount of extra precision is needed for the shortest times
    @m.extradps(1000) 
    def arFpDensity(Da, r0a, La, ta):
        D, r0, L, t = (m.mpmathify(str(x)) for x in (Da, r0a, La, ta))
        pref = m.exp(-r0 ** 2 / (4 * D * t)) * D / (4 * m.sqrt(m.pi) * (D * t) ** 1.5)
        z = (1j * L * r0) / (2 * D * t)
        q = m.exp(-L ** 2 / (D * t))
        return float(m.re(2 * pref * (r0 * m.jtheta(4, z, q) - 0 * 1j * L * m.djtheta(4, z, q))))
    runlength = slabThickness - 2 * stopThickness / 2 - 2 * rA
    initsep = distToAbsorbing - stopThickness / 2 - rA
    theolist = np.array([arFpDensity(DA, initsep, runlength, x)  for x in hi.binmids])
    renormlist = theolist / np.trapz(theolist, hi.binmids)
    return (hi.binmids, renormlist)


'''
The actual run
''' 

if __name__ == '__main__':
    run()    
    sampledhist = makeHist()
    calcdhist = calcHist()
    res = sampledhist[ - 1] - calcdhist[ - 1]
    chi2 = sum(res ** 2)
    reldev = np.sqrt(chi2) / len(hi.binmids)
    print 'bin timepoints, residues:'
    print (hi.binmids, res)
    print 'sum of squared deviations:'
    print chi2
    print 'rmsd averaged over sampled bins with equal weight:'
    print reldev
    if MAKEFIGURE: 
        makeFig()
        plotHist(sampledhist)
        plotHist(calcdhist)
        plt.show()
