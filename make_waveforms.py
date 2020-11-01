import lal
import lalsimulation as lalsim

import pylab as pl
import numpy as np


def make_tidal_waveform(approx='TaylorT4', T=4096, Lambda1=None, Lambda2=None,
                        mass1=1.4, mass2=1.3, inclination=0, distance=100,
                        eccentricity=0, meanPerAno=0, phiRef=0, f_min=30,
                        f_ref=0, longAscNodes=0, s1x=0, s1y=0, s1z=0, s2x=0,
                        s2y=0, s2z=0, eos=None):
    # Sanity check
    if (Lambda1 is None) and (Lambda2 is None) and (eos is None):
        # Assuming Lambdas to be zero is not provided
        print("Assuming tidal deformability is zero")
        print("Use arguments Lambda1=, and Lambda2= to provide deformabilities")
        Lambda1 = 0.0
        Lambda2 = 0.0
    mass1 = mass1*lal.MSUN_SI
    mass2 = mass2*lal.MSUN_SI
    distance = distance * 1e6 * lal.PC_SI
    deltaT = 1.0/T
    approximant = lalsim.GetApproximantFromString(approx)
    lal_pars = lal.CreateDict()
    lalsim.SimInspiralWaveformParamsInsertTidalLambda1(lal_pars, Lambda1)
    lalsim.SimInspiralWaveformParamsInsertTidalLambda2(lal_pars, Lambda2)
    hp, hc = lalsim.SimInspiralChooseTDWaveform(mass1, mass2, s1x=s1x, s1y=s1y, 
                                                s1z=s1z, s2x=s2x, s2y=s2y,
                                                s2z=s2z, distance=distance,
                                                inclination=inclination,
                                                phiRef=phiRef,
                                                longAscNodes=longAscNodes,
                                                eccentricity=eccentricity,
                                                meanPerAno=meanPerAno,
                                                deltaT=deltaT, f_min=f_min,
                                                f_ref=f_ref, params=lal_pars,
                                                approximant=approximant)

    return (hp, hc)

def make_waveform_plot(hp, hc):
    tstart = hp.epoch.gpsSeconds + hp.epoch.gpsNanoSeconds*1e-9
    t = np.arange(tstart, 0, hp.deltaT)
    t = t[:len(hp.data.data)]
    pl.plot(t, hp.data.data, 'r-')
    pl.xlabel('Time (s)')
    pl.ylabel('Gravitational wave strain')
    pl.show()




