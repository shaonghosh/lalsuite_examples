import lal
import lalsimulation as lalsim

import pylab as pl
import numpy as np


def make_tidal_waveform(approx='TaylorT4', rate=4096, Lambda1=None, Lambda2=None,
                        mass1=1.4, mass2=1.3, inclination=0, distance=100,
                        eccentricity=0, meanPerAno=0, phiRef=0, f_min=30,
                        f_ref=0, longAscNodes=0, s1x=0, s1y=0, s1z=0, s2x=0,
                        s2y=0, s2z=0, eos=None, save=False):
    # Sanity check
    if (Lambda1 is None) and (Lambda2 is None) and (eos is None):
        # Assuming Lambdas to be zero is not provided
        print("Assuming tidal deformability is zero")
        print("Use arguments Lambda1=, and Lambda2= to provide deformabilities")
        Lambda1 = 0.0
        Lambda2 = 0.0
    if eos:
        if Lambda1 or Lambda2:
            print("Warning: Both eos and Lambda1 and/or Lambda2 has been provided")
            print("Ignoring Lambdas in favor of the eos")
        e = lalsim.SimNeutronStarEOSByName(eos)
        fam = lalsim.CreateSimNeutronStarFamily(e)
        max_mass = lalsim.SimNeutronStarMaximumMass(fam)/lal.MSUN_SI
        assert mass1 < max_mass, "mass1 greater than maximum mass allowed for the neutron star"
        assert mass2 < max_mass, "mass2 greater than the maximum mass allowed for the neutron star"
        r1 = lalsim.SimNeutronStarRadius(mass1*lal.MSUN_SI, fam)
        k1 = lalsim.SimNeutronStarLoveNumberK2(mass1*lal.MSUN_SI, fam)

        r2 = lalsim.SimNeutronStarRadius(mass2*lal.MSUN_SI, fam)
        k2 = lalsim.SimNeutronStarLoveNumberK2(mass2*lal.MSUN_SI, fam)
        c2 = mass2*lal.MRSUN_SI/r2
        c1 = mass1*lal.MRSUN_SI/r1

        Lambda1 = (2/3)*k1/(c1**5)
        Lambda2 = (2/3)*k2/(c2**5)


    mass1 = mass1*lal.MSUN_SI
    mass2 = mass2*lal.MSUN_SI
    distance = distance * 1e6 * lal.PC_SI
    deltaT = 1.0/rate
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
    if save:
        plus_data = hp.data.data
        cross_data = hc.data.data
        tstart_p = hp.epoch.gpsSeconds + hp.epoch.gpsNanoSeconds*1e-9
        tstart_c = hc.epoch.gpsSeconds + hc.epoch.gpsNanoSeconds*1e-9
        tp= np.arange(tstart_p, 0, hp.deltaT)
        tp = tp[:len(hp.data.data)]
        tc = np.arange(tstart_c, 0, hc.deltaT)
        tc = tc[:len(hc.data.data)]
        output_plus = np.vstack((tp, plus_data)).T
        output_cross = np.vstack((tc, cross_data)).T

        np.savetxt("plus_polarization_data.txt", output_plus, fmt="%f\t%e")
        np.savetxt("cross_polarization_data.txt", output_cross, fmt="%f\t%e")


        

    return (hp, hc)

def make_waveform_plot(hp, hc, labels):
    pl.rcParams.update({'font.size': 16})
    pl.figure(figsize=(12,8))
    if type(hp) != list:
        hp = [hp]
    if type(hc) != list:
        hc = [hc]
    fig, (ax1, ax2) = pl.subplots(2,1)
    for x, y, thisLabel in zip(hp, hc, labels):
        tstart_p = x.epoch.gpsSeconds + x.epoch.gpsNanoSeconds*1e-9
        tstart_c = y.epoch.gpsSeconds + y.epoch.gpsNanoSeconds*1e-9
        tp= np.arange(tstart_p, 0, x.deltaT)
        tp = tp[:len(x.data.data)]
        tc = np.arange(tstart_c, 0, y.deltaT)
        tc = tc[:len(y.data.data)]
        ax1.plot(tp, x.data.data, '-', alpha=0.6, label=thisLabel)
        ax1.set_ylabel('Gravitational wave strain')
        ax1.set_title("Plus polarization")
        ax1.legend()

        ax2.plot(tc, y.data.data, '-', alpha=0.6, label=thisLabel)
        ax2.set_ylabel('Gravitational wave strain')
        ax2.set_title("Cross polarization")
        ax2.legend()
    pl.xlabel('Time (s)')
    pl.show()




