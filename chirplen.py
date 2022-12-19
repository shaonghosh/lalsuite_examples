import lal
import numpy as np


def chirplen(m1, m2, fstart):
    """
    METHOD: This a python adaptation of the lalapps_chirplen code.
            `https://lscsoft.docs.ligo.org/lalsuite/lalapps/chirplen_8c_source.html`
            This allows the user to pass the input parameters of the chirplen code
            either as floats or as numpy arrays to compute the length of the signals
    
    """

    mtot = m1 + m2
    eta = (m1*m2)/(mtot**2)
    mchirp = (eta**(0.6))*(mtot)
    fstop = 1.0/(6.0*np.sqrt(6.0)*np.pi*mtot*lal.MTSUN_SI)

    c0 = 5*mtot*lal.MTSUN_SI/(256*eta)
    c2 = 743.0/252.0 + eta*11.0/3.0
    c3 = -32.*np.pi/5.
    c4 = 3058673.0/508032.0 + eta*(5429.0/504.0 + eta*617.0/72.0)
    x  = (np.pi*mtot*lal.MTSUN_SI*fstart)**(1/3)
    x2 = x*x
    x3 = x*x2
    x4 = x2*x2
    x8 = x4*x4
    chirpTime = c0*(1 + c2*x2 + c3*x3 + c4*x4)/x8

    return chirpTime
