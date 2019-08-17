import sys
import numpy as np

import lal
import lalsimulation as lalsim


def getPSD(deltaF, npts, psd="SimNoisePSDaLIGOZeroDetHighPowerPtr"):
    """
    This function accepts a frequency step-size and number of frequency points. 
    It then computes from there the noise PSD for ligo. The default psd is 
    SimNoisePSDaLIGOZeroDetHighPowerPtr which can be changed using the argument
    psd.

    deltaF :: The frequency interval between data points.
    npts   :: Total number of date points
    psd    :: The noise model in lalsimulation. To get info use
              lalsim.__dict__ in interpretor and search string SimNoisePSD
    """
    lalseries = lal.CreateREAL8FrequencySeries('', lal.LIGOTimeGPS(0),
                                               0, deltaF, lal.DimensionlessUnit,
                                               npts)
    func = lalsim.__dict__[psd]
    lalsim.SimNoisePSD(lalseries, 0, func)
    f_end = lalseries.f0 + len(lalseries.data.data)*lalseries.deltaF
    freq = np.linspace(lalseries.f0, f_end, len((lalseries.data.data)))
    output = np.vstack((freq, lalseries.data.data)).T

    return output

