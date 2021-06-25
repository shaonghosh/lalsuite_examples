import make_waveforms
import pylab as pl
import numpy as np


hp,hc= make_waveforms.make_tidal_waveform(mass1= 10, mass2 = 8, s1z=0, approx = "SpinTaylorT4")

hp2,hc2=make_waveforms.make_tidal_waveform(mass1= 10, mass2 = 8, s1z = 0.5, approx = "SpinTaylorT4")

hp3,hc3=make_waveforms.make_tidal_waveform(mass1 = 10, mass2 = 8, s1z = -0.5, approx = "SpinTaylorT4")

def get_times(hp):
    t=[]
    for t in range (0,10):
        t=np.arange(0, -hp.epoch.gpsSeconds - hp.epoch.gpsNanoSeconds*1e-9 + hp.deltaT, hp.deltaT)
    t=t[:hp.data.length]
    return t
pl.plot(get_times(hp), hp.data.data, label= "Without Spin")
    
    #t=np.arange(0, -hp.epoch.gpsSeconds -hp.epoch.gpsNanoSeconds*1e-9 +hp.deltaT, hp.deltaT)
    #t=t[:hp.data.length]
    #return t
pl.plot(get_times(hp2), hp2.data.data, label= "With Spin")
pl.plot(get_times(hp3), hp3.data.data, label= "With Negative Spin")
pl.xlabel("Time (s)")
pl.ylabel("Strain")
pl.legend()
pl.savefig("MultiPlot1.png") 
