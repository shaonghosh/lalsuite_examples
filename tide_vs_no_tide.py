import numpy as np
import matplotlib.pyplot as plt
import make_waveforms

srate = 16384*2
tidal = make_waveforms.make_tidal_waveform(rate=srate, Lambda1=6000, Lambda2=8000)
nontidal = make_waveforms.make_tidal_waveform(rate=srate)

start_time_tidal = tidal[0].epoch.gpsSeconds + 1e-9*(tidal[0].epoch.gpsNanoSeconds)
start_time_nontidal = nontidal[0].epoch.gpsSeconds + 1e-9*(nontidal[0].epoch.gpsNanoSeconds)

times_tidal = np.arange(start_time_tidal, 0, tidal[0].deltaT)
times_nontidal = np.arange(start_time_nontidal, 0, nontidal[0].deltaT)

print("Non-Tidal")
print(len(times_nontidal))
print(len(nontidal[0].data.data))

print("Tidal")
print(len(times_tidal))
print(len(tidal[0].data.data))

plt.plot(times_tidal[:-1], tidal[0].data.data, label="Tidal effects ON")
plt.plot(times_nontidal[:-1], nontidal[0].data.data, label="Tidal effects OFF")
plt.legend()
plt.show()


