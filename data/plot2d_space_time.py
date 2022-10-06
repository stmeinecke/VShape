import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import gridspec


#plot parameters
DPI = 100
scale = 1.0
scaleCM = 1/2.54
fs = 10*scale
fig = plt.figure(figsize=(5*scale,3.5*scale), dpi=DPI)
plt.rcParams.update({'font.size': fs})



filename = 'out_TS'
data =  np.loadtxt(filename)


CCRTT = 1.0
dt = data[1,0]-data[0,0]
dt = np.around(dt,8)
print(dt)


OSTime = 0.3
OSSteps = int(OSTime / dt)

TSdata = data[:int(2.0/dt),:] 

TS = data[OSSteps:,1]

SPRT = int(CCRTT/dt)
print(SPRT)

Nmax = TS.size
print(Nmax)

RTs = np.floor(Nmax/SPRT)
print(RTs)

space_time_TS = np.reshape(TS[0:int(SPRT*RTs)], (int(RTs),int(SPRT))) 

space_time_TS /= space_time_TS.max()

extentarray = [0, CCRTT, 0,RTs]



#pl1 = plt.imshow(np.flipud(space_time_TS), aspect="auto", cmap="Reds", interpolation="none", extent=extentarray)
pl1 = plt.imshow(np.flipud(space_time_TS), aspect="auto", cmap="Blues", interpolation="none", extent=extentarray, norm=mpl.colors.LogNorm(vmin = 1e-4), rasterized=True)

cb0 = plt.colorbar()
cb0.set_label("Normalized Intensity")

plt.xlabel(r'Time $t$ in $T$')
plt.ylabel(r'Roundtrip Number')


plt.tight_layout()
#plt.savefig(filename+"space_time_TS.pdf")
#plt.savefig(filename+"_space_time_TS.png", dpi=300)
plt.show()
