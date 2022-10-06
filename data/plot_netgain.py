from matplotlib import rc
rc('text', usetex=True)

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import gridspec

import matplotlib.ticker as ticker


#plot parameters
DPI = 100
scale = 1.0
fs = 8*scale
fig = plt.figure(figsize=(5.74916*scale,3.0*scale))
plt.rcParams.update({'font.size': fs})
gs = gridspec.GridSpec(2,2, width_ratios=[1,1], height_ratios=[2.5,1])


gridcolor = '0.8'

labeltextx = 0.01
labeltexty = 0.98

tmin = 0
tmax = 3

pcolor='C0'
Gcolor='g'
Qcolor='C3'
NGcolor='purple'



##################################################################################################################

##################################################################################################################
ax00 = plt.subplot(gs[0,:])
plt.grid(c=gridcolor)


dataFML = np.loadtxt('out_simpleTS')
tau01 = 0.25


Qshift = 5.0
Gshift = Qshift * (tau01/0.5)


tshift = 0.0
dt = dataFML[1,0] - dataFML[0,0]
kmax = int(1.1/dt)
#print kmax
Imax = 0.0
kpos = 0
for k in range(kmax):
  if dataFML[k,1] > Imax:
    Imax = dataFML[k,1]
    kpos = k

#print k
#print Imax
tshift = kpos*dt+0.5
#print tshift



plt.plot(np.array([0.0,0.5,1.0,1.5,2.0,2.5,3.0]),[Qshift,0.0,Qshift,0.0,Qshift,0.0,Qshift], ls='--', color='C0', alpha=0.5)

plt.axhline(y=0.0, ls=':', c='k')
plt.plot(dataFML[:,0]-dataFML[0,0]-tshift,dataFML[:,1]/dataFML[:,1].max(),color=pcolor)

plt.axhline(y=Gshift, ls=':', c='k')
plt.plot(dataFML[:,0]-dataFML[0,0]-tshift,dataFML[:,2]/np.abs(dataFML[:,2]).max()+Gshift,color=Gcolor)

plt.axhline(y=Qshift, ls=':', c='k')
plt.plot(dataFML[:,0]-dataFML[0,0]-tshift,dataFML[:,3]/np.abs(dataFML[:,3]).max()+Qshift,color=Qcolor)




plt.xlim(tmin,tmax)
plt.ylim(-0.2,5.7)


#plt.ylabel(r'Intensity $|E|^2$')

def multicolor_ylabel(ax,list_of_strings,list_of_colors,axis='x',anchorpad=0,**kw):
    """this function creates axes labels with multiple colors
    ax specifies the axes object where the labels should be drawn
    list_of_strings is a list of all of the text items
    list_if_colors is a corresponding list of colors for the strings
    axis='x', 'y', or 'both' and specifies which label(s) should be drawn"""
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

    # x-axis label
    if axis=='x' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',**kw)) 
                    for text,color in zip(list_of_strings,list_of_colors) ]
        xbox = HPacker(children=boxes,align="center",pad=0, sep=5)
        anchored_xbox = AnchoredOffsetbox(loc=3, child=xbox, pad=anchorpad,frameon=False,bbox_to_anchor=(0.2, -0.09),
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_xbox)

    # y-axis label
    if axis=='y' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='right',va='bottom',rotation=90,**kw)) 
                     for text,color in zip(list_of_strings[::-1],list_of_colors[::-1]) ]
        ybox = VPacker(children=boxes,align="center", pad=0, sep=5)
        anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=(-0.08, 0.1), 
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_ybox)

#multicolor_ylabel(ax20,('Line1','and','Line2','with','extra','colors!'),('r','k','b','k','m','g'),axis='y',size=15,weight='bold')
multicolor_ylabel(ax00,('Intensity','Gain','Absorption'),(pcolor, Gcolor, Qcolor),axis='y',size=fs)


ax00.xaxis.set_ticks_position('both')
ax00.yaxis.set_ticks_position('both')
ax00.tick_params(which='both', direction='in')

ax00.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax00.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
#ax00.xaxis.set_major_formatter(ticker.NullFormatter())
#ax00.xaxis.set_minor_formatter(ticker.NullFormatter())

ax00.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax00.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
#ax00.yaxis.set_major_formatter(ticker.NullFormatter())
#ax00.yaxis.set_minor_formatter(ticker.NullFormatter())


plt.title(r'FML', fontsize=fs)
plt.text(labeltextx, labeltexty,r'(a)', horizontalalignment='left', verticalalignment='top', transform=ax00.transAxes)

###########################################

ax10 = plt.subplot(gs[1,:])
plt.grid(c=gridcolor)


plt.plot(dataFML[:,0]-dataFML[0,0]-tshift,dataFML[:,4],color=NGcolor)

plt.xlim(tmin,tmax)

plt.ylabel(r'Netgain $\mathcal{G}$')
plt.xlabel(r'Time $t$ in $T$')

ax10.xaxis.set_ticks_position('both')
ax10.yaxis.set_ticks_position('both')
ax10.tick_params(which='both', direction='in')

ax10.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax10.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
#ax10.xaxis.set_major_formatter(ticker.NullFormatter())
#ax10.xaxis.set_minor_formatter(ticker.NullFormatter())

#ax10.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
#ax10.yaxis.set_minor_locator(ticker.MultipleLocator(0.04))
#ax10.yaxis.set_major_formatter(ticker.NullFormatter())
#ax10.yaxis.set_minor_formatter(ticker.NullFormatter())

plt.text(labeltextx, labeltexty,r'(c)', horizontalalignment='left', verticalalignment='top', transform=ax10.transAxes)




###################################################################################################################
###################################################################################################################

plt.tight_layout()
#plt.subplots_adjust(left=0.08, bottom=0.13, right=0.98, top=0.94, hspace=0.1, wspace=0.15)


#strFigName = "netgain"
#plt.savefig(strFigName+".pdf")
#plt.savefig(strFigName+".png", dpi=600)

plt.show()
