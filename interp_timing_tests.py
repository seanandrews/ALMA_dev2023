import os
import sys
import time
from timeit import default_timer as timer
import importlib
import numpy as np
import matplotlib.pyplot as plt
from regrid_interpolators import *

# plotting style setups
_ = importlib.import_module('plot_setups')
plt.style.use(['default', '/home/sandrews/mpl_styles/nice_line.mplstyle'])


redo_calc = False


# fake line
def gfunc(x, pars):
    return pars[0] * np.exp(-0.5 * (x - pars[1])**2 / pars[2]**2)


Nchan = 2 ** np.arange(3, 20, 1)
imeth = ['nearest', 'linear', 'cubic', 'fftshift1d']
nms = ['$\mathtt{nearest}$', '$\mathtt{linear}$',
       '$\mathtt{cubic}$', '$\mathtt{fftshift}$']
cols = ['#136F63', '#E03AC3', '#F34213', '#3E25FB']
cols = ['#F55D3E', '#878E88', '#F7CB15', '#76BED0']

def_shift = 0.4
Niters = 150


if redo_calc:
    # iterate over number of channels
    ctime = np.empty((len(Nchan), len(imeth)))
    stime = np.empty((len(Nchan), len(imeth)))

    # iterate over interpolation methods
    for im in range(len(imeth)):

        # iterate over Nchan
        for i in range(len(Nchan)):
            print(im, i)

            # define the channel grid
            chix = np.arange(Nchan[i])
            xch = chix - np.mean(chix)
            xsh = xch + def_shift * (xch[1] - xch[0])

            # define the spectrum
            spec = gfunc(xch, [1., 0., 5.])

            # perform interpolation and time the calculation
            t_indiv = np.empty(Niters)
            if imeth[im] == 'fftshift1d':
                t_indiv = np.empty(Niters)
                for j in range(Niters):
                    t0 = time.process_time()
                    ispec = eval("interp_"+imeth[im]+"(xch, spec, def_shift)")
                    t_indiv[j] = time.process_time() - t0
            else:
                t_indiv = np.empty(Niters)
                for j in range(Niters):
                    t0 = time.process_time()	#timer()
                    ispec = eval("interp_"+imeth[im]+"(xch, spec, xsh)")
                    t_indiv[j] = time.process_time() - t0
                      
            # mean, dev of compute times
            ctime[i, im] = np.median(t_indiv)

    np.savez('timings.npz', ctime=ctime)#, stime=stime)



# load calculations
_ = np.load('timings.npz')
ctime = _['ctime'] #, _['stime']

fig, ax = plt.subplots(figsize=(3.5, 2.1))

for im in range(len(imeth)):

    ax.plot(Nchan, 1e3*ctime[:,im], '-', color=cols[im], label=nms[im])
    fint = interp1d(Nchan, 1e3*ctime[:,im], kind='cubic')
    print(fint(8e4))

ax.legend(loc='upper left', prop={'size':7}, framealpha=0)
ax.set_xlim([200, 5e5])
ax.set_ylim([2e-2, 50])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_yticks([0.1, 1, 10])
ax.set_yticklabels(['0.1', '1', '10'])
ax.set_xlabel('number of channels')
ax.set_ylabel('compute time per visibility (ms)')

fig.subplots_adjust(left=0.11, right=0.88, bottom=0.15, top=0.98)
plt.savefig('figs/timings.pdf')
fig.clf()
