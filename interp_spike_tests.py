import os
import sys
import importlib
import numpy as np
import matplotlib.pyplot as plt
from interpolators import *
from scipy.ndimage import convolve1d

# plotting style setups
_ = importlib.import_module('plot_setups')
plt.style.use(['default', '/home/sandrews/mpl_styles/nice_line.mplstyle'])


imeth = ['nearest', 'linear', 'cubic', 'fftshift']
nms = ['$\mathtt{nearest}$', '$\mathtt{linear}$',
       '$\mathtt{cubic}$', '$\mathtt{fftshift}$']
cols = ['#136F63', '#E03AC3', '#F34213', '#3E25FB']
cols = ['#F55D3E', '#878E88', '#F7CB15', '#76BED0']

Nshifts = 100
shifts = np.arange(-Nshifts, Nshifts) / Nshifts


# create the input model: a spike injected post-SRF
ch = np.arange(-100, 100, 1)
s0 = np.zeros_like(ch)
s0[ch == 0] = 1

fig, ax = plt.subplots(nrows=2, figsize=(3.5, 2.85))


# iterate over type of spike
for ii in range(2):

    if ii == 0:
        ispec0 = 1 * s0
    else:
        blc_SRF = np.array([0.25, 0.5, 0.25])
        ispec0 = np.convolve(s0, blc_SRF / np.sum(blc_SRF), mode='same')

    # iterate over interpolation methods
    for im in range(len(imeth)):

        # iterate over sub-channel shifts
        xch_ = []
        ispec_ = []
        for i in range(len(shifts)):
            # define the interpolate channel grid
            xch = ch + shifts[i] 

            # perform interpolation and time the calculation
            if imeth[im] == 'fftshift':
                ispec = eval("interp_"+imeth[im]+"(ch, ispec0, shifts[i])")
                xch -= 2 * shifts[i]
            else:
                ispec = eval("interp_"+imeth[im]+"(ch, ispec0, xch)")

            # append
            xch_ += [xch]
            ispec_ += [ispec]

        # concatenate for accumulated shifts
        xch_all = np.concatenate(xch_)
        ispec_all = np.concatenate(ispec_)
        xch_s = xch_all[np.argsort(xch_all)]
        ispec_s = ispec_all[np.argsort(xch_all)]
   
        ax[ii].plot(xch_s, ispec_s, '-', color=cols[im], label=nms[im])#, lw=0.8)
    ax[ii].plot(ch, ispec0, 'ok', ms=2)

    ax[ii].set_xlim([-10, 10])
    ax[ii].set_ylim([-0.3, 1.15])
    ax[ii].set_xticks([-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10]) 
    if ii < 1:
        ax[ii].set_xticklabels([])
        ax[ii].set_yticklabels([])

ax[0].text(0.95, 0.8, 'WSU', transform=ax[0].transAxes, ha='right', 
           va='center', color='k')
ax[1].text(0.95, 0.8, 'BLC', transform=ax[1].transAxes, ha='right', 
           va='center', color='k')
ax[1].set_xlabel('channel offsets in units of $\\Delta \\nu$')
ax[1].set_ylabel('interpolated spectra')#, labelpad=0)
ax[0].legend(loc='upper left', prop={'size':6}, framealpha=0)

fig.subplots_adjust(left=0.13, right=0.87, bottom=0.12, top=0.99, hspace=0.07)
plt.savefig('figs/spikes.pdf')
fig.clf()
