import os
import sys
import importlib
sys.path.append('/pool/asha0/SCIENCE/csalt/')
import numpy as np
import scipy.constants as sc
from csalt.helpers import *
import matplotlib.pyplot as plt

# plotting style setups
_ = importlib.import_module('plot_setups')
plt.style.use(['default', '/home/sandrews/mpl_styles/nice_line.mplstyle'])



# filenames
ddir = '/data/sandrews/ALMA_regridding/storage/'
truth_MS = 'TRUTH'                                      # oversampled truth
BLC_MS = 'ALMA-BLC_122kHz_SAMPLED'
WSU_MS = 'ALMA-WSU_native_SAMPLED.bin8x'

# colors
tru_col = 'xkcd:cement'

# velocity limits
vlims = [-2.2, 2.2]

# rest frequency
nu0 = 230.538e9

# choose a rough baseline length for a representative spectrum
uv_ = 250e3	# in lambdas (set to 0 if you want min, e.g.)


# load files
tru_dict = read_MS(ddir+truth_MS+'.ms')
BLC_dict = read_MS(ddir+BLC_MS+'.ms')
WSU_dict = read_MS(ddir+WSU_MS+'.ms')

# get (u, v)
u, v = tru_dict['0'].ulam, tru_dict['0'].vlam
uv = np.sqrt(u**2 + v**2)

# choose a spectrum
ix = np.abs(uv - uv_).argmin()

# identify the appropriate timestamp for this index
_stamp = int(tru_dict['0'].tstamp[ix])


# LSRK velocities
tru_vel = 1e-3 * sc.c * (1 - tru_dict['0'].nu_LSRK[_stamp,:] / nu0)
BLC_vel = 1e-3 * sc.c * (1 - BLC_dict['0'].nu_LSRK[_stamp,:] / nu0)
WSU_vel = 1e-3 * sc.c * (1 - WSU_dict['0'].nu_LSRK[_stamp,:] / nu0)

# true visibility spectrum
tru_vis = tru_dict['0'].vis[0,:,ix]
BLC_vis = BLC_dict['0'].vis[0,:,ix]
WSU_vis = WSU_dict['0'].vis[0,:,ix]


# map out plotting based on real, imag peaks
lo_re, hi_re = tru_vis.real.min(), tru_vis.real.max()
max_re = np.abs(tru_vis.real).max()
max_im = np.abs(tru_vis.imag).max()
grd_re = 0.15

re_range = [lo_re - grd_re * max_re, hi_re + grd_re * max_re]
im_range = [-(1 + grd_re) * max_im, (1 + grd_re) * max_im]
hrat = np.diff(re_range)[0] / np.diff(im_range)[0]


fig, ax = plt.subplots(nrows=2, figsize=(3.5, 2.85))

# BLC
ax[1].plot(tru_vel, tru_vis.real, '-', color=tru_col, lw=1.5, alpha=0.8)
ax[1].plot(BLC_vel, BLC_vis.real, 'o', ms=3, color='k')

ax[1].text(0.04, 0.82, 'BLC\n(122 kHz)', transform=ax[1].transAxes, ha='left',
           va='center', color='k', fontsize=7.5)

ax[1].set_xlim(vlims)
ax[0].set_xticklabels([])
ax[1].set_ylim(re_range)
ax[1].set_ylabel('real  $\\mathcal{V}\,(\\nu)$  (Jy)', labelpad=3)


# WSU
ax[0].plot(tru_vel, tru_vis.real, '-', color=tru_col, lw=1.5, alpha=0.8)
ax[0].plot(WSU_vel, WSU_vis.real, 'o', ms=3, color='k')

ax[0].text(0.04, 0.82, 'WSU, x8\n(108 kHz)', transform=ax[0].transAxes, 
           ha='left', va='center', color='k', fontsize=7.5)

ax[0].set_xlim(vlims)
ax[0].set_ylim(re_range)
ax[1].set_xlabel('LSRK  $v - v_{\\rm sys}$  (km s$^{-1}$)')
ax[0].set_ylabel('real  $\\mathcal{V}\,(\\nu)$  (Jy)', labelpad=3)

l, r, t, b = 0.13, 0.87, 0.99, 0.12
fig.subplots_adjust(left=l, right=r, bottom=b, top=t, hspace=0.07)
fig.savefig('figs/demo_SSP.pdf')
fig.clf()
