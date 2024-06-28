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
sdir = '/data/sandrews/ALMA_regridding/storage/'

truth_MS = 'TRUTH'                                     
dWSU_MS = 'ALMA-WSU_native_SAMPLED.bin8x'
tWSU_MS = 'TRUTH.ALMA-WSU_native.bin8x'
dBLC_MS = 'ALMA-BLC_122kHz_SAMPLED'
tBLC_MS = 'TRUTH.ALMA-BLC_122kHz'

imeth = ['nearest', 'linear', 'cubic', 'fftshift']
nms = ['$\mathtt{nearest}$', '$\mathtt{linear}$',
       '$\mathtt{cubic}$', '$\mathtt{fftshift}$']
cols = ['#F55D3E', '#878E88', '#F7CB15', '#76BED0']
ig = '04'

# colors
tru_col = 'xkcd:cement'

# velocity limits
vlims = [-2.2, 2.2]

# rest frequency
nu0 = 230.538e9

# load truth
tru_dict = read_MS(sdir+truth_MS+'.ms')

# load reference files
tWSU_dict = read_MS(sdir+tWSU_MS+'.REGRID'+ig+'.ms')
tBLC_dict = read_MS(sdir+tBLC_MS+'.REGRID'+ig+'.ms')

# choose a rough baseline length for a representative spectrum
uv_ = 250e3	# in lambdas (set to 0 if you want min, e.g.)

# get (u, v)
u, v = tru_dict['0'].ulam, tru_dict['0'].vlam
uv = np.sqrt(u**2 + v**2)

# choose a spectrum
ix = np.abs(uv - uv_).argmin()

# identify the appropriate timestamp for this index
_stamp = int(tru_dict['0'].tstamp[ix])


# LSRK velocities
tru_vel = 1e-3 * sc.c * (1 - tru_dict['0'].nu_LSRK[_stamp,:] / nu0)
tBLC_vel = 1e-3 * sc.c * (1 - tBLC_dict['0'].nu_LSRK[_stamp,:] / nu0)
tWSU_vel = 1e-3 * sc.c * (1 - tWSU_dict['0'].nu_LSRK[_stamp,:] / nu0)

# true visibility spectrum
tru_vis = tru_dict['0'].vis[0,:,ix]
tBLC_vis = tBLC_dict['0'].vis[0,:,ix]
tWSU_vis = tWSU_dict['0'].vis[0,:,ix]


# map out plotting based on real, imag peaks
lo_re, hi_re = tru_vis.real.min(), tru_vis.real.max()
max_re = np.abs(tru_vis.real).max()
max_im = np.abs(tru_vis.imag).max()
grd_re = 0.15

re_range = [lo_re - grd_re * max_re, hi_re + grd_re * max_re]
im_range = [-(1 + grd_re) * max_im, (1 + grd_re) * max_im]
hrat = np.diff(re_range)[0] / np.diff(im_range)[0]


fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(7, 2.85))

# BLC
ax[0,1].plot(tru_vel, tru_vis.real, '-', color=tru_col, lw=1, alpha=0.5)
ax[0,1].plot(tBLC_vel, tBLC_vis.real, '-', color='xkcd:gunmetal')

# baseline
ax[1,1].axhline(0, color='xkcd:gunmetal')

for im in range(len(imeth)):

    # load interpolated data
    int_dict = read_MS(sdir+dBLC_MS+'.REGRID'+ig+'_'+imeth[im]+'.ms')
    int_vel = 1e-3 * sc.c * (1 - int_dict['0'].nu_LSRK[_stamp,:] / nu0)
    int_vis = int_dict['0'].vis[0,:,ix]

    # plot the interpolates
    ax[0,1].plot(int_vel, int_vis.real, 'o', ms=3, color=cols[im])

    # plot the fractional interpolate residuals
    ax[1,1].plot(int_vel, int_vis.real - tBLC_vis.real, 
                 'o', ms=3, color=cols[im], label=nms[im])


# WSU
ax[0,0].plot(tru_vel, tru_vis.real, '-', color=tru_col, lw=1, alpha=0.5)
ax[0,0].plot(tWSU_vel, tWSU_vis.real, '-', color='xkcd:gunmetal')

# baseline
ax[1,0].axhline(0, color='xkcd:gunmetal')

for im in range(len(imeth)):

    # load interpolated data
    int_dict = read_MS(sdir+dWSU_MS+'.REGRID'+ig+'_'+imeth[im]+'.ms')
    int_vel = 1e-3 * sc.c * (1 - int_dict['0'].nu_LSRK[_stamp,:] / nu0)
    int_vis = int_dict['0'].vis[0,:,ix]

    # plot the interpolates
    ax[0,0].plot(int_vel, int_vis.real, 'o', ms=3, color=cols[im])

    # plot the fractional interpolate residuals
    ax[1,0].plot(int_vel, int_vis.real - tWSU_vis.real, 
                 'o', ms=3, color=cols[im])


# axes and labeling
ax[0,0].text(0.04, 0.82, 'WSU, x8\n(108 kHz)', transform=ax[0,0].transAxes, 
             ha='left', va='center', color='k', fontsize=7.5)
ax[0,1].text(0.04, 0.82, 'BLC\n(122 kHz)', transform=ax[0,1].transAxes,
             ha='left', va='center', color='k', fontsize=7.5)

ax[0,0].set_xlim(vlims)
ax[0,0].set_xticklabels([])
ax[0,0].set_ylim(re_range)
ax[0,0].set_ylabel('real  $\\mathcal{V}\,(\\nu)$  (Jy)', labelpad=6)

ax[1,0].set_xlim(vlims)
ax[1,0].set_ylim([-0.24, 0.24])
ax[1,0].set_xlabel('LSRK  $v - v_{\\rm sys}$  (km s$^{-1}$)')
ax[1,0].set_ylabel('residuals  (Jy)', labelpad=7)

ax[0,1].set_xlim(vlims)
ax[0,1].set_xticklabels([])
ax[0,1].set_ylim(re_range)

ax[1,1].set_xlim(vlims)
ax[1,1].set_xlabel('LSRK  $v - v_{\\rm sys}$  (km s$^{-1}$)')
ax[1,1].set_ylim([-0.24, 0.24])


l, r, t, b = 0.07, 0.93, 0.99, 0.12
fig.subplots_adjust(left=l, right=r, bottom=b, top=t, hspace=0.07, wspace=0.20)
fig.savefig('figs/demo_sregrid_top.pdf')
fig.clf()
