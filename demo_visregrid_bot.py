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
truth_MS = 'TRUTH'                                      # oversampled truth
dWSU_MS = 'ALMA-WSU_native_SAMPLED.bin8x'
dBLC_MS = 'ALMA-BLC_122kHz_SAMPLED'
tWSU_MS = 'TRUTH.ALMA-WSU_native.bin8x'
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

# load reference files
tWSU_dict = read_MS(sdir+tWSU_MS+'.REGRID'+ig+'.ms')
tBLC_dict = read_MS(sdir+tBLC_MS+'.REGRID'+ig+'.ms')

# calculate reference velocities
tBLC_vel = 1e-3 * sc.c * (1 - tBLC_dict['0'].nu_LSRK / nu0)
tWSU_vel = 1e-3 * sc.c * (1 - tWSU_dict['0'].nu_LSRK / nu0)

# identify frequency indices corresponding to specified vlims
ixl_BLC = np.min(np.argmin(np.abs(tBLC_vel - vlims[1]), axis=1))
ixh_BLC = np.max(np.argmin(np.abs(tBLC_vel - vlims[0]), axis=1))
ixl_WSU = np.min(np.argmin(np.abs(tWSU_vel - vlims[1]), axis=1))
ixh_WSU = np.max(np.argmin(np.abs(tWSU_vel - vlims[0]), axis=1))

# select subset of reference visibility spectra
tBLC_vis = tBLC_dict['0'].vis[0,ixl_BLC:ixh_BLC+1,:]
tWSU_vis = tWSU_dict['0'].vis[0,ixl_WSU:ixh_WSU+1,:]

duv = np.arange(0, 2000, 100)


fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(7, 1.65))

for im in range(len(imeth)):

    # load interpolated data
    int_dict = read_MS(sdir+dBLC_MS+'.REGRID'+ig+'_'+imeth[im]+'.ms')
    int_vis = int_dict['0'].vis[0,ixl_BLC:ixh_BLC+1,:]
    int_uv = np.tile(np.sqrt(int_dict['0'].ulam**2 + int_dict['0'].vlam**2), 
                     (int_vis.shape[0], 1)).flatten()
    int_resid = np.absolute((int_vis - tBLC_vis) / tBLC_vis).flatten()
    xx = np.digitize(1e-3 * int_uv, duv)
    yy = np.empty(len(duv))
    for ix in range(len(duv)):
        yy[ix] = np.median(int_resid[xx == ix])
            
    # plot the scatter as a function of baseline length
    ax[1].plot(duv, yy, '-', color=cols[im], label=nms[im])

    # load interpolated data
    int_dict = read_MS(sdir+dWSU_MS+'.REGRID'+ig+'_'+imeth[im]+'.ms')
    int_vis = int_dict['0'].vis[0,ixl_WSU:ixh_WSU+1,:]
    int_uv = np.tile(np.sqrt(int_dict['0'].ulam**2 + int_dict['0'].vlam**2), 
                     (int_vis.shape[0], 1)).flatten()
    int_resid = np.absolute((int_vis - tWSU_vis) / tWSU_vis).flatten()
    xx = np.digitize(1e-3 * int_uv, duv)
    yy = np.empty(len(duv))
    for ix in range(len(duv)):
        yy[ix] = np.median(int_resid[xx == ix])

    # plot the scatter as a function of baseline length
    ax[0].plot(duv, yy, '-', color=cols[im], label=nms[im])


ax[0].legend(loc='upper left', prop={'size':6}, framealpha=0, 
             handletextpad=0.5, ncols=2)

ax[0].set_xlim([0, 2000])
ax[0].set_xticks([0, 500, 1000, 1500, 2000])
ax[0].set_xlabel('baseline length  (k$\\lambda$)')
ax[0].set_ylim([0, 0.55])
ax[0].set_ylabel('$\\langle | (\\mathcal{V}_{\\rm obs} - \\mathcal{V}_{\\rm ref}) \,\, / \,\, \\mathcal{V}_{\\rm ref} | \\rangle$', labelpad=10)

ax[1].set_xlim([0, 2000])
ax[1].set_xticks([0, 500, 1000, 1500, 2000])
ax[1].set_xlabel('baseline length  (k$\\lambda$)')
ax[1].set_ylim([0, 0.55])


l, r, t, b = 0.07, 0.93, 0.94, 0.20
fig.subplots_adjust(left=l, right=r, bottom=b, top=t, wspace=0.20)
fig.savefig('figs/demo_sregrid_bot.pdf')
fig.clf()
