import os
import sys
import importlib
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt
from matplotlib.colorbar import Colorbar
from matplotlib.patches import Ellipse
import matplotlib.colors as mcolors
from astropy.io import fits
import cmasher as cmr

# plotting style setups
_ = importlib.import_module('plot_setups')
plt.style.use(['default', '/home/sandrews/mpl_styles/nice_img.mplstyle'])

# set color map
c2 = plt.cm.Reds(np.linspace(0, 1, 32))
c1 = plt.cm.Blues_r(np.linspace(0, 1, 32))
c1 = np.vstack([c1, np.ones((6, 4))])
colors = np.vstack((c1, c2))
mymap = mcolors.LinearSegmentedColormap.from_list('eddymap', colors)
cmap = 'cmr.neutral'


# assign interpolation methods
ig0 = 4
imeth = ['cubic', 'linear', 'nearest']
nms = ['$\mathtt{cubic}$', '$\mathtt{linear}$', '$\mathtt{nearest}$']
col_nms = 'k'

# assign files
sdir = 'MAPS_data/images/'
prefix = 'HD_163296_220GHz_spw6.bin_30s.2x_spw13.REGRID04_'

# frequency selections
nch = 10
ix0, dch = 80, 4
nu0 = 230.538e9

ix0, dch = 42, 2

# set up plotting grids
gs_dat0 = plt.GridSpec(ncols=nch, nrows=1, wspace=0.05, 
                       left=0.01, right=0.99, bottom=0.715, top=0.895)
gs_res0 = plt.GridSpec(ncols=nch, nrows=len(imeth), wspace=0.05, hspace=0.05,
                       left=0.01, right=0.99, bottom=0.01, top=0.575)
fig = plt.figure(figsize=(7, 3.6))
xlims, ylims = [4.2, -4.2], [-4.2, 4.2]

# colorscales
vmin, vmax = 0, 35
rmin, rmax = -3.3, 3.3


# load the reference cube (linear interpolation; as standard)
hdur = fits.open(sdir+prefix+'fftshift.cube.image.fits')
ref_I, ref_h = np.squeeze(hdur[0].data), hdur[0].header
hdur.close()
    
# set the coordinate grids
dx = 3600 * ref_h['CDELT1'] * \
     (np.arange(ref_h['NAXIS1']) - (ref_h['CRPIX1'] - 1))
dy = 3600 * ref_h['CDELT2'] * \
     (np.arange(ref_h['NAXIS2']) - (ref_h['CRPIX2'] - 1))
ddx, ddy = np.meshgrid(dx, dy)
ext = (dx.max(), dx.min(), dy.min(), dy.max())

# set the PSF dimensions
bmaj, bmin, bpa = 0.217 / 3600, 0.151 / 3600, -74.036
bm = (np.pi * bmaj * bmin / (4 * np.log(2))) * (np.pi / 180)**2
print(bm)

# set the frequencies and velocities
nu = ref_h['CRVAL3'] + ref_h['CDELT3'] * \
     (np.arange(ref_h['NAXIS3']) - (ref_h['CRPIX3'] - 1))
vel = sc.c * (1 - nu / nu0)

# plot reference (truth) channel maps
for i in range(nch):
    # desired cube-frame index
    j = ix0 - dch * i 

    # calculate brightness temperatures
    ref_Tb = (1e-26 * ref_I[j,:,:] / bm) * sc.c**2 / (2 * sc.k * nu[j]**2)

    # make the channel map plots
    ax = fig.add_subplot(gs_dat0[i])
    im = ax.imshow(ref_Tb, origin='lower', cmap=cmap, extent=ext, 
                   aspect='equal', vmin=vmin, vmax=vmax)
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_xticks([])
    ax.set_yticks([])
    if i == 0:
        ax.text(0.04, 0.88, '$\mathtt{fftshift}$', ha='left', va='center', 
                color='w', transform=ax.transAxes, fontsize=5)

        beam = Ellipse((xlims[0] + 0.1*np.diff(xlims),
                        -xlims[0] - 0.1*np.diff(xlims)),
                       3600 * bmaj, 3600 * bmin, angle=90-bpa)
        beam.set_facecolor('w')
        ax.add_artist(beam)

    v = vel[j]
    if np.abs(v) < 0.001: v = 0.0
    if np.logical_or(np.sign(v) == 1, np.sign(v) == 0):
        pref = '+'
    else:
        pref = ''
    vstr = pref+'%.2f' % (1e-3 * v)
    ax.text(0.98, 0.05, vstr, transform=ax.transAxes, ha='right',
            va='center', fontsize=5, color='w')


# iterate over interpolation methods
for ia in range(len(imeth)):
    # load the regridded cubes
    hdur = fits.open(sdir+prefix+imeth[ia]+'.cube.image.fits')
    int_I = np.squeeze(hdur[0].data)
    hdur.close()

    for i in range(nch):
        # desired cube-frame index
        j = ix0 - dch * i

        # calculate brightness temperatures
        rTb = (1e-26 * (int_I[j,:,:] - ref_I[j,:,:]) / bm) * \
              sc.c**2 / (2 * sc.k * nu[j]**2)

        # make the channel map plots
        ax = fig.add_subplot(gs_res0[ia,i])
        rim = ax.imshow(rTb, origin='lower', cmap=mymap, extent=ext,
                        aspect='equal', vmin=rmin, vmax=rmax)
        #__ = ax.contour(ddx, ddy, rTb, levels=[-2, 2], colors='k')
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        ax.set_xticks([])
        ax.set_yticks([])
        if i == 0:
            ax.text(0.04, 0.88, nms[ia], ha='left', va='center',
                    color=col_nms, transform=ax.transAxes, fontsize=5)


# colorbars
cbax = fig.add_axes([0.15, 0.905, 0.70, 0.018])
cb = Colorbar(ax=cbax, mappable=im, orientation='horizontal', 
              ticklocation='top', extend='both', extendfrac=0.02)
cbax.tick_params(which='both', labelsize=6, pad=0)
cb.set_label('$T_{\\rm b}$ (K)', fontsize=6, rotation=0, labelpad=3)

cbax = fig.add_axes([0.15, 0.585, 0.70, 0.018])
cb = Colorbar(ax=cbax, mappable=rim, orientation='horizontal',
              ticklocation='top', extend='both', extendfrac=0.02)
cbax.tick_params(which='both', labelsize=6, pad=0)
cb.set_label('residual $T_{\\rm b}$ (K)', fontsize=6, rotation=0, labelpad=3)

figo = 'figs/HD163'
fig.savefig(figo+'.pdf')
print('figure at '+figo+'.pdf')
fig.clf()
