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
ig = '04'
imeth = ['nearest', 'linear', 'cubic', 'fftshift']
nms = ['$\mathtt{nearest}$', '$\mathtt{linear}$',
       '$\mathtt{cubic}$', '$\mathtt{fftshift}$']
col_nms = 'k'

# assign files
sdir = '/data/sandrews/ALMA_regridding/storage/images/'
tfile_BLC = sdir+'TRUTH.ALMA-WSU_native.bin8x.REGRID'+ig
dfile_BLC = sdir+'ALMA-WSU_native_SAMPLED_postavg.bin8x.REGRID'+ig
tfile_WSU = sdir+'TRUTH.ALMA-WSU_native.bin8x.REGRID'+ig
dfile_WSU = sdir+'ALMA-WSU_native_SAMPLED.bin8x.REGRID'+ig

# frequency selections
nch = 5
ix0_blc, dch_blc = 52, 3
ix0_wsu, dch_wsu = 52, 3
nu0 = 230.538e9

# set up plotting grids
gs_dat1 = plt.GridSpec(ncols=nch, nrows=1, wspace=0.05, 
                       left=0.01, right=0.485, bottom=0.755, top=0.905)
gs_dat0 = plt.GridSpec(ncols=nch, nrows=1, wspace=0.05,
                       left=0.515, right=0.99, bottom=0.755, top=0.905)
gs_res1 = plt.GridSpec(ncols=nch, nrows=len(imeth), wspace=0.05, hspace=0.05,
                       left=0.01, right=0.485, bottom=0.01, top=0.63)
gs_res0 = plt.GridSpec(ncols=nch, nrows=len(imeth), wspace=0.05, hspace=0.05,
                       left=0.515, right=0.99, bottom=0.01, top=0.63)
fig = plt.figure(figsize=(7, 4.2))
xlims, ylims = [2.2, -2.2], [-2.2, 2.2]

# colorscales
vmin, vmax = 0, 50
rmin, rmax = -5, 5
smin, smax = -0.5, 0.5


# load the TRUTH cubes
hdut = fits.open(tfile_BLC+'.cube.image.fits')
tI_blc, h_blc = np.squeeze(hdut[0].data), hdut[0].header
hdut.close()
    
hdut = fits.open(tfile_WSU+'.cube.image.fits')
tI_wsu, h_wsu = np.squeeze(hdut[0].data), hdut[0].header
hdut.close()

# set the coordinate grids
dx = 3600 * h_blc['CDELT1'] * \
     (np.arange(h_blc['NAXIS1']) - (h_blc['CRPIX1'] - 1))
dy = 3600 * h_blc['CDELT2'] * \
     (np.arange(h_blc['NAXIS2']) - (h_blc['CRPIX2'] - 1))
ext = (dx.max(), dx.min(), dy.min(), dy.max())

# set the PSF dimensions
bmaj_blc, bmin_blc, bpa_blc = h_blc['BMAJ'], h_blc['BMIN'], h_blc['BPA']
bm_blc = (np.pi * bmaj_blc * bmin_blc / (4 * np.log(2))) * (np.pi / 180)**2

bmaj_wsu, bmin_wsu, bpa_wsu = h_wsu['BMAJ'], h_wsu['BMIN'], h_wsu['BPA']
bm_wsu = (np.pi * bmaj_wsu * bmin_wsu / (4 * np.log(2))) * (np.pi / 180)**2
    
# set the frequencies and velocities
nu_blc = h_blc['CRVAL3'] + h_blc['CDELT3'] * \
         (np.arange(h_blc['NAXIS3']) - (h_blc['CRPIX3'] - 1))
vel_blc = sc.c * (1 - nu_blc / nu0)

nu_wsu = h_wsu['CRVAL3'] + h_wsu['CDELT3'] * \
         (np.arange(h_wsu['NAXIS3']) - (h_wsu['CRPIX3'] - 1))
vel_wsu = sc.c * (1 - nu_wsu / nu0)

# plot reference (truth) channel maps
for i in range(nch):
    # desired cube-frame index
    j_blc = ix0_blc - dch_blc * i 
    j_wsu = ix0_wsu - dch_wsu * i

    # calculate brightness temperatures
    tTb_blc = (1e-26 * tI_blc[j_blc,:,:] / bm_blc) * \
              sc.c**2 / (2 * sc.k * nu_blc[j_blc]**2)
    tTb_wsu = (1e-26 * tI_wsu[j_wsu,:,:] / bm_wsu) * \
              sc.c**2 / (2 * sc.k * nu_wsu[j_wsu]**2)

    # make the channel map plots
    # - BLC
    ax = fig.add_subplot(gs_dat0[i])
    im = ax.imshow(tTb_blc, origin='lower', cmap=cmap, extent=ext, 
                   aspect='auto', vmin=vmin, vmax=vmax)
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_xticks([])
    ax.set_yticks([])
    if i == 0:
        ax.text(0.04, 0.88, '$\mathtt{reference}$', ha='left', va='center', 
                color='w', transform=ax.transAxes, fontsize=5)

        beam = Ellipse((xlims[0] + 0.1*np.diff(xlims),
                        -xlims[0] - 0.1*np.diff(xlims)), 
                       3600 * bmaj_blc, 3600 * bmin_blc, angle=90-bpa_blc)
        beam.set_facecolor('w')
        ax.add_artist(beam)

    v_blc = vel_blc[j_blc]
    if np.abs(v_blc) < 0.001: v_blc = 0.0
    if np.logical_or(np.sign(v_blc) == 1, np.sign(v_blc) == 0):
        pref = '+'
    else:
        pref = ''
    vstr_blc = pref+'%.2f' % (1e-3 * v_blc)
    ax.text(0.98, 0.05, vstr_blc, transform=ax.transAxes, ha='right',
            va='center', fontsize=5, color='w')

    # - WSU
    ax = fig.add_subplot(gs_dat1[i])
    im = ax.imshow(tTb_wsu, origin='lower', cmap=cmap, extent=ext, 
                   aspect='auto', vmin=vmin, vmax=vmax)
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_xticks([])
    ax.set_yticks([])
    if i == 0:
        ax.text(0.04, 0.88, '$\mathtt{reference}$', ha='left', va='center',
                color='w', transform=ax.transAxes, fontsize=5)

        beam = Ellipse((xlims[0] + 0.1*np.diff(xlims),
                        -xlims[0] - 0.1*np.diff(xlims)), 
                       3600 * bmaj_wsu, 3600 * bmin_wsu, angle=90-bpa_wsu)
        beam.set_facecolor('w')
        ax.add_artist(beam)

    v_wsu = vel_wsu[j_wsu]
    if np.abs(v_wsu) < 0.001: v_wsu = 0.0
    if np.logical_or(np.sign(v_wsu) == 1, np.sign(v_wsu) == 0):
        pref = '+'
    else:
        pref = ''
    vstr_wsu = pref+'%.2f' % (1e-3 * v_wsu)
    ax.text(0.98, 0.05, vstr_wsu, transform=ax.transAxes, ha='right',
            va='center', fontsize=5, color='w')


# iterate over interpolation methods
for ia in range(len(imeth)):
    # load the regridded cubes
    hdur = fits.open(dfile_BLC+'_'+imeth[ia]+'.cube.image.fits')
    rI_blc = np.squeeze(hdur[0].data)
    hdur.close()

    hdur = fits.open(dfile_WSU+'_'+imeth[ia]+'.cube.image.fits')
    rI_wsu = np.squeeze(hdur[0].data)
    hdur.close()

    for i in range(nch):
        # desired cube-frame index
        j_blc = ix0_blc - dch_blc * i
        j_wsu = ix0_wsu - dch_wsu * i

        # calculate brightness temperatures
        rTb_blc = (1e-26 * (rI_blc[j_blc,:,:] - tI_blc[j_blc,:,:]) / bm_blc) * \
                  sc.c**2 / (2 * sc.k * nu_blc[j_blc]**2)
        rTb_wsu = (1e-26 * (rI_wsu[j_wsu,:,:] - tI_wsu[j_wsu,:,:]) / bm_wsu) * \
                  sc.c**2 / (2 * sc.k * nu_wsu[j_wsu]**2)

        # make the channel map plots
        # - BLC
        ax = fig.add_subplot(gs_res0[ia,i])
        sim = ax.imshow(rTb_blc, origin='lower', cmap=mymap, extent=ext,
                        aspect='auto', vmin=smin, vmax=smax)
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        ax.set_xticks([])
        ax.set_yticks([])
        if i == 0:
            ax.text(0.04, 0.88, nms[ia], ha='left', va='center',
                    color=col_nms, transform=ax.transAxes, fontsize=5)

        # - WSU
        ax = fig.add_subplot(gs_res1[ia,i])
        rim = ax.imshow(rTb_wsu, origin='lower', cmap=mymap, extent=ext,
                        aspect='auto', vmin=rmin, vmax=rmax)
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        ax.set_xticks([])
        ax.set_yticks([])
        if i == 0:
            ax.text(0.04, 0.88, nms[ia], ha='left', va='center',
                    color=col_nms, transform=ax.transAxes, fontsize=5)


# colorbars

# WSU side
cbax = fig.add_axes([0.255, 0.915, 0.23, 0.018])
cb = Colorbar(ax=cbax, mappable=im, orientation='horizontal', 
              ticklocation='top', extend='both')
cbax.tick_params(which='both', labelsize=6, pad=0)
cb.set_label('$T_{\\rm b}$ (K)', fontsize=6, rotation=0, labelpad=3)

cbax = fig.add_axes([0.255, 0.64, 0.23, 0.018])
cb = Colorbar(ax=cbax, mappable=rim, orientation='horizontal',
              ticklocation='top', extend='both', ticks=[-6, -4, -2, 0, 2, 4, 6])
cbax.tick_params(which='both', labelsize=6, pad=0)
cb.set_label('residual $T_{\\rm b}$ (K)', fontsize=6, rotation=0, labelpad=3)


# BLC side
cbax = fig.add_axes([0.76, 0.915, 0.23, 0.018])
cb = Colorbar(ax=cbax, mappable=im, orientation='horizontal',
              ticklocation='top', extend='both')
cbax.tick_params(which='both', labelsize=6, pad=0)
cb.set_label('$T_{\\rm b}$ (K)', fontsize=6, rotation=0, labelpad=3)

cbax = fig.add_axes([0.76, 0.64, 0.23, 0.018])
cb = Colorbar(ax=cbax, mappable=sim, orientation='horizontal',
              ticklocation='top', extend='both')#, ticks=[-6, -4, -2, 0, 2, 4, 6])
cbax.tick_params(which='both', labelsize=6, pad=0)
cb.set_label('residual $T_{\\rm b}$ (K)', fontsize=6, rotation=0, labelpad=3)


lax = fig.add_axes([0, 0.92, 0.2, 0.8], alpha=0)
lax.text(0.07, 0.015, 'WSU, $pre$-averaged', ha='left', va='center',
         transform=lax.transAxes, fontsize=8.5, color='k')
lax.axis('off')

lax = fig.add_axes([0.505, 0.92, 0.2, 0.8], alpha=0)
lax.text(0.07, 0.015, 'WSU, $post$-averaged', ha='left', va='center',
         transform=lax.transAxes, fontsize=8.5, color='k')
lax.axis('off')



figo = 'figs/demo_residcubes_avgoop'
fig.savefig(figo+'.pdf')
print('figure at '+figo+'.pdf')
fig.clf()
