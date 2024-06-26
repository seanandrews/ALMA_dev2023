import os
import sys
import importlib
import numpy as np
import matplotlib.pyplot as plt
from interpolators import *
from scipy.interpolate import interp1d
from scipy.ndimage import convolve1d

# plotting style setups
_ = importlib.import_module('plot_setups')
plt.style.use(['default', '/home/sandrews/mpl_styles/nice_line.mplstyle'])

gen_noise = False



# Define a snippet of channels and upsampled sub-channels
nch, upsample, nedge = 31, 20, 5
ch = np.arange(0, nch, 1)
ch_ = np.interp(np.arange((nch - 1) * upsample + 1),
                np.arange(0, nch * upsample, upsample), ch)
if (len(ch_) % 2 == 1): ch_ = ch_[:-1]
nch_ = len(ch_)


if gen_noise:
    # Random (Gaussian) draws for oversampled spectra
    ndraws = 100000
    raw_draws = np.random.normal(0, 1, (nch_, ndraws))

    # Convolve draws with appropriate SRF, downsample, clip
    # - BLC
    xch = ch_ - np.mean(ch_)
    blc_SRF = 0.5 * np.sinc(xch) + 0.25 * (np.sinc(xch - 1) + np.sinc(xch + 1))
    blc_draws_ = convolve1d(np.sqrt(upsample * 8/3) * raw_draws, 
                            blc_SRF / np.sum(blc_SRF), axis=0, mode='nearest')
    blc_draws = blc_draws_[::upsample,:]
    blc_draws = blc_draws[nedge:-nedge,:]
    blc_ch = np.arange(0, blc_draws.shape[0], 1)

    # - WSU
    _ = np.load('/pool/asha0/SCIENCE/csalt/csalt/data/WSU_SRF.npz')
    wint = interp1d(_['chix'], _['srf'], kind='cubic', 
                    fill_value=0, bounds_error=False)
    wsu_SRF = wint(xch)
    wsu_draws_ = convolve1d(np.sqrt(upsample * 1.15) * raw_draws, 
                            wsu_SRF / np.sum(wsu_SRF), axis=0, mode='nearest')
    wsu_draws = wsu_draws_[::upsample,:]
    wsu_draws = wsu_draws[nedge:-nedge,:]
    wsu_ch = np.arange(0, wsu_draws.shape[0], 1)

    # save results
    np.savez('noise_data.npz', ndraws=ndraws, 
                               raw_ch=ch_, raw_draws=raw_draws, 
                               blc_ch=blc_ch, blc_draws=blc_draws, 
                               wsu_ch=wsu_ch, wsu_draws=wsu_draws) 
else:
    _ = np.load('../ALMA-regridding/noise_data.npz')
    ndraws = _['ndraws']
    blc_ch, blc_draws = _['blc_ch'], _['blc_draws']
    wsu_ch, wsu_draws = _['wsu_ch'], _['wsu_draws']



# sub-channel interpolates
subchan = 0.1
chi = np.arange(0, blc_ch[-1], subchan)
print(chi)

# nearest neighbor interpolation of the draws; "noise" of the interpolates
blc_nearest = np.std(interp_nearest(blc_ch, blc_draws, chi, axis=0), axis=1)
wsu_nearest = np.std(interp_nearest(wsu_ch, wsu_draws, chi, axis=0), axis=1)

# linear interpolation of the draws; "noise" of the interpolates
blc_linear = np.std(interp_linear(blc_ch, blc_draws, chi, axis=0), axis=1)
wsu_linear = np.std(interp_linear(wsu_ch, wsu_draws, chi, axis=0), axis=1)

# cubic interpolation of the draws; "noise" of the interpolates
blc_cubic = np.std(interp_cubic(blc_ch, blc_draws, chi, axis=0), axis=1)
wsu_cubic = np.std(interp_cubic(wsu_ch, wsu_draws, chi, axis=0), axis=1)

# fftshift interpolation of the draws; "noise" of the interpolates
blc_fftshift = np.zeros(chi.shape)
wsu_fftshift = np.zeros(chi.shape)
for i in range(int(1 / subchan)):
    dch = chi[i] - blc_ch[0]
    _ = np.std(np.array([interp_fftshift(blc_ch, blc_draws[:,j], dch) 
                         for j in range(ndraws)]).T, axis=1)
    e = np.std(np.array([interp_fftshift(wsu_ch, wsu_draws[:,j], dch)
                         for j in range(ndraws)]).T, axis=1)
    for j in range(len(_)-1):
        blc_fftshift[i + j * 10] = _[j]
        wsu_fftshift[i + j * 10] = e[j]




# make the figure
fig, ax = plt.subplots(nrows=4, figsize=(3.5, 3.5))
l, r, t, b = 0.12, 0.88, 0.99, 0.09

# shift to move away from edge effects
ax[0].plot(chi-nedge, blc_nearest, '--', color='#F55D3E')
ax[0].plot(chi-nedge, wsu_nearest, '-', color='#F55D3E')
ax[1].plot(chi-nedge, blc_linear, '--', color='#878E88')
ax[1].plot(chi-nedge, wsu_linear, '-', color='#878E88')
ax[2].plot(chi-nedge, blc_cubic, '--', color='#F7CB15')
ax[2].plot(chi-nedge, wsu_cubic, '-', color='#F7CB15')
ax[3].plot(chi-nedge, blc_fftshift, '--', color='#76BED0')
ax[3].plot(chi-nedge, wsu_fftshift, '-', color='#76BED0')

ax[0].plot([-10, -9], [-10, -9], '-', color='k', label='WSU')
ax[0].plot([-10, -9], [-10, -9], '--', color='k', label='BLC')

ax[0].legend(loc='lower right', prop={'size':7}, framealpha=0, 
             markerfirst=True)

nms = ['$\mathtt{nearest}$', '$\mathtt{linear}$', 
       '$\mathtt{cubic}$', '$\mathtt{fftshift}$']

[ax[j].text(0.03, 0.75, nms[j], transform=ax[j].transAxes, 
            fontsize=7, ha='left') for j in range(4)]

[ax[j].set_xlim([0, 5]) for j in range(4)]
[ax[j].set_ylim([0.5, 1.25]) for j in range(4)]
[ax[j].set_yticks([0.6, 0.8, 1.0, 1.2]) for j in range(4)]
[ax[j].set_xticklabels([]) for j in range(3)]
[ax[j].set_yticklabels([]) for j in range(3)]
ax[-1].set_xlabel('channel offsets in units of $\\Delta \\nu$')
ax[-1].set_ylabel('$\\sigma_{\\rm eff} \,\, / \,\, \\sigma_{\\rm true}$',
                  labelpad=6)

fig.subplots_adjust(left=l, right=r, bottom=b, top=t, hspace=0.07)
plt.savefig('figs/noise.pdf')
fig.clf()
