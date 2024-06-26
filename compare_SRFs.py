import os
import sys
import importlib
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

# plotting style setups
_ = importlib.import_module('plot_setups')
plt.style.use(['default', '/home/sandrews/mpl_styles/nice_line.mplstyle'])


# load tabulated WSU SRF
_ = np.load('/pool/asha0/SCIENCE/csalt/csalt/data/WSU_SRF.npz')
chix, SRF_wsu = _['chix'], _['srf']

# compute BLC SRF on these channel indices
SRF_blc = 0.5 * np.sinc(chix) + \
          0.25 * np.sinc(chix - 1) + 0.25 * np.sinc(chix + 1)


# plot both
fig, ax = plt.subplots(figsize=(3.5, 2.1))

ax.plot(chix, SRF_wsu / np.trapz(SRF_wsu, chix), '-k', label='WSU')
ax.plot(chix, SRF_blc / np.trapz(SRF_blc, chix), '--', color='xkcd:gray', label='BLC')


# labeling
ax.set_xlim([-8, 8])
ax.set_xlabel('channel offsets in units of $\\Delta \\nu$')
ax.set_ylim([-0.08, 1.08])
ax.set_ylabel(r'$\Phi(\nu)$')

ax.legend(prop={'size': 7}, edgecolor='w', markerfirst=False)


fig.subplots_adjust(left=0.115, right=0.885, bottom=0.15, top=0.98)


plt.savefig('figs/compare_SRFs.pdf')
fig.clf()
