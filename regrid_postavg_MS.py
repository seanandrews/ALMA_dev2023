import os
import sys
sys.path.append('/pool/asha0/SCIENCE/csalt/')
from csalt.model import *
from csalt.helpers import *
import numpy as np
import scipy.constants as sc
from casatasks import mstransform
from interpolators import *


# storage data directory
sdir = '/data/sandrews/ALMA_regridding/storage/'

# example for binned WSU
name = sdir+'ALMA-WSU_native_SAMPLED'
binn = '.bin8x'
binfactor = 8


# load input data and information on input channel grid
i_dict = read_MS(name+'.ms')
_nu = i_dict['0'].nu_LSRK[round(i_dict['0'].nstamps / 2),:] 
restfreq = 230.538e9


# interpolation methods
imeth = ['nearest', 'linear', 'cubic', 'fftshift']

# offsets (and labeling strings) in native channel units 
df = [0.4]        
nm_f = ['04']


# iterate over interpolation methods and interpolate offsets
for ig in range(len(df)):

    for im in range(len(imeth)):

        # make a copy of the 'target' (output) MS structure
        out_MS = name+'_postavg'+binn+'.REGRID'+nm_f[ig]+'_'+imeth[im]
        os.system('rm -rf '+out_MS+'.ms*')
        os.system('cp -r '+name+binn+'.REGRID'+nm_f[ig]+'_'+imeth[im]+'.ms '+\
                  out_MS+'.ms')

        # load in the regridded MS shell
        o_dict = read_MS(out_MS+'.ms')

        # cycle over timestamps to interpolate onto regridded LSRK channels
        for j in range(i_dict['0'].nstamps):
            ixl = np.min(np.where(i_dict['0'].tstamp == j))
            ixh = np.max(np.where(i_dict['0'].tstamp == j)) + 1

            # bin factor
            tbf = int(round(np.mean(np.diff(o_dict['0'].nu_LSRK[j,:])) / \
                            np.mean(np.diff(i_dict['0'].nu_LSRK[j,:]))))

            # find the right offset for the pre-binned frequencies
            _nu = i_dict['0'].nu_LSRK[j,:]
            _bnu = _nu[:(_nu.size // tbf) * tbf].reshape(-1, tbf).mean(axis=1)
            doff = np.mean(_bnu - o_dict['0'].nu_LSRK[j,:])

            # interpolate true visibilities onto pre-binned frequencies
            tinterp = interp1d(i_dict['0'].nu_LSRK[j,:],
                               i_dict['0'].vis[:,:,ixl:ixh],
                               kind='cubic', fill_value='extrapolate', axis=1)
            svis = tinterp(i_dict['0'].nu_LSRK[j,:] - doff)

            # now bin the visibilities in the spectral dimension
            t_ = np.transpose(svis, (1, 0, 2))
            _ = [np.take(t_, np.arange(i*tbf, (i+1)*tbf), 0).mean(axis=0) \
                 for i in np.arange(svis.shape[1] // tbf)]
            o_dict['0'].vis[:,:,ixl:ixh] = np.array(_).transpose((1, 0, 2))

        # write out the updated MS file
        write_MS(o_dict, outfile=out_MS+'.ms', direct_file=True)
