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
name = 'ALMA-WSU_native'
binn = '.bin8x'
binfactor = 8

# example for BLC
#name = 'ALMA-BLC_122kHz'
#binn = ''
#binfactor = 1


# input (un-regridded) data filename
in_MS = sdir+name+'_SAMPLED'+binn

# corresponding truth filename
t_MS = sdir+'TRUTH.'+name

# interpolation methods
imeth = ['nearest', 'linear', 'cubic', 'fftshift']

# load input data and information on input channel grid
i_dict = read_MS(in_MS+'.ms')
_nu = i_dict['0'].nu_LSRK[round(i_dict['0'].nstamps / 2),:] 
restfreq = 230.538e9

# load the truth
ti_dict = read_MS(t_MS+'.ms')

# offsets (and labeling strings) in native channel units 
df = [0.49]
nm_f = ['049']


# iterate over interpolation methods and interpolate offsets
for ig in range(len(df)):

    ### "data" 
    for im in range(len(imeth)):

        # channelization instructions
        chanstart = f'{_nu[0] + df[ig] * np.mean(np.diff(_nu)): .3f}'+'Hz'
        chanwidth = f'{np.mean(np.diff(_nu)): .3f}'+'Hz'

        # output data file
        out_MS = in_MS+'.REGRID'+nm_f[ig]+'_'+imeth[im]

        # mstransform to get an MS file on desired output grid
        os.system('rm -rf '+out_MS+'.ms*')
        mstransform(vis=in_MS+'.ms', outputvis=out_MS+'.ms',
                    datacolumn='data', regridms=True, mode='frequency',
                    nchan=len(_nu), start=chanstart, width=chanwidth, 
                    restfreq=str(restfreq / 1e9)+'GHz', outframe='LSRK')

        # read in the regridded MS for this EB into a dictionary
        o_dict = read_MS(out_MS+'.ms')

        # cycle over timestamps to interpolate onto regridded LSRK channels
        for j in range(i_dict['0'].nstamps):
            ixl = np.min(np.where(i_dict['0'].tstamp == j))
            ixh = np.max(np.where(i_dict['0'].tstamp == j)) + 1

            if imeth[im] == 'fftshift':
                y = np.transpose(i_dict['0'].vis[:,:,ixl:ixh], (1, 0, 2))
                x = i_dict['0'].nu_LSRK[j,:]
                chshift = np.mean(x - o_dict['0'].nu_LSRK[j,:])
                ysample = eval("interp_"+imeth[im]+"(x, y, chshift)")
                o_dict['0'].vis[:,:,ixl:ixh] = np.transpose(ysample, (1, 0, 2))
            else:
                cmd = "interp_"+imeth[im]+\
                      "(i_dict['0'].nu_LSRK[j,:], "+\
                      "i_dict['0'].vis[:,:,ixl:ixh], "+\
                      "o_dict['0'].nu_TOPO, axis=1)"
                o_dict['0'].vis[:,:,ixl:ixh] = eval(cmd)

        # write out the MS file
        write_MS(o_dict, outfile=out_MS+'.ms', direct_file=True)


    ### "truth"
    # copy a regridded MS into a new file for the truth
    tout_MS = t_MS+binn+'.REGRID'+nm_f[ig]
    os.system('rm -rf '+tout_MS+'.ms*')
    os.system('cp -r '+out_MS+'.ms '+tout_MS+'.ms')
    
    # load in the regridded MS shell
    to_dict = read_MS(tout_MS+'.ms')

    # cycle over timestamps to interpolate onto regridded LSRK channels
    for j in range(ti_dict['0'].nstamps):
        ixl = np.min(np.where(ti_dict['0'].tstamp == j))
        ixh = np.max(np.where(ti_dict['0'].tstamp == j)) + 1

        # interpolate if no binning
        if binfactor > 1:
            tbf = int(round(np.mean(np.diff(to_dict['0'].nu_LSRK[j,:])) / \
                            np.mean(np.diff(ti_dict['0'].nu_LSRK[j,:]))))

            # find the right offset for the pre-binned frequencies
            _nu = ti_dict['0'].nu_LSRK[j,:]
            _bnu = _nu[:(_nu.size // tbf) * tbf].reshape(-1, tbf).mean(axis=1)
            doff = np.mean(_bnu - to_dict['0'].nu_LSRK[j,:])

            # interpolate true visibilities onto pre-binned frequencies
            tinterp = interp1d(ti_dict['0'].nu_LSRK[j,:],
                               ti_dict['0'].vis[:,:,ixl:ixh],
                               kind='cubic', fill_value='extrapolate', axis=1)
            svis = tinterp(ti_dict['0'].nu_LSRK[j,:] - doff)

            # now bin the visibilities in the spectral dimension
            t_ = np.transpose(svis, (1, 0, 2))
            _ = [np.take(t_, np.arange(i*tbf, (i+1)*tbf), 0).mean(axis=0) \
                 for i in np.arange(svis.shape[1] // tbf)]
            to_dict['0'].vis[:,:,ixl:ixh] = np.array(_).transpose((1, 0, 2))

        else:
            _ = interp_cubic(ti_dict['0'].nu_LSRK[j,:], 
                             ti_dict['0'].vis[:,:,ixl:ixh], 
                             to_dict['0'].nu_TOPO, axis=1)
            to_dict['0'].vis[:,:,ixl:ixh] = _

    # write out the updated MS file
    write_MS(to_dict, outfile=tout_MS+'.ms', direct_file=True)
