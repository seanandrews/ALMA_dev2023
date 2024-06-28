import os
import sys
sys.path.append('/pool/asha0/SCIENCE/csalt/')
from csalt.model import *
from csalt.helpers import *
import numpy as np
import scipy.constants as sc
import casatools
from casatasks import (split, mstransform, concat)
from interpolators import *


# input data file prefix
in_MS = 'MAPS_data/HD_163296_220GHz_spw6.bin_30s.2x_spw13'

# design output channel grid
i_dict = read_MS(in_MS+'.ms')
_nu = i_dict['2'].nu_LSRK[round(i_dict['2'].nstamps / 2), :]
restfreq = 230.538e9
ig, Nshifts = 4, 10

# interpolation methods
imeth = ['nearest', 'linear', 'cubic', 'fftshift']

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

        # cycle over EBs + stamps to interpolate onto regridded LSRK channels
        for i in range(o_dict['Nobs']):
            for j in range(i_dict[str(i)].nstamps):
                ixl = np.min(np.where(i_dict[str(i)].tstamp == j))
                ixh = np.max(np.where(i_dict[str(i)].tstamp == j)) + 1

                if imeth[im] == 'fftshift':
                    y = np.transpose(i_dict[str(i)].vis[:,:,ixl:ixh], (1,0,2))
                    x = i_dict[str(i)].nu_LSRK[j,:]
                    chshift = np.mean(x - o_dict[str(i)].nu_LSRK[j,:])
                    ys = eval("interp_"+imeth[im]+"(x, y, chshift)")
                    o_dict[str(i)].vis[:,:,ixl:ixh] = np.transpose(ys, (1,0,2))
            else:
                cmd = "interp_"+imeth[im]+\
                      "(i_dict[str(i)].nu_LSRK[j,:], "+\
                      "i_dict[str(i)].vis[:,:,ixl:ixh], "+\
                      "o_dict[str(i)].nu_TOPO, axis=1)"
                o_dict[str(i)].vis[:,:,ixl:ixh] = eval(cmd)

        # write out the MS file
        write_MS(o_dict, outfile=out_MS+'.ms', direct_file=True)
