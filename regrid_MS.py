import os
import sys
sys.path.append('/pool/asha0/SCIENCE/csalt/')
from csalt.model import *
from csalt.helpers import *
import numpy as np
import scipy.constants as sc
from casatasks import mstransform
from interpolators import *


# input data file (not regridded)
sdir = '/data/sandrews/ALMA_regridding/storage/'
in_MS = sdir+'ALMA-WSU_native_SAMPLED.bin8x'

# interpolation methods
imeth = ['nearest', 'linear', 'cubic', 'fftshift']

# information on input channel grid
i_dict = read_MS(in_MS+'.ms')
_nu = i_dict['0'].nu_LSRK[round(i_dict['0'].nstamps / 2),:] 
restfreq = 230.538e9

# offsets (and labeling strings) in native channel units 
df = [0.5]        
nm_f = ['05']


# iterate over interpolation methods and interpolate offsets
for im in range(len(imeth)):

    for ig in range(len(df)):  

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
