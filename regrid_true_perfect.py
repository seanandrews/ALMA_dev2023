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

# example for BLC
#name = 'ALMA-BLC_122kHz'
#binn = ''


# input (un-regridded) data filename
in_MS = sdir+name+'_SAMPLED'+binn

# corresponding truth filename
t_MS = sdir+'TRUTH'

# load the truth
ti_dict = read_MS(t_MS+'.ms')

# offsets (and labeling strings) in native channel units 
nm_f = ['04']


# iterate over interpolation methods and interpolate offsets
for ig in range(len(nm_f)):

    ### "truth"
    # copy a regridded MS into a new file for the truth
    tout_MS = t_MS+'_PERFECT.'+name+binn+'.REGRID'+nm_f[ig]
    os.system('rm -rf '+tout_MS+'.ms*')
    os.system('cp -r '+in_MS+'.REGRID'+nm_f[ig]+'_nearest''.ms '+tout_MS+'.ms')
    
    # load in the regridded MS shell
    to_dict = read_MS(tout_MS+'.ms')

    # cycle over timestamps to interpolate onto regridded LSRK channels
    for j in range(ti_dict['0'].nstamps):
        ixl = np.min(np.where(ti_dict['0'].tstamp == j))
        ixh = np.max(np.where(ti_dict['0'].tstamp == j)) + 1

        _ = interp_cubic(ti_dict['0'].nu_LSRK[j,:], 
                         ti_dict['0'].vis[:,:,ixl:ixh], 
                         to_dict['0'].nu_TOPO, axis=1)
        to_dict['0'].vis[:,:,ixl:ixh] = _

    # write out the updated MS file
    write_MS(to_dict, outfile=tout_MS+'.ms', direct_file=True)
