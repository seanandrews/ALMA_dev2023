import os
import sys
sys.path.append('/pool/asha0/SCIENCE/csalt/')
import numpy as np
from csalt.model import *
from csalt.helpers import *
from casatasks import split
import matplotlib as mpl
mpl.rcParams['backend'] = 'TkAgg'


# setup
sdir = '/data/sandrews/ALMA_regridding/storage/'

# example for demo WSU setup
name = 'ALMA-WSU_native'    # base-name for files (placed in sdir)
SRF_ = 'ALMA-WSU'           # type of SRF (ALMA or ALMA-WSU)
dnu_ = 13.5e3               # native channel spacing (in Hz)
online_bin = True           # emulate on-line (pre-)averaging?
binfactor = 8               # if True, binning factor (# of native channels)
                            # note that files will have a '.bin8x' suffix if,
                            # e.g., binfactor = 8

# example for demo BLC setup
name = 'ALMA-BLC_122kHz'    # base-name for files (placed in sdir)
SRF_ = 'ALMA'               # type of SRF (ALMA or ALMA-WSU)
dnu_ = 122e3                # native channel spacing (in Hz)
online_bin = False          # emulate on-line (pre-)averaging?
binfactor = 8               # if True, binning factor (# of native channels)
                            # note that files will have a '.bin8x' suffix if,
                            # e.g., binfactor = 8


#------------------------------------------------------------------------------

# Instantiate a csalt model
cm = model('CSALT0', path='/pool/asha0/SCIENCE/csalt/')

# Create an empty MS from scratch
cdir = '/pool/asha0/casa-release-5.7.2-4.el7/data/alma/simmos/'
cm.template_MS(sdir+'templates/template_'+name+'.ms',
               config=[cdir+'alma.cycle8.6.cfg'],
               t_total=['15min'], t_integ='30s',
               observatory='ALMA', date=['2025/04/20'], HA_0=['0h'],
               restfreq=230.538e9, dnu_native=dnu_, V_tune=0.0e3, V_span=7.5e3,
               RA='16:00:00.00', DEC='-30:00:00.00')

# Get the data dictionary from the empty MS
ddict = read_MS(sdir+'templates/template_'+name+'.ms')

# Set the CSALT model parameters
pars = np.array([
                   30,  # incl (deg)
                   60,  # PA (deg), E of N to redshifted major axis
                  1.0,  # Mstar (Msun)
                  300,  # R_out (au)
                  0.3,  # emission height z_0 (") at r = 1"
                  1.0,  # phi for z(r) = z_0 * (r / 1)**phi
                  150,  # Tb_0 at r = 10 au (K)
                 -0.5,  # q for Tb(r) = Tb_0 * (r / 10 au)**q
                   20,  # maximum Tb for back surface of disk (~ Tfreezeout)
                  297,  # linewidth at r = 10 au (m/s)
                  3.0,  # log(tau_0) at 10 au
                   -1,  # p for tau(r) = tau_0 * (r / 10 au)**p
                  0e3,  # systemic velocity (m/s)
                    0,  # RA offset (")
                    0   # DEC offset (")
                     ])

# Scale the noise (based on online ALMA sensitivity calculator)
if SRF_ == 'ALMA-WSU':
    x_scale = 1.2
else:
    x_scale = 1.0
noise = 7.3 * np.sqrt(30.5e3 / dnu_) / x_scale

# Generate the SAMPLED (noiseless) and NOISY visibility spectra
fixed_kw = {'FOV': 5.11, 'Npix': 128, 'dist': 150, 'Nup': 10, 
            'doppcorr': 'exact', 'SRF': SRF_, 'noise_inject': noise}
sampl_mdict, noisy_mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)
write_MS(sampl_mdict, outfile=sdir+name+'_SAMPLED.ms')
write_MS(noisy_mdict, outfile=sdir+name+'_NOISY.ms')

# on-line binning
if online_bin:
    os.system('rm -rf '+sdir+name+'_SAMPLED.bin'+str(binfactor)+'x.ms*')
    split(vis=sdir+name+'_SAMPLED.ms',
          outputvis=sdir+name+'_SAMPLED.bin'+str(binfactor)+'x.ms',
          datacolumn='data', width=binfactor)

    os.system('rm -rf '+sdir+name+'_NOISY.bin'+str(binfactor)+'x.ms*')
    split(vis=sdir+name+'_NOISY.ms',
          outputvis=sdir+name+'_NOISY.bin'+str(binfactor)+'x.ms',
          datacolumn='data', width=binfactor)
