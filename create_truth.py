import os
import sys
sys.path.append('/pool/asha0/SCIENCE/csalt/')
import numpy as np
from csalt.model import *
from csalt.helpers import *


# setup
name = 'TRUTH'
dnu_ = 1e3


# Instantiate a csalt model
cm = model('CSALT0', path='/pool/asha0/SCIENCE/csalt/')

# Create an empty MS from scratch
sdir = '/data/sandrews/ALMA_regridding/storage/'
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

# Generate the TRUE (2 kHz-sampled) visibility spectra
fixed_kw = {'FOV': 5.11, 'Npix': 128, 'dist': 150, 'Nup': None,
            'doppcorr': 'approx', 'SRF': None, 'noise_inject': None}
sampl_mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)
write_MS(sampl_mdict, outfile=sdir+name+'.ms')
