import os
import sys
sys.path.append('/pool/asha0/SCIENCE/csalt/')
import numpy as np
from csalt.model import *
from csalt.helpers import *
from scipy.ndimage import convolve1d
from scipy.interpolate import interp1d


# setup
sdir = '/data/sandrews/ALMA_regridding/storage/'

# over-sampled, high-resolution "truth" for this model
dnu_ = 1e3

# convolve the "truth" with these SRFs
convolve_SRF = True
SRFs = ['ALMA-WSU', 'ALMA-BLC']
dnu_srf = [13.5e3, 122e3]
nms = ['native', '122kHz']


#-------------------------------------------------------------------------------

# Instantiate a csalt model
cm = model('CSALT0', path='/pool/asha0/SCIENCE/csalt/')

# Create an empty MS from scratch
cdir = '/pool/asha0/casa-release-5.7.2-4.el7/data/alma/simmos/'
cm.template_MS(sdir+'templates/template_TRUTH.ms',
               config=[cdir+'alma.cycle8.6.cfg'],
               t_total=['15min'], t_integ='30s',
               observatory='ALMA', date=['2025/04/20'], HA_0=['0h'],
               restfreq=230.538e9, dnu_native=dnu_, V_tune=0.0e3, V_span=7.5e3,
               RA='16:00:00.00', DEC='-30:00:00.00')

# Get the data dictionary from the empty MS
ddict = read_MS(sdir+'templates/template_TRUTH.ms')

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

# Generate the TRUE (1 kHz-sampled) visibility spectra (no SRF convolution!)
fixed_kw = {'FOV': 5.11, 'Npix': 128, 'dist': 150, 'Nup': None,
            'doppcorr': 'approx', 'SRF': None, 'noise_inject': None}
sampl_mdict = cm.modeldict(ddict, pars, kwargs=fixed_kw)
write_MS(sampl_mdict, outfile=sdir+'TRUTH.ms')



if convolve_SRF:
    def SRF_kernel(dnu_native, SRF_type='None', dnu_sampled=1e3, N_over=25):
        # the over-sampling factor for input spectra
        f_oversample = dnu_native / dnu_sampled

        # number of over-sampled channels to sample N_over native channels
        nchan = int(np.round(N_over * f_oversample))

        # channel frequency grid
        chix = np.arange(nchan)
        xch = (chix - np.mean(chix)) / f_oversample

        # compute the SRF on the channel frequency grid
        if SRF_type == 'ALMA-BLC':
            srf = 0.50 * np.sinc(xch) + \
                  0.25 * (np.sinc(xch - 1) + np.sinc(xch + 1))
        elif SRF_type == 'ALMA-WSU':
            _wsu = np.load('/pool/asha0/SCIENCE/csalt/csalt/data/WSU_SRF.npz')
            wint = interp1d(_wsu['chix'], _wsu['srf'],
                            fill_value='extrapolate', kind='cubic')
            srf = wint(xch)
        else:
            print('I do not know that SRF.')
            return 0

        return srf / np.sum(srf)


    # Load the truth MS file into a dictionary
    t_dict = read_MS(sdir+'TRUTH.ms')

    for ik in range(len(SRFs)):
        # message
        print('\nConvolving "truth" with the '+SRFs[ik]+' response function.')

        # assign output MS filename
        truconv_MS = sdir+'TRUTH.'+SRFs[ik]+'_'+nms[ik]

        # copy the truth MS into an output SRF-convolved truth file
        os.system('rm -rf '+truconv_MS+'.ms')
        os.system('cp -r '+sdir+'TRUTH.ms '+truconv_MS+'.ms')
        ot_dict = read_MS(truconv_MS+'.ms')

        # compute the SRF kernel
        dnu_truth = np.diff(t_dict['0'].nu_TOPO)[0]
        ker = SRF_kernel(dnu_srf[ik], SRF_type=SRFs[ik], dnu_sampled=dnu_truth)

        # cycle over the timestamps and process the *true* visibility spectra
        for j in range(t_dict['0'].nstamps):
            print(str(j+1)+' / '+str(t_dict['0'].nstamps))

            # indices corresponding to this timestamp
            tixl = np.min(np.where(t_dict['0'].tstamp == j))
            tixh = np.max(np.where(t_dict['0'].tstamp == j)) + 1

            # convert true visibilities into a useable array format
            tvis = np.empty((t_dict['0'].npol, t_dict['0'].nchan,
                             tixh - tixl, 2))
            tvis[0,:,:,0] = t_dict['0'].vis[0,:,tixl:tixh].real
            tvis[1,:,:,0] = t_dict['0'].vis[1,:,tixl:tixh].real
            tvis[0,:,:,1] = t_dict['0'].vis[0,:,tixl:tixh].imag
            tvis[1,:,:,1] = t_dict['0'].vis[1,:,tixl:tixh].imag

            # convolve the true spectra with the SRF for this timestamp
            tvis_conv = convolve1d(tvis, ker, axis=1, mode='nearest')

            # revert back to complex format
            _tvis = tvis_conv[:,:,:,0] + 1j * tvis_conv[:,:,:,1]

            # assign to regridded TRUTH
            ot_dict['0'].vis[:,:,tixl:tixh] = _tvis

        # write out the MS file
        write_MS(ot_dict, outfile=truconv_MS+'.ms', direct_file=True)
        print('SRF-convolved truth is stored in '+truconv_MS+'.ms')
