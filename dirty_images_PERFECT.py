import os
import sys
from casatasks import tclean, exportfits

# setup
sdir = '/data/sandrews/ALMA_regridding/storage/'

# example for demo WSU setup
name = 'ALMA-WSU_native'
binn = '.bin8x'

# example for demo BLC setup
#name = 'ALMA-BLC_122kHz'
#binn = ''

nm_f = ['04']

for ig in range(len(nm_f)):

    ### regridded TRUTH
    prefix = 'TRUTH_PERFECT.'+name+binn+'.REGRID'+nm_f[ig]

    # prepare / cleanup
    ext = ['image', 'model', 'pb', 'psf', 'residual', 'sumwt']
    for i in range(len(ext)):
        if os.path.exists(sdir+'images/'+prefix+'.cube.'+ext[i]):
            os.system('rm -rf '+sdir+'images/'+prefix+'.cube.'+ext[i])

    # dirty cube
    tclean(vis=sdir+prefix+'.ms', 
           imagename=sdir+'images/'+prefix+'.cube',
           specmode='cubedata', imsize=512, cell='0.02arcsec', niter=0,
           restoringbeam='common')

    # make a fits cube
    exportfits(sdir+'images/'+prefix+'.cube.image', 
               fitsimage=sdir+'images/'+prefix+'.cube.image.fits',
               overwrite=True)
