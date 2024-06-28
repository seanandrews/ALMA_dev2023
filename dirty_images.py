import os
import sys
from casatasks import tclean, exportfits

# setup
sdir = '/data/sandrews/ALMA_regridding/storage/'

# example for demo WSU setup
name = 'ALMA-WSU_native'
binn = '.bin8x'

# example for demo WSU post-averaging setup
name = 'ALMA-WSU_native'
binn = '_postavg.bin8x'

# example for demo BLC setup
#name = 'ALMA-BLC_122kHz'
#binn = ''

img_truth = True
nm_f = ['049']
imeth = ['nearest', 'linear', 'cubic', 'fftshift']


for ig in range(len(nm_f)):

    if img_truth:
        ### regridded TRUTH
        prefix = 'TRUTH.'+name+binn+'.REGRID'+nm_f[ig]

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


    # iterate over interpolation methods
    for im in range(len(imeth)):
        
        ### regridded data
        prefix = name+'_SAMPLED'+binn+'.REGRID'+nm_f[ig]+'_'+imeth[im]

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
