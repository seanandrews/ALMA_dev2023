import os
import sys
from casatasks import tclean, exportfits

# naming
nm_f = ['04']
imeth = ['nearest', 'linear', 'cubic', 'fftshift']


for ig in range(len(nm_f)):

    for im in range(len(imeth)):
         
        ### regridded data
        prefix = 'HD_163296_220GHz_spw6.bin_30s.2x_spw13.REGRID'+\
                 nm_f[ig]+'_'+imeth[im]

        # prepare / cleanup
        ext = ['image', 'model', 'pb', 'psf', 'residual', 'sumwt']
        for i in range(len(ext)):
            if os.path.exists('MAPS_data/images/'+prefix+'.cube.'+ext[i]):
                os.system('rm -rf MAPS_data/images/'+prefix+'.cube.'+ext[i])

        # dirty cube
        tclean(vis='MAPS_data/'+prefix+'.ms',
               imagename='MAPS_data/images/'+prefix+'.cube',
               specmode='cubedata', imsize=512, cell='0.02arcsec', 
               weighting='natural', niter=0)

        # make a fits cube
        exportfits('MAPS_data/images/'+prefix+'.cube.image',
                   fitsimage='MAPS_data/images/'+prefix+'.cube.image.fits',
                   overwrite=True)
