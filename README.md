# ALMA_dev2023
repository for 2023 ALMA Development Study on spectral regridding

create_data.py - Generates disk simulations using csalt package.

create_truth.py - Generates oversampled (1 kHz) true input simulation.

interpolators.py - Functions/wrappers for interpolating visibility spectra.

regrid_MS.py - Regrids a MS (roundtripping outside CASA).

regrid_postavg_MS.py - Regrids a MS with post-averaging (outside CASA).

regrid_true_perfect.py - Matches oversampled inputs onto regridded channels.

regrid_HD163.py - Regrids a MAPS dataset (outside CASA).

dirty_images.py - Makes dirty channel maps.

dirty_images_PERFECT.py - Makes the reference channel maps for Appendix.

dirty_images_HD163.py - Makes the dirty channel maps for MAPS dataset.
