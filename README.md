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

compare_SRFs.py - Makes Figure 1.

demo_SSP.py - Makes Figure 2.

noise.py - Makes Figure 3.

demo_visregrid_top.py - Makes top two rows of Figure 4.

demo_visregrid_bot.py - Makes bottom row of Figure 4.

demo_residcubes.py - Makes Figure 5.

demo_residcubes_avgoop.py - Makes Figure 6.

demo_residHD163.py - Makes Figure 7.

interp_timing_tests.py - Makes Figure 8.

interp_spike_tests.py - Makes Figure 9.

demo_visregrid_PERFECT_top.py - Makes top two rows of Figure 10 in Appendix.

demo_visregrid_PERFECT_bot.py - Makes bottom row of Figure 10 in Appendix.

demo_residcubes_PERFECT.py - Makes Figure 11 in Appendix.
