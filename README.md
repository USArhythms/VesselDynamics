# VesselDynamics
Code for processing widefield and two photon vascular imaging data
Many files in this repository use functions from the CHRONUX toolbox: http://chronux.org/ and/or functions developed by Xiang Ji: https://github.com/xiangjiph/VCRA

## Widefield
Extract and save df/f (x,t) from tiff images. Compute the peak vasomotor frequency, preform space-frequency SVD, and plot the dominant mode's magnitude and phase. 

## Pial2P
Process two photon frame scan tiff images. Perform background subtraction, create vessel masks, extract GCaMP df/f (x,t), and programatically draw cross-lines for FWHM calculations. Includes parallelized code to extract FWHM, utilizing many cross-lines to reduce noise in diameter estimation.

## Penetrating Ves
Process interleaved two photon frame scan tiff images. Estimate cross sectional diameter using the TiRS method (Developed by Dr. Patrick Drew), and perform spectral calculations on the resulting diameter time series. Calculate spectra, phase, and coherence for each trial and plot these results.

## Statistics
Helper functions for various calculations. Perform regression through the origin and calculate the associated uncertainty in slope. Perform line subtraction (developed by Partha Mitra and others) to detect and remove sinusoidal components from complex timeseries data. Calculate 2D cross correlation using MATLAB's fft2 to detect motion artifacts in two photon data. Plot segmentation of vessels on top of the associated mask or raw image for visualization purposes.

## Plotting
Code to produce figures. Load each dataset individually, plot and save individual graphs and figures across different imaging modalities.
