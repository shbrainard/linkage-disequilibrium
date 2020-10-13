# linkage-disequilibrium
R functions to visualize the r2 output from vcftools

## r2Decay.R 
* calculates the r2 between pairs of markers separated by increasing genetic distance (up to 1 Mb).
* fits a self-starting non-linear decay function to the average of each bin of marker pairs
* plot this decay function on a log-scale, along with intercepts with the lines y=0.1 and y=0.2

## r2HeatMap.R
* plots the r2 values between all pairs of markers on a selected chromosome
* VCF files should be thinned prior to plotting (example is thinned to 1 marker / 1 kb)

## r2SlidingWindow.R
* calculates the average r2 between a given marker, and the next nearest 100 SNPs
* plots the resulting r2 values as a Manhattan plot
