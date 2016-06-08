Performs enrichment plots of GWAS p-value lists in DHSs (like Fig. 5 of Maurano, Humbert, et al. Science 2012)

The two key scripts are in src. The pipeline is divided into a Perl script which processes the overlap, and then a second script which does the plotting in R.
1) You can run the overlap as follows:
cd src
perl ./gwasVsRegions.pl -p pvalues.hg19.bed5 -s ../hg19/namedFDR5pctHotspots.starch -r ../results_hotspots_nocoding

The P-value file is provided by the user. It includes the dbSNP ID in  column 4, and the P-value in column 5; only the latter is actually used.

namedFDR5pctHotspots.starch is a starch archive (see BEDOPS) containing the DHS master list.

The intermediate results go into results_hotspots_nocoding/crohns. There is one text file per cell type, though you can see at the end of a Perl script I put all of these together into a single file before plotting.

Samples listed in excluded_samples.txt will be ignored

(This script was written by Eric Haugen, UW)


2) doGWASPlot.R is the plotting script. It is invoked automatically at the end of the overlap script. You can see that right now only the x-axis upper limit and of the number of cell types to label are parameterized on the command line. If you look inside, you'll see that the legendSamples list maps samples to group names and colors using regexp.


This could be easily parallelized by chromosome.

Key variables for each plot need to be optimized by the user:
* which samples to display on the plot
* what the labels to give each sample and how to color them (only the top N samples get labels)
* The sample groups and labels (i.e. endothelia, muscle, fetal, etc.)
* adjusting the range of the X and Y axes
* plot title
