First we need to decide the genomic region to capture in the haplotype clustering analysis. 
We do this by using a modelling approach that fits an exponential model to find the centre of the H12 peak
From this we decide on the % decay from the centre. In my analysis i used 5% decay on each side to get the start and end position
To model the centre of the peaks you first run the [./haplotype_clustering/peak_centre/h12_calibration.ipynb] notebook which generate h12 GWSS calibration data for different window sizes
For each chromosome arm, we can now determine the centre of the H12 peaks using the h12_signal_detection notebooks. 
Peak Untils contains functions required by the signal detection notebooks. 

Now we have the genomic boundary positions. We can plot the dendograms

# [./plot_haplotypes_clusters.ipynb] A notebook to plot dendogram of hapotypes within specified genomic regions. It adds SNPS and Dups that tags/associated with each of the haplotype clusters
# To find snps associated with each haplotypes, you need to save the output of the  find_clusters function in this notebook, download snps within a region of interest using 
#haplotype_snps_genotype.ipynb and use cluster_snps_association.r to find the snps that are strongly associated with each haplotype clusters. 