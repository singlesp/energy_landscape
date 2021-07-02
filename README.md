# energy_landscape
Code to reproduce analysis in Singleton et al. 2021 ("LSD flattens the brain's energy landscape: insights from receptor-informed network control theory", BioRxiv.).

These scripts require that you have first downloaded the brain_states repo by Cornblath et al 2020 and added these to that directory (https://github.com/ejcorn/brain_states). Some of the code here is slightly modified Cornblath code, or at least relys heavily on that repo's core functions. 

The most simple way to re-implement this work is to start by looking at the script `example.m` in the Cornblath repo before digging into the full set of scripts.

## Requirements:
  - MATLAB R2017a or later
  - R 3.2.5 or later, with packages:
    - ggplot2
    - R.matlab
    - RColorBrewer
    - lm.beta
    - reshape2
    - viridis
    - plotrix
  
  - Friendship and positivity
  
## PRO-TIP: 
- k-means appears to have a bug when handling singles as data-type where it will have difficulty converging compared to the same data converted to a double. Use doubles for your TxnParc time-series data.
  


## General procedure and order of operations:

Please refer to the Cornblath repo README for a detailed overview of the vast majority of the code used for this analysis and how to use the .sh files provided in that analysis.

Download the Cornblath repo (https://github.com/ejcorn/brain_states) and add the scripts from this repo to the BASEDIR. (or however you'd like to get it all in the same place). Total install time <5 min.

For larger projects, you may need cluster access - with 30x 15 min fMRI scans, I was able to perform all of this analysis locally on a 2020 MacBook Pro with a total run time of <1 day. The permutation test requires the most time (12 hrs), followed by replications of k-means (1hr).

You may start by specifying a range of k over which to perform your preliminary analysis. At least 3 independent analyses have now found this ideal range to be 4-6 but you may find something different. 
To check - run `repeatkmeans_sps.m` over a range of k [2:max] (see Cornblath et al 2020 for choosing max k) followed by `elbow_sps.m` to view the variance explained plot. 

After a range or particular value of k is chosen, use `ami_calc.m` to assess clustering stability and choose the partition with the highest amount of adjusted mutual information shared with all other partitions.

You may find you wans to reorder your clusters for simplicity of comparing across choices of k or processing/data streams. To do this you can manually change parition values, centroid order, and cluster names with `comb_clusters.m`

Next, `transProbsEnergy.m` will generate figure 4a, (i) from manuscript.

`countclusters.m` will produce the violin data used in figure 3b, while `violinplots_ICLdata.R` will produce the figures themselves.

`subcentroids.m` will generate subject-specific centroids for energy calculations and compute subj-specific energies for a T you provide. This is a bit backwards, but after saving the centroids you can use `T_sweep_sps.m` to find the optimal T.

`permutePLavg.m` will scramble a receptor input vector over `nperms` and recalulated weighted energies

`E_corrs.m` will compute correlations between overall energy reduction and other values

## Input specification clarification

The user can also specify certain parameters in finalmain.sh:

- `split`: This variable was used in this analysis for simply keeping track of output files from different data/processing streams, whereas it is used for something else in Cornblath .sh files. 


## Sample data
To simulate the analysis you can replace <concTS> and the <connectivity> variables with random numbers. The BOLD data used for this analysis is available at: https://openneuro.org/datasets/ds003059/versions/1.0.0


Please contact Parker Singleton (sps253@cornell.edu) with any questions regarding this code.
