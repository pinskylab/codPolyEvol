# codPolyEvol

Analysis scripts and processed data used in Reid, Star, and Pinsky 2023 (in revision) to evaluate evidence for polygenic selection in Atlantic cod (*Gadus morhua*) using temporal genomic data.

Building upon previous work by [Pinsky et al. (2021)](https://doi.org/10.1073/pnas.2025453118), which used temporally-sampled whole genome data from cod to evaluate whether cod showed a loss of genetic diversity or evidence of adaptation at loci of large effect in response to heavy fishing pressure in the 20th century, this project applies a novel covariance-based method for detecting polygenic selection developed by [Buffalo and Coop (2020)](https://doi.org/10.1073/pnas.1919039117) to evaluate whether cod display evidence for parallel polygenic adaptation across the Atlantic. Simulations were also used to evaluate whether neutral processes (migration and drift) could replicate the covariance signal observed in the empirical data and whether background selection on shared deleterious variation or selection on shared quantitative trait loci (QTL) were plausible explanations for the patterns observed.

The scripts provided here can be used to replicate the analysis. The basic workflow for this project consists of:

1) Filtering the SNP dataset and calculating allele frequencies for each time point and location using vcftools and plink (link script).

2) Calculating allele frequency change and covariance in allele frequency change, bootstrapping over the whole genome and individual linkage groups to obtain confidence intervals, examining specific chromosomal subsets or regions (coding sequences / inversions), and plotting results using R (link R script). 

3) Performing neutral simulations and a simulation incorporating deleterious mutations and QTL in SLiM (link SLiM scripts).

4) Analyzing and visualizing simulated data using R (link sim R scripts).

Please feel free to get in touch with me if you have any questions!

```
br450[at]rutgers.edu
** or **
nerdbrained[at]gmail.com
```
