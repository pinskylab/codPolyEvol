# codPolyEvol

Analysis scripts and processed data used in Reid, Star, and Pinsky 2023 (in revision) to evaluate evidence for polygenic selection in Atlantic cod (*Gadus morhua*) using temporal genomic data.

Building upon previous work by [Pinsky et al. (2021)](https://doi.org/10.1073/pnas.2025453118), which used temporally-sampled whole genome data from cod to evaluate whether cod showed a loss of genetic diversity or evidence of adaptation at loci of large effect in response to heavy fishing pressure in the 20th century, this project applies a novel covariance-based method for detecting polygenic selection developed by [Buffalo and Coop (2020)](https://doi.org/10.1073/pnas.1919039117) to evaluate whether cod display evidence for parallel polygenic adaptation across the Atlantic. Simulations were also used to evaluate whether neutral processes (migration and drift) could replicate the covariance signal observed in the empirical data and whether background selection on shared deleterious variation or selection on shared quantitative trait loci (QTL) were plausible explanations for the patterns observed.

[empirical_data](https://github.com/pinskylab/codPolyEvol/tree/main/empirical_data) contains allele frequency data from historic and contemporary cod populations and some quality control data.

[slim_code](https://github.com/pinskylab/codPolyEvol/tree/main/slim_code) contains code used to run simulations, and [slim_outputs] (https://github.com/pinskylab/codPolyEvol/tree/main/slim_outputs) contains summary outputs (allele frequencies and genetic divergences) from these simulations.

[analysis_scripts](https://github.com/pinskylab/codPolyEvol/tree/main/analysis_scripts) contains scripts that can be used to replicate the analysis. 

The basic workflow for this project consists of:

1) Filtering the SNP dataset and calculating allele frequencies for each time point and location using vcftools v.0.1.17 and plink v.2.0 using [these scripts](https://github.com/pinskylab/codPolyEvol/blob/main/analysis_scripts/filtering_afreq_scripts.txt).

2) Calculating allele frequency change and covariance in allele frequency change, bootstrapping over the whole genome and individual linkage groups to obtain confidence intervals, examining specific chromosomal subsets or regions (coding sequences / inversions), and plotting results using [R scripts](https://github.com/pinskylab/codPolyEvol/blob/main/analysis_scripts/freqchange_convcor.R) in R v.4.2.0. 

3) Performing neutral simulations and a simulation incorporating deleterious mutations and QTL in SLiM v.4.0 using [these scripts](https://github.com/pinskylab/codPolyEvol/tree/main/slim_code).

4) Analyzing and visualizing simulated data using R v.4.2.0 (separate scripts for [calculating covariance](https://github.com/pinskylab/codPolyEvol/blob/main/analysis_scripts/slim_covcalc_maf.R), [Fst values](https://github.com/pinskylab/codPolyEvol/blob/main/analysis_scripts/slim_fst_maf.R), and [mutation age vs prevalence](https://github.com/pinskylab/codPolyEvol/blob/main/analysis_scripts/slim_mutationages.R)).

Please feel free to get in touch with me if you have any questions!

```
br450[at]rutgers.edu
** or **
nerdbrained[at]gmail.com
```
