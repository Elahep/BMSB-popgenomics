# Genome-wide scan analysis to detect outlier SNPs associated with invasiveness

Here we will use the recently developed C2 statistics <a href="https://academic.oup.com/mbe/article/37/8/2369/5821433" title="(Olazcuaga et al 2020)">(Olazcuaga et al 2020)</a> implemented in <a href="http://www1.montpellier.inra.fr/CBGP/software/baypass/" title="BayPass v2.3">BayPass v2.3</a> to identify outlier SNPs potentially associated with invasion success in the BMSB. 

***

## summary of the approach
For analysing loci associated with a binary trait such as invasive vs. non-invasive, we will use the contrast analysis in BayPass v2.3.
For such an analysis, the new version of BayPass has developed a nonparametric counterpart for the association model implemented in older versions of BayPass.
This new nonparametric model relies on a contrast statistic, C2, that compares the standardized population allele frequencies (i.e., allele frequencies corrected for the
population structure) between the two groups of populations specified by the binary covariable of interest. However, previous versions of BayPass relied on a paramteric models to calculate Bayes Factor (BF). Using the core model of Baypass, we can simultaneously estimate both BF and C2.
 
 
## 1-preparing input files

Using the <a href="https://gitlab.com/YDorant/Toolbox/-/blob/master/reshaper_baypass.py" title="reshaper_baypass.py">reshaper_baypass.py</a> script by <a href="https://gitlab.com/YDorant/Toolbox" title="Yann Dorant">Yann Dorant</a> we can convert VCF to the appropriate BayPass genotype file:

`python ./H1_bialSNP_MAF_geno_LD_reordered.vcf popmap.txt BMSB.geno`

The number of rows in the .geno file is equal to the number of SNPs and the number of columns is twice the number of populations (because we have two alleles for each site). The first number in each pairs of columns represents the number of the copies of the first allele in that population and the second number represents the number of the copies of the alternative allele. 0 indicates no copies or missing alleles.

We will also need a "contrast file" to assign each population to a binary trait, namely "invasive" and "native" in our analysis. The contrast file has one line and n columns (space seperated), where n denotes the number of sampled populations. In this file, we will allocate native populations to the first group (1) and invasive populations to the second/alternative group (-1). We can also exclude any population from pairwise comparisons by allocating them to the 0 group.
For the BMSB dataset, we have previously reordered individuals in the finel VCF using VCFtools so that we first have our native populations (Japan and China) which are followed by all invasive ones (10 invasive populations):

1 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1


## 2- running BayPass

On nesi, we need to load ifort.

```
 module load ifort
 ./i_baypass -gfile BMSB.geno -contrastfile contrast.ecotype -outprefix BMSB_ -nthreads 2
 ```
