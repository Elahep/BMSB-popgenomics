# Genome-wide scan analysis to detect outlier SNPs associated with invasiveness

Here we will use the recently developed C2 statistics <a href="https://academic.oup.com/mbe/article/37/8/2369/5821433" title="(Olazcuaga et al 2020)">(Olazcuaga et al 2020)</a> implemented in <a href="http://www1.montpellier.inra.fr/CBGP/software/baypass/" title="BayPass v2.3">BayPass v2.3</a> to identify outlier SNPs potentially associated with invasion success in the BMSB. 

***

## summary of the approach
For analysing loci associated with a binary trait such as invasive vs. non-invasive, we will use the contrast analysis in BayPass v2.3.
For such an analysis, the new version of BayPass has developed a nonparametric counterpart for the association model implemented in older version of BayPass.
This new nonparametric model relies on a contrast statistic, C2, that compares the standardized population allele frequencies (i.e., the allele frequencies corrected for the
population structure) between the two groups of populations specified by the binary covariable of interest. However, previous versions of BayPass relied on a paramteric models to calculate Bayes Factor (BF). Using the core model of Baypass, we can actually estimate both BF and C2.
 
 
## 1-preparing input files

Using the <a href="https://gitlab.com/YDorant/Toolbox/-/blob/master/reshaper_baypass.py" title="reshaper_baypass.py">reshaper_baypass.py</a> script by <a href="https://gitlab.com/YDorant/Toolbox" title="Yann Dorant">Yann Dorant</a> we can convert VCF to the appropriate BayPass genotype file:

`python ./H1_bialSNP_MAF_geno_LD.vcf popmap.txt BMSB.geno`


We will also need a "contrast file" to assign each population to a binary trait, namely "invasive" and "native" in our analysis. The contrast file has one line and n columns, where n denotes the number of sampled populations. 
