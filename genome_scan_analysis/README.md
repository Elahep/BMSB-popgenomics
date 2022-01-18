# Genome-wide scan analysis to detect outlier SNPs associated with invasiveness

Here we will use the recently developed C2 statistics <a href="https://academic.oup.com/mbe/article/37/8/2369/5821433" title="(Olazcuaga et al 2020)">(Olazcuaga et al 2020)</a> implemented in <a href="http://www1.montpellier.inra.fr/CBGP/software/baypass/" title="BayPass v2.3">BayPass v2.3</a> to identify outlier SNPs potentially associated with invasion success in the BMSB. 



## 1-preparing input files

Using the <a href="https://gitlab.com/YDorant/Toolbox/-/blob/master/reshaper_baypass.py" title="reshaper_baypass.py">reshaper_baypass.py</a> script (by <a href="https://gitlab.com/YDorant/Toolbox" title="Yann Dorant"<Yann Dorant</a> ) we can convert VCF to the appropriate BayPass genotype file:

`python ./H1_bialSNP_MAF_geno_LD.vcf popmap.txt BMSB.geno`


We will also need a "contrast file" to assign each population into our binary trait, namely "invasive" and "native" in our analysis. The contrast file has one line and n columns, where n denotes the number of sampled populations. 
