# Genome-wide scan analysis to detect outlier SNPs associated with invasiveness

Here we will use the recently developed C2 statistics <a href="https://academic.oup.com/mbe/article/37/8/2369/5821433" title="(Olazcuaga et al 2020)">(Olazcuaga et al 2020)</a> implemented in <a href="http://www1.montpellier.inra.fr/CBGP/software/baypass/" title="BayPass v2.3">BayPass v2.3</a> to identify outlier SNPs potentially associated with invasion success in the BMSB. 

***

## Summary of the approach
For analysing loci associated with a binary trait such as invasive vs. non-invasive, we will use the contrast analysis in BayPass v2.3.
For such an analysis, the new version of BayPass has developed a nonparametric counterpart for the association model implemented in older versions of BayPass.
This new nonparametric model relies on a contrast statistic, C2, that compares the standardized population allele frequencies (i.e., allele frequencies corrected for the
population structure) between the two groups of populations specified by the binary covariable of interest. However, previous versions of BayPass relied on paramteric models to calculate Bayes Factor (BF). Using the core model of Baypass, we can simultaneously estimate both BF and C2. In this tutorial we will be estimating C2.
 
 
## 1-Preparing input files

Using the <a href="https://gitlab.com/YDorant/Toolbox/-/blob/master/reshaper_baypass.py" title="reshaper_baypass.py">reshaper_baypass.py</a> script by <a href="https://gitlab.com/YDorant/Toolbox" title="Yann Dorant">Yann Dorant</a> we can convert VCF to the appropriate BayPass genotype file:

```
python ./reshaper_baypass.py H1_bialSNP_MAF_geno_LD_reordered.vcf popmap.txt BMSB.geno
```

The number of rows in the .geno file is equal to the number of SNPs and the number of columns is twice the number of populations (because we have two alleles for each site). The first number in each pairs of columns represents the number of the copies of the first allele in that population and the second number represents the number of the copies of the alternative allele. 0 indicates no copies or missing alleles.

We will also need a "contrast file" to assign each population to a binary trait, namely "invasive" and "native" in our analysis. The contrast file has one line and n columns (space seperated), where n denotes the number of sampled populations. In this file, we will allocate native populations to the first group (1) and invasive populations to the second/alternative group (-1). We can also exclude any population from pairwise comparisons by allocating them to the 0 group.
For the BMSB dataset, we have previously reordered individuals in the finel VCF using VCFtools so that we first have our native populations (Japan and China) which are followed by all invasive ones (10 invasive populations):

1 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1


## 2- Running BayPass

On nesi, we need to load ifort.

```
 module load ifort
 ./i_baypass -gfile BMSB.geno -contrastfile contrast.ecotype -outprefix BMSB_ -nthreads 2
 ```



As BayPass runs a MCMC analysis we need to perform several independent MCMC runs (e.g., 3-5), using different seeds. We then need to compare the estimates of parameters such as Ω or statistics such as BF or C2 across various runs to empirically check the convergence of chains.


## 3- Sanity check

Before choosing the candidate SNPs, we must check the behaviour of p-values to make sure our results are not affected by problems associated with multiple testing. We will check the distribution of p-values using a histogram. A well-behaved p-value histogram must be close to uniform for higher p-values (see the plot below). Check this <a href="http://varianceexplained.org/statistics/interpreting-pvalue-histogram/" title="tutorial">tutorial</a> to see what kinds of distributions you can expect in your p-value histograms. 

```
BMSB.C2=read.table("BMSB_summary_contrast.out",h=T)
hist(10**(-1*BMSB.C2$log10.1.pval.),freq=F,breaks=40)
abline(h=1)
```

![image](https://user-images.githubusercontent.com/13001264/178163848-e74149c4-0829-46e3-a40e-9130264a6873.png)


Another step in checking the reliability of our results is to assess convergence and the reproducibility of the MCMC estimates. In the previuos part, we recommended running BayPass for three more independent runs using different seeds. Here, we will calculate Forstner and Moonen distance (FMD) between the Ω matrix of these independent runs using the R function fmd.dist (distributed within the BayPass package). Small FMD distance denotes convergence of the MCMC chains. We can also plot correlation between the Ω parameter and the XtX statistic across different independent BayPass runs to confirm the convergence of MCMC chains in each run.


```
##calculating FMD
source("./baypass_utils.R")

mat0 <- as.matrix(read.table("BMSB_mat_omega.out"))   ##the original run
mat1 <- as.matrix(read.table("MCMC convergence check/run1/BMSB_run1_mat_omega.out"))   ##first replicate
mat2<- as.matrix(read.table("MCMC convergence check/run2/BMSB_run2_mat_omega.out"))  ##second replicate
mat3 <- as.matrix(read.table("MCMC convergence check/run3/BMSB_run3_mat_omega.out"))  ##third replicate


FMD01 <- fmd.dist(mat0,mat1)
FMD02 <- fmd.dist(mat0,mat2)
FMD03 <- fmd.dist(mat0,mat3)


##correlation plots

library("ggpubr")
library(plyr)

run0_omega <- read.table("../BMSB_summary_lda_omega.out", header = TRUE)
run1_omega <- read.table("./run1/BMSB_run1_summary_lda_omega.out", header = TRUE)
run2_omega <- read.table("./run2/BMSB_run2_summary_lda_omega.out", header = TRUE)
run3_omega <- read.table("./run3/BMSB_run3_summary_lda_omega.out", header = TRUE)

o_runs <- as.data.frame(cbind(run0_omega[,5], run1_omega[,5], run2_omega[,5], run3_omega[,5]))

o_runs <- rename(o_runs, c("V1"="main_run", "V2"="replicate_1", "V3"="replicate_2",
                       "V4"="replicate_3"))

o_mainvsrun1 <- ggscatter(o_runs, x = "main_run", y = "replicate_1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "main run omega", ylab = "replicate_1 omega") +
          geom_point(size = 2) +
         theme(axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10)) +
         theme(axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) 

o_mainvsrun2 <- ggscatter(o_runs, x = "main_run", y = "replicate_2", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "pearson",
                        xlab = "main run omega", ylab = "replicate_2 omega") +
                        geom_point(size = 2) +
                        theme(axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10)) +
                        theme(axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) 

o_mainvsrun3 <- ggscatter(o_runs, x = "main_run", y = "replicate_3", 
                        add = "reg.line", conf.int = TRUE, 
                        cor.coef = TRUE, cor.method = "pearson",
                        xlab = "main run omega", ylab = "replicate_3 omega") +
                        geom_point(size = 2) +
   theme(axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10)) +
   theme(axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) 

run0_xtx <- read.table("../BMSB_summary_pi_xtx.out", header = TRUE)
run1_xtx <- read.table("./run1/BMSB_run1_summary_pi_xtx.out", header = TRUE)
run2_xtx <- read.table("./run2/BMSB_run2_summary_pi_xtx.out", header = TRUE)
run3_xtx <- read.table("./run3/BMSB_run3_summary_pi_xtx.out", header = TRUE)

x_runs <- as.data.frame(cbind(run0_xtx[,4], run1_xtx[,4], run2_xtx[,4], run3_xtx[,4]))

x_runs <- rename(x_runs, c("V1"="main_run", "V2"="replicate_1", "V3"="replicate_2",
                           "V4"="replicate_3"))

x_mainvsrun1 <- ggscatter(x_runs, x = "main_run", y = "replicate_1", 
                          add = "reg.line", conf.int = TRUE, 
                          cor.coef = TRUE, cor.method = "pearson",
                          xlab = "main run xtx", ylab = "replicate_1 xtx") +
  geom_point(size = 2) +
  theme(axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) 

x_mainvsrun2 <- ggscatter(x_runs, x = "main_run", y = "replicate_2", 
                          add = "reg.line", conf.int = TRUE, 
                          cor.coef = TRUE, cor.method = "pearson",
                          xlab = "main run xtx", ylab = "replicate_2 xtx") +
  geom_point(size = 2) +
  theme(axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) 

x_mainvsrun3 <- ggscatter(x_runs, x = "main_run", y = "replicate_3", 
                          add = "reg.line", conf.int = TRUE, 
                          cor.coef = TRUE, cor.method = "pearson",
                          xlab = "main run xtx", ylab = "replicate_3 xtx") +
  geom_point(size = 2) +
  theme(axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) 

ggarrange(o_mainvsrun1,o_mainvsrun2,o_mainvsrun3,
          x_mainvsrun1,x_mainvsrun2,x_mainvsrun3,ncol = 3, nrow = 2)
          
          
## 4- Choosing candidate SNPs

`BMSB__summary_contrast.out` contains posterior mean of the C2 contrast statistics (M_C2), standard deviation of C2 contrast statistics (SD_C2), calibrated estimator of C2 statistics (C2_std) and its corrected p value (log10(1/pval)). We can use the 0.001 p-value threshold (0.1%, recommended in the tutorial) or 0.01 threshold (1%, reported in Olazcuaga et al 2020) as a cut-off for considering significantly differentiated outlier SNPs. Any SNP with a p-value (in -log10 scale) higher than this threshold is considered as an outlier. Additionally, we need to calculate false discovery rate and control for multiple testing. To do that, we will calculate q-value for each inividual SNP and will consider SNPs with q-values < 0.1 (10% false discovery rate) as the best candidates for being subject to selective sweeps.

Let's use R to extract the outlier SNPs using the 0.1% p value threshold and calculate q-value for choosing the candidate SNPs:

```
##import the BayPass output file and scaffold list
BMSB.C2=read.table("BMSB_summary_contrast.out",h=T)
scaffolds = read.table("scaffold_list.txt") ##I previuosly extracted scaffold names from VCF using bash commands: cat H1_bialSNP_MAF_geno_LD_reordered.vcf | grep -v "#" | cut -f3 > scaffold_list.txt
BMSB.C2 = as.data.frame(cbind(BMSB.C2, scaffolds))

##make a simple Manhattan plot to check the distribution of outlier SNPs
plot(BMSB.C2$log10.1.pval.)
abline(h=3,lty=2) _#0.001 p--value theshold_

##calculate q-value for each SNP and choose SNPs with q-values < 0.1
library(qvalue)
pvalues <- as.data.frame(10^-BMSB.C2[,6]) ##convert -log p values to normal p values
BMSB.C2 <- cbind(pvalues,BMSB.C2)  ##append it to the df
pvalues <- as.vector(BMSB.C2$`10^-BMSB.C2[, 6]`)
qval <- qvalue(p = pvalues)
plot(qval$qvalues)
abline(h=0.1)  ##qvalues < 0.1 are highly sig.

qvalues = as.data.frame(qval$qvalues)

BMSB.C2 <- cbind(ChvsAll.C2,qvalues)

##extract those SNPs with q value < 0.1
selected_SNPs = BMSB.C2[BMSB.C2$qvalues < 0.1, ]
write.table(selected_SNPs,"BMSB_C2SNPs_qval01.txt", sep = "\t")
```


We can do different pairwise comparisons using BayPass, for example comparing only Japan aginast all invasive populations or only China versus all invasive populations, and use a Venn diagram to check the common SNPs among different comparisons. We will use R to create the Venn diagram and export the list of common SNPs:

```
##draw the Venn diagram
library(ggvenn)
WWvsAll = read.delim("./WWvsInvas_C2SNPsBiggerthan3.txt") ##both Japan and China versus invasive populations
ChvsAll = read.delim("./CHvsInvas_C2SNPsBiggerthan3.txt") ##China versus invasive populations
JPvsAll = read.delim("./JPvsALL_C2SNPsBiggerthan3.txt") ##Japan ersus invasive populations
x= list(ChJPvsAll=WWvsAll$names, ChvsAll=ChvsAll$names, JPvsAll=JPvsAll$names)
ggvenn(x, fill_color = c("#00b4d8", "#c1121f", "#fdf0d5"), fill_alpha = 0.5, stroke_size = 0.2, set_name_size = 4, stroke_color = "navy")

#get the list of common SNPs
library(gplots)
common_snps = venn(x, show.plot = FALSE)
common_snps_list = attributes(common_snps)$intersections$`WwvsAll:ChvsAll:JPvsAll`
write.table(common_snps_list, "BayPass_3comparisons_commonSNPs.txt", sep = "\t")
```

Annotating the outlier SNPs can tell us potential genes/proteins associated with the invasive status of our studied species.





