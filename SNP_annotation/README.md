Here we will use snpEff to annotate SNPs and predict their effect on genes and proteins. In simple words, snpEff checks whether a SNP is located in a gene.

We can directly use VCF files (containing either full dataset or just the BayPass outliers) as an input for snpEff. However, as we study non-model species, the annotation information of our studied species is not availabe in snpEff and we need to manually add our database to the program.

**snpEff can be installed through bioconda:

`conda install -c bioconda snpeff -p ~`



## Creating SnpEff database for non-model species

We first need to create our own config file that has genomic information of our species of interest. We need to find the main config file of snpEff, copy it to our working directory and then add our species genome. Let's start in the current working directory which I call _snpEff_analysis_. In this directory we need to create a couple of more directories and copy the config file:

```
$mkdir data
$cd data
$mkdir mygenome
$cd ../
$find ~ -name snpEff.config
/home/eparvizi/.conda/pkgs/snpsift-4.3.1t-hdfd78af_3/share/snpsift-4.3.1t-3/snpEff.config
$cp /home/eparvizi/.conda/pkgs/snpsift-4.3.1t-hdfd78af_3/share/snpsift-4.3.1t-3/snpEff.config ./
```

Now check the config file using `nano`. Make sure the path to data directory in the config file is: `data.dir = ./data/` 

Under the header lines of the _Database & Genomes_ section in the config file add the genome of your species of interest:

```
# my genome
mygenome.genome : BMSB
```

Now we need to copy the genome assembly (fasta file) and associated genome annotation (GFF file) into the _mygenome_ directory that we created in the first step. You must rename the files to "mygenome.fa" and "genes.gff".

Now we are ready to create our own snpEff database. Go back to the _snpEff_analysis_ directory and run the following:

```
snpEff build -c snpEff.config -gff3 -v mygenome > snpEff.stdout 2> snpEff.stderr
```

## Runnig SnpEff

Now we are all sorted to run snpEff and do SNP annotation.

Let's do the annotation on the BayPass outlier SNPs that were obtained from comparing China agianst all invasive populations.

We first need to extract the outlier SNPs from our main VCF file which includes all SNPs. By having the list of outliers in a txt file we can use VCFtools to extract the outlier SNPs:

```
vcftools --vcf H1_bialSNP_MAF_geno_LD_387_renamed_reordered.vcf --snps outlierSNPs_list.txt --recode --recode-INFO-all
mv recode.out.vcf outlierSNPs.vcf
```
Now we can annotate the outlier SNPs:

```
snpEff -c snpEff.config mygenome ./outlierSNPs.vcf > outlierSNPs.annot.vcf
```

The resulting annotated vcf has ANN filed which gives you information on properties of each annotated snp.

Now use snpSift to get a nice txt file that has all the information about the interesting variants from snpEff vcf output:

```
cat outlierSNPs.annot.vcf | /home/eparvizi/share/snpeff-5.0-1/scripts/vcfEffOnePerLine.pl | java -Xmx8g -jar /home/eparvizi/share/snpsift-4.3.1t-3/SnpSift.jar extractFields - CHROM POS ID REF ALT AF "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].BIOTYPE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" > outlierSNPs.annot.sift.txt
```

In the resulting txt file you can see if any of the variants have caused major protein change or is just silent using the IMPACT column.
You can see HIGH, MODERATE, LOW and MODIFIER in this column. The definition of each (and other genomics/transcriptomics terms) can be found at:https://asia.ensembl.org/Help/Glossary.


