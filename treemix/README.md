Treemix uses allele frequency data of SNPs to infer maximum likelihood trees representing evolutionary relationships between populations. Treemix can allow a number of migration events in the tree.


The input file for treemix consists of a header with a space-delimited list of the names of populations, followed by lines containing the allele counts at each SNP (very similar to BayPass input file).
We can use Stacks to convert VCF to treemix input file:


```
module load Stacks
populations -V ./H2_BMSB.vcf - O ./ -M ./popmap.txt --treemix
```

Note that the resulting file from Stacks has a header that we need to remove before using it in treemix.

This input file needs to be in _gz_ format so zip the file using:

```
gzip H2_final.treemix
```

We will consider 1-7 migration events in the tree and generating bootstrap replicates by resampling blocks of 500 SNPs. The ML tree is rooted with Japan's population:

```
module load Miniconda3
for i in {1..6};
do
treemix -i H2_final.treemix.gz -root aJapan -m ${i} -bootstrap -k 500 -o BMSB_${i} > treemix_${i}_log 
done
```

