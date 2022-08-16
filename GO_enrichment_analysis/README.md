In this tutorial, we will learn how to extract exons (or genes) that are located at flanking regions of the outlier SNPs. We will then use these exons in a gene ontology (GO) analysis. GO term analysis helps us identify biological processes associated with the genes identified as candidates for invasion success.



## 1- Preparing input files from outlier SNPs


We first need to extract the outlier SNPs from our full SNP VCF file.  We need a .txt file containing the list of these outliers with the scaffold and SNP position information for each SNP. The .txt file should look like this:


```
head -n 4 BayPass_outliers.txt
NW_020110192.1:859087
NW_020110206.1:72475
NW_020110209.1:1769865
NW_020110236.1:179165
```

We will use VCFtools to extract these outliers:


```
module load VCFtools
vcftools --vcf BMSB_FullData.vcf --snps BayPass_outliers.txt --recode --recode-INFO-all
mv out.recode.vcf BayPass_outliers.vcf
```

We now need to convert .vcf to .bed file. "bed" format is a text file format that represents genomic data (scaffold or chromosome names and SNP coordinates) in the form of columns separated by spaces or tabs.


```
module load BEDOPS
vcf2bed < BayPass_outliers.vcf > BayPass_outliers.bed
```

The .bed file can now be used to extarct flanking sequence information using the reference genome. BEDTools can be used for this purpose, however, it needs a specific "genome file" to do the job. I used this bash script to create the **BEDTools genome file** from the reference genome fasta index file:


```
awk -v OFS='\t' {'print $1,$2'} Hhal20.fa.fai > Hhal10.txt   ###Hhal20.txt is the special bedtools genome format (to get fai.fa file use: samtools faidx Hhal20.fa)
```

Now all necessary input files are ready to use BEDTools to extract flanking regions (e.g. 10,000 base pairs) around the outlier SNPs:


```
module load BEDTools
bedtools flank -i BayPass_outliers.bed -g Hhal20.txt -b 10000 | cut -f1-3 > BayPass_outliers_10kb_flank.bed
```


## 2- Extracting exons using UCSC Genome Browser

We can use Genome Browser (UCSC) to extract exons or genes around outlier SNPs using the bed file created above. We will use the <a href="https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1427809609_lBBYXSD3L9xru48CQoF5X7Eoar06" title="Table Browser tool" >Table Browser tool</a>.  




![image](https://user-images.githubusercontent.com/13001264/184762373-5015cfea-f3d2-41ab-960e-5253f073155f.png)





After choosing the right genome and assembly for your taxon of interest under the *Select dataset* section, we can upload the BayPass_outliers_10kb_flank.bed file by clicking on *define regions*. We then need to choose *sequence* for our *output format* and give the *output filename* a suitable name. By selecting *get output* we will be directed to a new page which allows us to choose sequence retrieval regions. As we are interested in extracting exonic regions, we will choose *5 UTR Exons*, *CDS Exons*, and *3 UTR Exons*. Lets also obtain *One FASTA record per gene*. After choosing the appropriate options, we can click on *get sequence* which will give us a fasta file containg all exonic regions at 10000 bp upstream and downstream of our outlier SNPs.

The output of the Table Browser tool contains information about the assembly ID used for extracting exons, the transcript ID and range (scaffold and coordinates) for each exon. We only need the **transcript ID** for each exon. In the case of the BMSB, transcript IDs start with "M_". So we will use simple bash scripts to get only the lines that contain the transcript IDs and use notepad+ to manually edit the final file and remove any extra character/information from the list of transcript IDs. We will call the final file "outliers_transcriptIDs.txt".


## 3- GO term analysis using InterProScan and the R package topGO


InterProScan helps with functional characterization of novel nucleotide or protein sequences. Using the DNA sequence of the flanking regions of the SNPs in .fasta format, it can search various protein family databases such as Pfam and PANTHER to obtain GO terms associated with each DNA sequence. The resulting GO terms will then be used in the R package topGO to formally assess which of the GO terms are significantly more represented in our candidate gene set.


To begin with, in addition to the transcipt IDs of the exonic regions around our outlier SNPs (above), we need the full list of all transcript IDs in our BMSB dataset, including candidates and non-candidates, to be able to extract GO annotations for all BMSB genes. This full dataset is needed as an input for the GO enrichment analysis using topGO to calculate whether the GO terms associated with our candidate genes are significantly more represented compared to the other genes in our full dataset. To obtain the full list of all transcripts I have used snpEff (check out the tutorial <a href="https://github.com/Elahep/BMSB-popgenomics/tree/main/SNP_annotation" title="here" >here</a>). I extract transcript IDs of all BMSB SNPs (from the BMSB_FullData.vcf) using the "snpeff_genes.txt" output of snpEff. 

```
snpEff -c snpEff.config mygenome ./BMSB_FullData.vcf > BMSB_FullData.anno.vcf
cat snpEff_genes.txt | grep -v "#" | cut -f3 | uniq > BMSB_FullData_transIDfromsnpeff.txt
```

In order to get the DNA sequences associated with all transcript IDs obtained from snpEff, we can use the <a href="https://www.ncbi.nlm.nih.gov/sites/batchentrez" title="Batch Entrez tool" >Batch Entrez tool</a> of NCBI. We need to choose the *nucleotide* database and upload our full transcript ID list and export .fasta file format containing the DNA sequences associated with the transcript IDs. 




![image](https://user-images.githubusercontent.com/13001264/184761835-18e60d0d-5609-4cc3-a860-f2a998cdb56f.png)





Note that NCBI Batch Entrez only accepts up to 2500 transcript IDs at one time. So we need to create smaller files each containing up to 2500 transcript IDs:

```
cat BMSB_FullData_transIDfromsnpeff.txt | sed -n '1,2500p' > 1-2500.txt
cat BMSB_FullData_transIDfromsnpeff.txt | sed -n '2501,5000p' > 2501-5000.txt
cat BMSB_FullData_transIDfromsnpeff.txt | | sed -n '5001,7500p' > 5001-7500.txt
```

After downloading DNA sequences for each of the above transcript ID lists, we can append the separate .fasta files together:

```
cat *.fasta > FullData_batchentrez.fasta
```

This .fasta file can now be used to obtain GO terms for each transcript ID using InterProScan. This analysis would normally take up to 24 hours to finish using 8 cpus.

```
module load InterProScan
interproscan.sh -i FullData_batchentrez.fasta -t n --goterms
```

Using the above command, InterProScan will output multiple files (such as .gff3, .json and .tsv). We only need the .tsv file to extract GO terms for different transcript IDs. When we look at the .tsv output file, we can see that InterProScan has not found GO terms for some transcript IDs. Also, there can be multiple GO terms suggested for one transcript ID. Lets first remove the lines that do not contain any GO terms. In bash:

```
grep -w "GO" FullData_batchentrez.fasta.tsv | cut -f1,14 > FullData_GOlist.txt 
```

The resulting file contains two columns: in the first column we have transcript IDs of all of our genes (candidate and non-candidate) and in the second column we have GO terms associated with each transcript ID. Note that multiple GO terms in each row are separated by "|".


Now we have all the input files ready to run Go term enrichment analysis using topGo in R. Here I have used the weight01 algorithm and Fisher's exact test, and retained enriched GO terms with a p-value of < 0.05.


```
library(topGO)


##import tab delimited file of all genes and their GO terms
geneID2GO <- readMappings("./InterProScan/FullData_GOlist.txt")  
geneID2GO$XM_014420108.2  #check the GO terms for some of the transcript IDs
str(head(geneID2GO))

geneNames <- names(geneID2GO)
head(geneNames)

##import transcript IDs of outlier SNPs:
interesting_genes = read.table("./outliers_transcriptIDs.txt")
myInterestingGenes <- as.vector(interesting_genes$V1) ##import the vector of your genes of interest. This should be just the ID of the genes.
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)

##create a topGo object for the "Biological Processes" terms (BP):
GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
##get the list of significant genes
sig_genes = sigGenes(GOdata_BP) 
## return the GO graph
graph(GOdata_BP)
##GO enrichment test
resultFisher <- runTest(GOdata_BP, algorithm="weight01", statistic="fisher")  ##resultFisher is a TopGoresult objetc. You can see p values in this object using the "score" function: 
resultFisher  #this shows how many GO terms are significant
allRes <- GenTable(GOdata_BP, raw.p.value = resultFisher, classicFisher = resultFisher, ranksOf = "classicFisher", Fis = resultFisher, topNodes = length(resultFisher@score)) 
allRes
write.table(allRes, "GOresults_BP.txt", sep = "\t")  ##save the results for the BP terms.


##now we can do GO enrichment term analysis using different ontology. For example, we will analyze for the "Molecular Function" terms here (MF):
GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_MF

##GO enrichment test
resultFisher_MF <- runTest(GOdata_MF, algorithm="weight01", statistic="fisher")
resultFisher_MF  #this shows how many GO terms are significant
allRes_MF <- GenTable(GOdata_MF, raw.p.value = resultFisher_MF, classicFisher = resultFisher_MF, ranksOf = "classicFisher", Fis = resultFisher_MF, topNodes = length(resultFisher_MF@score)) 
allRes_MF
write.table(allRes_MF, "GOresults_MF.txt", sep = "\t")

##Similarly, we will do the analysis for the "Cellular Component" terms as well (CC):
GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_CC
resultFisher_CC <- runTest(GOdata_CC, algorithm="weight01", statistic="fisher")
resultFisher_CC  #this shows how many GO terms are significant
allRes_CC <- GenTable(GOdata_CC, raw.p.value = resultFisher_CC, classicFisher = resultFisher_CC, ranksOf = "classicFisher", Fis = resultFisher_CC, topNodes = length(resultFisher_CC@score)) 
allRes_CC
write.table(allRes_CC, "GOresults_CC.txt", sep = "\t")

```

Results of GO enrichment analysis of the BMSB outlier SNPs and candidate genes showed enriched pathways in a range of molecular functions, including oxidoreductase activity (GO:0016491), odorant binding (GO:0005549) and protein binding (GO:0005515), as well as in biological processes such as cellular protein modification (GO:0006464). I did not obtain any significant GO terms for the Cellular Component analysis.

