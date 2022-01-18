# Population genomic analysis of the brown marmorated stink bug, _Halyomorpha halys_

Here we will analyze population genomic structure and identify putative loci associated with invasion success in the brown marmorated stink bug (BMSB).
We start by using the ddRAD data of the BMSB generated in <a href="https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07678-z" title="Yan et al 2020.">Yan et al 2020.</a> Also, we will use the whole genome assembly of the BMSB from <a href="https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6510-7" title="Sparks et al 2020.">Sparks et al 2020.<a>
	
  ## Downloading genomic data from NCBI
 
Demultiplexed ddRAD reads are available in the NCBI under the BioProject ID: <a href="https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA675311." title="PRJNA675311.">PRJNA675311.</a> (SRA files).
Whole genome assembly data are available in the NCBI as Genbank assembly accession <a href="https://www.ncbi.nlm.nih.gov/assembly/GCF_000696795.1/" title="GCA_000696795.1.">GCA_000696795.1.</a>
  
Whole genome assembly can be downloaded as .tar file. Use the `tar -xf` command to extract .tar file and get .fna.gz which is basically a nucleic acid fasta file. To get .fa file simply do: `mv REFERENCE_GENOME.fna REFERENCE_GENOME.fa`

 Count the number of scaffolds:
  
`grep -c "^>" REFERENCE_GENOME.fa`
  
Using <a href="https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump" title="SRAtoolkit">SRAtoolkit</a> we can download SRA files and extract FASTQ files from them. We know the SRR accession numbers for the BMSB samples are SRR13005202-SRR13005590. We can use a simple bash loop to get all the 389 files. As these are paired-end reads, we will get two FASTQ files per sample.
 
```
module load sratoolkit
for i in {202..590};
do
fasterq-dump --skip-technical SRR13005${i}
done
```

  
  ## Aligning ddRAD reads against reference genome
  
We will use <a href="http://bio-bwa.sourceforge.net/bwa.shtml" title="BWA">BWA</a> to map reads to the BMSB reference genome. 

 First, we need to index the reference genome (for more information about the importance of indexing see <a href="https://www.nature.com/articles/nbt0509-455" title="Trapnell & Salzberg 2009">Trapnell & Salzberg 2009</a>).
 
```
module load BWA
bwa index REFERENCE_GENOME.fa
```
  
As our ddRAD sequences are 150 bp we will use BWA-MEM to align the reads to the reference genome. For shorter sequences (up to 100 bp) use BWA-backtrack.
I will run the program with default features (mostly recommended) and adding the -M flag for marking secondary reads (for Picard compatibility) and the -R flag for including the read group header line (@RG).
  
```
module load BWA
for i in {202..590};
do
bwa mem -t 8 -M -R "@RG\tID:SRR_${i}\tLB:SRR_${i}\tPL:ILLUMINA\tPM:HISEQ\tSM:SRR_${i}" REFERENCE_GENOME.fa \ 
	SRR13005${i}_1.fastq.gz SRR13005${i}_2.fastq.gz > align_${i}.sam
done
  ```
  
 Now one .sam file (Sequence Alignment Map) will be generated for each sample. Use `less` to see the header. The lines starting with @ followed by a two-letter tag shows different parts of header. Examples:
 ``` 
  @HD  The header line
	VN: format version
	SO: Sorting order of alignments

	@SQ  Reference sequence dictionary
	SN: reference sequence name
	LN: reference sequence length
	SP: species

	@PG  Program
	PN: program name
	VN: program version
```
 
To check the alignment section of .sam file we can use SAMtools:
  
 ```
  module load SAMtools
  samtools view  align_202.sam.gz | head -4
  ```
  
  ```
SRR13005202.1   113     NW_014468270.1  21080   0       138M    NW_014467082.1  262233  0       ATTCCACAAACGAGCGGTCTCTCGGCATTCCCTCTTTCTTCCCTAGTTCTTTTGCCAGCTGAATCTTCAGTACCAGGGATAGGAATTGCATTCGGCACTGTACTTCGATCGTCTCAAACATCTCGGTGTATTTGTCTT   F
                                                                                                                                                NM:i:6  MD:Z:57A17T18A1C27A3A9  MC:Z:10S135M    AS:i:108        XS:i:108        RG:Z:SRR_202 XA:Z:NW_014466443.1,+9800,138M,6;NW_014468263.1,+9566,138M,8;NW_014467409.1,-138553,123M15S,6;
SRR13005202.1   177     NW_014467082.1  262233  0       10S135M NW_014468270.1  21080   0       AACCCAGAAACAATGGTCGGCGGGAGGACCGCAAGTGTCATGGTTGCGGGAGGGTTGGGCATCTAGTGGCGCATTACCCGCGCACGAGATATTTTGAATGTGGGGCAGAGGGACATATAGCCCGACAGTGTCCCTACATGTA
        FFFFFFFF::F::,F,,FFFFFFFFFF::FFFFFF,FFFFFFFFFFFFFFFFFF:FFFFFFFF:F:FFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF       NM:i:3  MD:Z:65G6T7G54  MC:Z:138M       AS:i:120        XS:i:117        RG:Z:S
        XA:Z:NW_014466550.1,-2027608,7S138M,5;NW_014466682.1,+109695,135M10S,4;NW_014468472.1,-25362,145M,7;NW_014467966.1,+97520,135M10S,6;
SRR13005202.2   97      NW_014466417.1  2512789 0       118M20S NW_014466990.1  295431  0       TATTTTTTCTTTAATTTATTGTAATTGCAATTATGATGGAGAGCTAAATAAGGTTTAAAATATTTTCATCTGCTCAACCTGGTCTCGATTCATCAAACTTCAGGAGGACACGCGCCACTGGACCACCAAGTACCCCCT   F
                                                                                                                                                NM:i:7  MD:Z:9C27A0T0G10T31C29T5        MC:Z:61M2D84M   AS:i:83 XS:i:81 RG:Z:SRR_202
SRR13005202.2   145     NW_014466990.1  295431  56      61M2D84M        NW_014466417.1  2512789 0       AGCAAAAGATTGCCCACGTGGTAGACATTCCTTTCTTCCCATCTCTTTCGCCTACTGTACAGAGGGTTGTAGCCAGGAAACCTATTCTTCTATTGACCTCTTTCTCATCCACTATTTCTCTCATATTCTTCATG
                FFFFFFFFFFF:FFFFFFF:FFFFFFFF,FF:FFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFF,FFFFFFFFFFF       NM:i:5  MD:Z:1A59^TT30C37A15    MC:Z:118M20S    AS:i:125        XS:i:9
```

 Each line corresponds to alignment information for a single read. Check out <a href="https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/04_alignment_quality.html" title="this tutorial">this tutorial</a> to find out the meaning of each field in the alignment section.
  
  

In order to do downstream analysis (e.g. SNP calling) we need sorted BAM files (Binary Alignment Map) and mark duplicates to be further discarded by the SNP caller. We can  use SAMtools or GATK to sort SAM files and mark duplicates. 
If using SAMtools these steps should be followed:
	
```
##first convert sam to bam and sort by name:
for i in {206..590};                                       
do                                                         
samtools sort -n -@ 8 align_${i}.sam.gz -o nSorted_${i}.bam ##@ is the number of threads
done
	
##now fixmate with -m option:
for i in {206..590};
do
samtools fixmate -m nSorted_${i}.bam fixmate_${i}.bam
done

## sort by coordinate:
for i in {206..590};
do
samtools sort -@ 8 fixmate_${i}.bam -o coordSort_${i}.bam
done

## mark duplicates:
for i in {206..590};
do
samtools markdup -@ 8 coordSort_${i}.bam markdup_${i}.bam
done
``` 
	
If using GATK the procedure is simpler:
	
```
module load GATK

for i in {208..306};
do
gatk MarkDuplicatesSpark -I align_${i}.sam -O sorted_dedup_${i}.bam; 
done
```
	
If some downstream applications needed indexed .bam files you can create index using:
	
`samtools index SRR13005202.bam`
  
We can explore the properties of sequence data (reads) in .bam files using samtools flagstats. 
	
`samtools flagstat SRR13005202.bam`

Or the following command for more comprehensive stats:

`samtools stats  SRR13005202.bam`

Runnig the help command for samtools flagstats will give information on what each FLAGS argument represents:
	
```
Each FLAGS argument is either an INT (in decimal/hexadecimal/octal) representing
a combination of the following numeric flag values, or a comma-separated string
NAME,...,NAME representing a combination of the following flag names:

0x1     1  PAIRED         paired-end / multiple-segment sequencing technology
0x2     2  PROPER_PAIR    each segment properly aligned according to aligner
0x4     4  UNMAP          segment unmapped
0x8     8  MUNMAP         next segment in the template unmapped
0x10    16  REVERSE        SEQ is reverse complemented
0x20    32  MREVERSE       SEQ of next segment in template is rev.complemented
0x40    64  READ1          the first segment in the template
0x80   128  READ2          the last segment in the template
0x100   256  SECONDARY      secondary alignment
0x200   512  QCFAIL         not passing quality controls or other filters
0x400  1024  DUP            PCR or optical duplicate
0x800  2048  SUPPLEMENTARY  supplementary alignment	
```

By knowing what each flag means (see above), we can edit .bam files, for example we can remove PCR duplicates or filter out unmapped reads:

`samtools view -F 4 -b SRR13005202.bam > SRR13005202_mapped.bam # this new bam only contains mapped reads, if using -f (instead of -F) we will only retain unmapped reads`

We can also use <a href="http://qualimap.conesalab.org/" title="Qualimap">Qualimap</a> to check the number of reads and mapping quality in each .bam file.
	
```
module load Miniconda3
qualimap multi-bamqc -d Qualimap_bamList.txt -r  ##the .txt file contains two columns, first the name of sample (without the .bam extension) and second the path to the sample (with the name and .bam extension)
```
	

## Variant calling

We will use <a href="https://samtools.github.io/bcftools/howtos/variant-calling.html" title="bcftools mpileup">bcftools mpileup</a> to call variants. 
We will use default parameters and only include sites with minimum mapping quality (--min-MQ) and minimum base quality (--min-BQ) higher than 20.
bcftools mpileup needs the list of all .bam files:
	
```
head bam_list_all.txt
sorted_dedup_202.bam
sorted_dedup_203.bam
sorted_dedup_204.bam
.
.
.
```
	
	
```
module load BCFtools
bcftools mpileup -Ou -f ../BMSB/ref_Hhal10/Hhal10.fa --bam-list bam_list_all.txt --min-MQ 20 --min-BQ 20 | bcftools call -mv -Ob -o calls_all_samples.bcf
```

Variant quality information can be extracted from the output of mpileup using BCFtools. This <a href="https://www.reneshbedre.com/blog/vcf-fields.html" title="tutorial">tutorial</a> explains what each field means.
	
```
bcftools query calls.bcf -f '%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\t%QUAL\n' > calls_bcf__FS.SOR.MQRS.RPRS.QD.MQ.DP.QUAL.txt

	##note that in my dataset, I did not have information for some of these fileds. My current dataset only has FS, MQ, DP, QUAL
```
	
Having information about variant information fields, we can produce histograms or report summary statistics for each filed and obtain information to set thresholds for further filtering of variants.

	
Now using BCFtools we can remove non-biallelic SNPs and indels:
	
`bcftools view -m2 -M2 -v snps,indels -O b -o H1_bialSNP.bcf`
	
We can also filter SNPs based on a MAF threshold:
	
`bcftools filter -i 'MAF > 0.03' -O v -o H1_bialSNP_MAF.vcf H1_bialSNP.bcf`
	
	
Now we will use PLINK to filter SNPs based on a missing genotype threshold. The following code filters out all SNPs with missing call rate exceeding 0.3.
		
`plink --vcf H1_bialSNP_MAF.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --geno 0.3 --recode vcf --out H1_bialSNP_MAF_geno`

PLINK can also be used to filter highly linked (correlated) SNPs:
	
`plink --vcf H1_bialSNP_MAF_geno.vcf --allow-extra-chr --set-missing-var-ids @:# --make-bed --indep-pairwise 10 10 0.8 --recode vcf --out H1_bialSNP_MAF_geno_LD`
	
	
This final VCF contains high quality SNPs. We can use VCFtools to check missing SNP per sample and exclude those samples with high missing data.
VCFtools can also be used to rename and reorder samples in the final VCF.
	
PLINK can be used to do a PCA and have a preliminary view of any population structure in the dataset:
`plink --vcf H1_bialSNP_MAF_geno_LD.vcf --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out H1_bialSNP_MAF_geno_LD_pca`

