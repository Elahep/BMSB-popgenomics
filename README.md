# Population genomic analysis of the brown marmorated stink bug, _Halyomorpha halys_

Here we will analyze population genomic structure and identify putative loci associated with invasion success in the brown marmorated stink bug (BMSB).
We start by using the ddRAD data of the BMSB generated in <a href="https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07678-z" title="Yan et al 2020.">Yan et al 2020.</a> Also, we will use the whole genome assembly of the BMSB from <a href="https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6510-7" title="Sparks et al 2020.">Sparks et al 2020.<a>
  
 
  ## Downloading genomic data from NCBI

  
Demultiplexed ddRAD reads are available in the NCBI under the BioProject ID: <a href="https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA675311." title="PRJNA675311.">PRJNA675311.</a> (SRA files).
Whole genome assembly data are available in the NCBI as Genbank assembly accession <a href="https://www.ncbi.nlm.nih.gov/assembly/GCF_000696795.1/" title="GCA_000696795.1.">GCA_000696795.1.</a>
  
Whole genome assembly can be downloaded as .tar file. Use the `tar -xf` command to extract .tar file and get .fna.gz which is basically a nucleic acid fasta file. To get .fa file simply do: `mv REFERENCE_GENOME.fna REFERENCE_GENOME.fa`

 Count the number of scaffolds:
  
`grep -c "^>" REFERENCE_GENOME.fa`
  
  
Using <a href="https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump" title="SRAtoolkit">SRAtoolkit</a> we can download SRA files and extract FASTQ files from them. We know the SRR accession number for the BMSB samples are SRR13005202-SRR13005590. We can use a simple bash loop to get all the 389 files. As these are paired-end reads, we will get two FASTQ files per sample.
 
```
module load sratoolkit
for i in {202..590};
do
fasterq-dump --skip-technical SRR13005${i}
done
```

  
  ## Aligning ddRAD reads against the reference genome
  
We will use <a href="http://bio-bwa.sourceforge.net/bwa.shtml" title="BWA">BWA</a> to map reads to the BMSB reference genome. 

 First, we need to index the reference genome (for more information in indexing see <a href="https://www.nature.com/articles/nbt0509-455" title="Trapnell & Salzberg 2009">Trapnell & Salzberg 2009</a>).
 
```
module load BWA
bwa index REFERENCE_GENOME.fa
```
  
As our ddRAD sequences are 150 bp we will use BWA-MEM to align the reads to the reference genome. For shorted sequences (up to 100 bp) use BWA-backtrack.
I will run the program with default features (mostly recommended).
  
```
module load BWA
for i in {202..590};
do
bwa mem -t 8 REFERENCE_GENOME.fa SRR13005${i}_1.fastq.gz SRR13005${i}_2.fastq.gz | gzip > SRR13005${i}.sam.gz
done
  ```
  
 Now one .sam file (Sequence Alignment Map) will be generated for each sample. Use less XX.sam.gz to see the header. The lines starting with @ followed by a two-letter tag shows different parts of header. Examples:
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
  samtools view SRR13005202.sam.gz | head -4
  ```
  
  ```
  SRR13005202.1	65	NW_020111013.1	207614	31	138M	NW_020110848.1	222837	0	AAGACAAATACACCGAGATGTTTGAGACGATCGAAGTACAGTGCCGAATGCAATTCCTATCCCTGGTACTGAAGATTCAGCTGGCAAAAGAACTAGGGAAGAAAGAGGGAATGCCGAGAGACCGCTCGTTTGTGGAAT	FFFFFFF:FFFFFFFFFFFF,,FF:FFF:FFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFF:FFFF:FFFFFFFFFFFFFFFFF,FFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF,FFFFFFF	NM:i:3	MD:Z:43T18A42T32	MC:Z:135M10S	AS:i:123	XS:i:108	XA:Z:NW_020112830.1,-22765,138M,6;
SRR13005202.1	129	NW_020110848.1	222837	5	135M10S	NW_020111013.1	207614	0	CGATACATGTAGGGACACTGTCGGGCTATATGTCCCTCTGCCCCACATTCAAAATATCTCGTGCGCGGGTAATGCGCCACTAGATGCCCAACCCTCCCGCAACCATGACACTTGCGGTCCTCCCGCCGACCATTGTTTCTGGGTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF:F:FFFFFFFF:FFFFFFFFFFFFFFFFFF,FFFFFF::FFFFFFFFFF,,F,::F::FFFFFFFF	NM:i:3	MD:Z:54C7A6C65	MC:Z:138M	AS:i:120	XS:i:115	XA:Z:NW_020113259.1,+11437,135M10S,4;NW_020114712.1,+18591,145M,7;NW_020113353.1,+54411,135M10S,6;NW_020113948.1,-35551,2S78M4I61M,10;
SRR13005202.2	97	NW_020110846.1	5749	0	116M22S	NW_020110616.1	493707	0	TATTTTTTCTTTAATTTATTGTAATTGCAATTATGATGGAGAGCTAAATAAGGTTTAAAATATTTTCATCTGCTCAACCTGGTCTCGATTCATCAAACTTCAGGAGGACACGCGCCACTGGACCACCAAGTACCCCCT	FFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:7	MD:Z:11A27G11A19T0T9C24G8	MC:Z:61M2D84M	AS:i:81	XS:i:81
SRR13005202.2	145	NW_020110616.1	493707	56	61M2D84M	NW_020110846.1	5749	0	AGCAAAAGATTGCCCACGTGGTAGACATTCCTTTCTTCCCATCTCTTTCGCCTACTGTACAGAGGGTTGTAGCCAGGAAACCTATTCTTCTATTGACCTCTTTCTCATCCACTATTTCTCTCATATTCTTCATGATTTGGTTATA	FFFFFFFFFFF:FFFFFFF:FFFFFFFF,FF:FFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFF,FFFFFFFFFFF	NM:i:5	MD:Z:1A59^TT30C37A15	MC:Z:116M22S	AS:i:125	XS:i:96
```

 Each line corresponds to alignment information for a single read. Check out <a href=" https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/04_alignment_quality.html" title="this tutorial">this tutorial</a> to find out the meaning of each field in the alignment section.
  
  
  



