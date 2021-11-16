# Population genomic analysis of the brown marmorated stink bug, _Halyomorpha halys_

Here we will analyze population genomic structure and identify putative loci associated with invasion success in the brown marmorated stink bug (BMSB).
We start by using the ddRAD data of the BMSB generated in <a href="https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07678-z" title="Yan et al 2020.">Yan et al 2020.</a> Also, we will use the whole genome assembly of the BMSB from <a href="https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6510-7" title="Sparks et al 2020.">Sparks et al 2020.<a>
  
 
  ## Downloading genomic data from NCBI
 
Demultiplexed ddRAD reads are available in the NCBI under the BioProject ID: <a href="https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA675311." title="PRJNA675311.">PRJNA675311.</a> (SRA files)
Whole genome assembly data are available in the NCBI as Genbank assembly accession <a href="https://www.ncbi.nlm.nih.gov/assembly/GCF_000696795.1/" title="GCA_000696795.1.">GCA_000696795.1.</a>
  
Using <a href="https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump" title="SRAtoolkit">SRAtoolkit</a> we can download SRA files and extract FASTQ files from them. We know the SRR accession number for the BMSB samples are SRR13005202-SRR13005590. We can use a simple bash loop to get all the 389 files. As these are paired-end reads, we will get two FASTQ files per sample.
  
  `module load sratoolkit`
  `for i in {202..590};`
  `do`
  `fasterq-dump --skip-technical SRR13005${i}`
  `done`

  
  

  




