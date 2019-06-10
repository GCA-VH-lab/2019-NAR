# Supplementary methods

### *"A role for the _Saccharomyces cerevisiae_ ABCF protein New1 during translation termination"*  
> preliminary version in [bioRxiv](https://www.biorxiv.org/content/10.1101/638064v1)  

### Introduction  

Following pipeline and scripts were used to process Ribo-Seq data published in https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7763.
The pipeline `Pipeline.py` performs data downloading; preprocessing including quality filtering, removing ncRNA; aligning RPFs to genome; and finalizes with corrected P-site assignment. Output of ribosome densities is stored in HDF5 format by default. Python scripts (sub-folder `scripts/`) can read HDF5 format and convert ribosome coverage to bedgraph files `hdf_2_bedgraph.py`; calculate relative 3' UTR coverage `compute_relative_3utr_coverage.py`; and queuing score `compute_queuing.py`. Genes coverage around stop codon for metagene plots are summed up by the script `hdf_2_metagene_tables_stop.py` The script can be controlled by gene list; or let it split data based on stop codon.  Input files, R script and output of differential expression (DE) 
analysis is in sub-folder __`Differential_Expression/`__.

### Prerequisite
1) Python v.3.6 from [Continuum](https://www.continuum.io/downloads) (or some other python). Continuum's Anaconda comes with a bunch of libraries and have a nice package manager `conda`.  Add the _bioconda_ channel ::  `conda config --add channels bioconda`
2) [wget](http://www.gnu.org/software/wget/) There is a detail [installation guide](https://coolestguidesontheplanet.com/install-and-configure-wget-on-os-x/) for OS-X. I would recommend precompiled binaries from [Rudix](http://rudix.org/packages/wget.html) site.  If you download fastq files manually and placed them under folder `1-Raw/` then the step 1, wget, is not necessary. Don't forget rename Fastq files as they are referred in the `Param.in` 
3) [cutadapt](https://cutadapt.readthedocs.io/en/stable/) :: `conda install cutadapt`   
4) `pigz` - optional if not installed cutadapt falls back to single core mode  
5) `HISAT2` v. 2.0.5 or higher `ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads`. `Hisat2` trims ends of reads with bad quality by default. Starting from version 2.0.5 there is an option to turn it off.  
6) [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 
7) [samtools](https://github.com/samtools/samtools/)  :: `conda install samtools`  
8) [pysam](https://github.com/pysam-developers/pysam) :: `conda install pysam`  

### Additional files
#### Sequences and annotation
are placed under sub-folder `0-References/`  
  * Genome.fa  - genome sequence in FastA format  
  * ncRNA.fa   - non coding RNA in FastA format  
  * Genome.gtf - genome annotation v.88 in GTF (gff2) format from _Ensembl_   

Other data files are derived from those and commands for that are listed in the file  `build_index.sh`.
```bash
# create index files for aligners and readjust GTF file for tabix and pysam
cd 0-References/
sh build_index.sh
```
#### Data tables  
`readlength_offsets.txt` - read length specific offsets are used for P-Site assignment.  
`E-MTAB-7763-riboseq.sdrf.txt` -  Samples table from ArrayExpress, contains only Ribo-Seq samples and needed for downloading fastq files (Step 1)  
`Param.in` - Contains different parameters used by Pipeline.py
### Usage
```bash
python  Pipeline.py
```

### Remarks
Annotation in GTF format is originated from _Ensembl_ and it contains features like `stop_codon` and `start_codon` what are essential for the pipeline. In versions 84 - 88  `stop_codon` annotation was missing approximately for 145 ORFs what was corrected at least in the version v.95. We rerunned some analyses with the annotation v.95 and saw subtle or no difference. We keep here annotation version 88 to be consistent with the paper.

### References
The pipelines backbone is based on a code used by [Radhakrishnan, A., et al. Cell (2016)](https://github.com/greenlabjhmi/2016-Cell-Dhh1) (Green's lab.)
