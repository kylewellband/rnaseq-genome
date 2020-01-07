# rnaseq-genome

RNAseq pipeline for species with an annotated genome. Accomodates both single and paired-end data using [fastp](https://github.com/OpenGene/fastp) for quality trimming, [STAR](https://github.com/alexdobin/STAR) for read-mapping and [stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml) for transcript assembly and quantification.

## Pre-requisites

A computer with sufficient memory resources for STAR (~31G for human- / salmon-sized genomes).

A `conda` package and environment management system such as Miniconda v3: https://conda.io/projects/conda/en/latest/user-guide/install/index.html

The repository comes with a conda environment YAML file that will create an environment and install all dependencies.

## Installation

```
git clone "https://github.com/kylewellband/rnaseq-genome.git"
cd ./rnaseq-genome
conda env create -f rnaseqgenome_env.yml
```

## Usage

- Place uncompressed genome fasta file in 02_reference and rename it 'genome.fasta'

- Place uncompressed annotation in GFF3 format in 02_reference and rename it 'gene.gff'

- Load the conda environment `conda activate rnaseqgenome`

- Run the script `./01_scripts/00_index_STAR_reference.sh` to generate the STAR genome reference index.

>*__NOTE:__ It may be necessary to manually adjust (reduce) the values of `--genomeChrBinNbits` and `--genomeSAindexNbases` depending on the memory available on your computer system.*

- Place de-multiplexed RNAseq data in folder: `./03_raw_data/`

>*__NOTE:__ Paired-end data should be suffixed with "\_R1.fastq.gz" and "\_R2.fastq.gz".*

- Execute the remaining scripts in order using the `*_se.sh` scripts for single end data and the `*_pe.sh` scripts for paired-end data.

- The *.ctab files located in `./06_stringtie_gtf/` are suitable for import into R using the [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) package. The files and sample names are listed in `./06_stringtie_gtf/stringtie_ctab_files.txt` which can be edited using a spreadsheet to contain sample and/or treatment specific information.

## License

CC share-alike

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">rnaseq-genome</span> by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Kyle Wellband</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
