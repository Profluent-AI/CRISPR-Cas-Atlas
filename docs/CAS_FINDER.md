# CRISPR-Cas Prediction Tool
The CRISPR-Cas discovery pipeline used to build the CRISPR-Cas Atlas

### Install

Install the package
```bash
git clone https://github.com/Profluent-AI/CRISPR-Cas-Atlas
cd CRISPR-Cas-Atlas
mamba env create -f environment.yml
source activate casfinder
bash install.sh
```

### Run example

```bash
cas_finder -s sample001 -t 32 cas_finder/test/input.fna output
```

### Usage

```
Usage: cas_finder [options] <contigs> <outdir> 

positional arguments:
 contigs   input seqs in fasta format
 outdir    output directory

options:
      -n   sample or assembly name
      -m   minimum contig length
      -f   number of flanking bases around each crispr array
      -t   number of threads to use
```

### Additional modules

- For prediction of tracrRNAs, see: https://github.com/skDooley/TRACR_RNA
- For prediction of PAMs, see: https://github.com/Matteo-Ciciani/PAMpredict
