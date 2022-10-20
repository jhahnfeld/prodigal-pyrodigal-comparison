# prodigal-pyrodigal-comparison
A Python script for the comparison of Prodigal and Pyrodigal predictions.

## Contents

- [Background](#background)
- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Input & Output](#input-and-output)
- [License](#license)

## Background

[Prodigal](https://github.com/hyattpd/Prodigal) is a gene prediction tool for prokaryotic genomes.
It is widely used in bioinformatics, however the latest version is from 2016 and is missing unreleased bug fixes.
[Pyrodigal](https://github.com/althonos/pyrodigal) provides Cython bindings and a Python interface to Prodigal, 
including the unpublished bug fixes and other optimizations.

## Description
The Python script `compare.py` compares the predictions of Prodigal and Pyrodigal to find mismatches or missing 
predictions. The differences in the predictions are stored in a TSV file.

## Installation
In order to compare Prodigal and Pyrodigal correctly, a current version of Prodigal is required, which must be compiled 
by the user, as there is no new version available. 

```bash
git clone https://github.com/hyattpd/Prodigal.git
cd Prodigal
make install
```
If you want to install Prodigal in a custom directory use:
```bash
make install INSTALLDIR=/where/i/want/prodigal/
```

`compare.py` requires the additional Python packages:
* [Biopython](https://github.com/biopython/biopython)
* [xopen](https://github.com/pycompression/xopen)
* [pyrodigal](https://github.com/althonos/pyrodigal)

The packages can be easily installed with Pip or Conda.
### Pip
```bash
python3 -m pip install --user biopython xopen pyrodigal
```

### Conda
```bash
conda install -c conda-forge -c bioconda biopython xopen pyrodigal
```

## Usage
```bash
usage: compare.py [-h] [--genome GENOME [GENOME ...]] [--prodigal PRODIGAL] [--output OUTPUT]

Compare CDS predictions of Prodigal to the predictions of Pyrodigal and save the differences in a TSV file.

optional arguments:
  -h, --help            show this help message and exit
  --genome GENOME [GENOME ...], -g GENOME [GENOME ...]
                        Input genomes (/some/path/*.fasta)
  --prodigal PRODIGAL, -p PRODIGAL
                        Path to a newly compiled Prodigal binary.
  --output OUTPUT, -o OUTPUT
                        Output path (default="./comparison")
```

## Examples
* Download 50 test genomes
```bash
wget -i test_genomes.txt
```

* Input: single compressed genome
```bash
./compare.py --genome /path/to/genome.fasta.gz --prodigal /path/to/binary/prodigal
```

* Input: multiple genomes
```bash
./compare.py --genome /path/to/*.fasta --prodigal /path/to/binary/prodigal
```


## Input and Output
### Input
`compare.py` can use a single or several (compressed) genomes in the FASTA format as input.
### Output
* A short summary of each comparison is printed to stdout: `Hits genome=GCF_000006765: prodigal=5681, pyrodigal=5681, equal=True`
* The `comparison` output directory contains a mismatches.tsv TSV file with differing predictions.
* The `comparison/tmp` directory contains the Prodigal train file and GFF output for each used genome.


## License

GNU General Public License v3.0