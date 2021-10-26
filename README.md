# GenomeAssembly-Statistics

A script to calculate statistics for genome/draft assemblies coming from a FASTA file, additionally to common N50 - L50, alongside with auto-calculated total assembly length. The statistics are saved in a .csv file.


## Arguments

| Argument | Description |
| --- | --- |
| -i, --input | full path to file, in FASTA format (.fasta) containing contigs or scaffolds. |
| -n | values of statistics you want to calculate. Example: 90 for calculation of N90 and L90. |
| -o, --output | output filename that contains the desired statistics without the CSV extension (.csv). |


## Usage

```bash
python GA_StatisticsCalculator.py -i /users/ckitsoulis/draft_assemblies/genome_assembly.fasta -n 50 85 90 85 -o statistics
```
&nbsp;

### Python libraries needed: [Numpy](https://numpy.org), [Pandas](https://pandas.pydata.org/), [Biopython](https://biopython.org/)

&nbsp;

&nbsp;

&nbsp;

&nbsp;

Who:

*Christos Kitsoulis [Github](https://github.com/ckitsoulis)*

Where:

*Genomics & Bioinformatics group (Genome Nerds) at IMBBC, HCMR*

When:

*During master thesis, October 2021*
