# GenomeAssembly-Statistics

A script to calculate basic statistics for genome/draft assemblies with optional ability for calculating additional ones. The results are saved in a .csv file.


## Arguments

| Argument | Description |
| --- | --- |
| -i, --input | full path to file, in FASTA format (.fasta) containing contigs or scaffolds. |
| -n | values for additional statistics you want to calculate. Example: 90 for calculation of N90 and L90. |
| -p | optional argument for calculating the percentage of each base. yes or no (default) |
| -o, --output | output filename that contains the desired statistics without the CSV extension (.csv). |


## Usage

GA_StatisticsCalculator [-h] -i path [-n int [int ...]] [-p yes OR no] -o output_filename

```bash
python GA_StatisticsCalculator.py -i /users/ckitsoulis/draft_assemblies/genome_assembly.fasta -n 50 85 90 85 [-p yes] -o results
```
&nbsp;

#### Python libraries used: [Numpy](https://numpy.org), [Pandas](https://pandas.pydata.org/), [Biopython](https://biopython.org/)

*The same script, without the need of the above libraries, is going to be uploaded soon.*

&nbsp;

&nbsp;

&nbsp;


Who:

*[Christos Kitsoulis](https://github.com/ckitsoulis)*

Where:

*Genomics & Bioinformatics group (Genome Nerds) at IMBBC, HCMR*

When:

*During master thesis, October 2021*
