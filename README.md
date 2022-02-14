## ELDAR - Genome Assembly Evaluation Tool

ELDAR is a simple evalution tool for genome assemblies which computes basic metrics. The current script includes options for the computation of additional metrics alongside to common ones: i) N50, ii) total assembly size, iii) fragments number etc. The results are saved in a CSV file while interactive plots of cumulative length distribution as HTML (report.html). Recently, it was added the option to retrieve a sorted version (descending order) of FASTA file, based on fragments' length.


### Arguments

| Argument | Description |
| --- | --- |
| -i, --input | full path to file, in FASTA format (.fasta) containing contigs or scaffolds. |
| -n | values for additional statistics you want to calculate, e.g. 90 for calculation of N90 and L90. |
| -p | optional argument for calculating the percentage of each base. yes or no = default |
| -s, --sorted | option argument for returning a sorted version of file, based on lengths (sorted.fasta) |
| -o, --output | output's filename that contains the results, without the CSV extension (.csv). |


### Usage

ELDAR [-h] -i path [-n int [int ...]] [-p yes OR no] [-s [SORT]] -o output_filename

```bash
python ELDAR.py -i /users/ckitsoulis/draft_assemblies/genome_assembly.fasta -n 85 90 95 -s -p yes -o results
```

### Dependencies 

1. [Numpy](https://numpy.org)

2. [Pandas](https://pandas.pydata.org/)

3. [Biopython](https://biopython.org/)

4. [Argparse](https://pypi.org/project/argparse/) 

5. [Plotly](https://plotly.com/) 

### Notes

The `requirements.txt` file lists all Python libraries that the tool depends on and they will be installed in the environment using the following command:

```bash
pip install -r requirements.txt
```

### Contact

If you have any questions or suggestions, please feel free to contact me via email at chriskits@hotmail.com

&nbsp;

&nbsp;



Who:

*Christos Kitsoulis*

Where:

*Genomics & Bioinformatics group at IMBBC, HCMR*

When:

*October 2021*
