# GenomeAssembly-Statistics

A scripts to calculate statistics for genome/draft assemblies coming from a FASTA file, additionally to common N50 - L50, defined by the user. The statistics are saved in a .csv file.

##Arguments

| Argument | Description |
| --- | --- |
| -i, --input | path to file, in FASTA format (.fasta) containig contigs or scaffolds. |
| -n | values of statistics you want to calculate. Example: 90 for calculation of N90 and L90. |
| -o, --output | output filename that contains the desired statistics without the CSV extension (.csv). |
