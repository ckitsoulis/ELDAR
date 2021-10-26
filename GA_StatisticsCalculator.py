from Bio import SeqIO
import pandas as pd
import numpy as np
import argparse

usage = "A program to calculate additional statistics (N and L) for a set of contigs or scaffolds based on their lengths."
tool_name = "draft_stats"
footer = "Who:\n Christos Kitsoulis (https://github.com/ckitsoulis) \n \nWhere: \n Genomics and Bioinformatics Group (Genome Nerds) at IMBBC, HCMR"

parser = argparse.ArgumentParser(description = usage, prog = tool_name, epilog = footer, formatter_class = argparse.RawDescriptionHelpFormatter,)
parser.add_argument('-i',
'--input', 
metavar = 'path',
dest = 'filename',  
required = True, 
help = 'path to file, in FASTA format (.fasta) containig contigs or scaffolds.')

parser.add_argument('-n', 
metavar = 'int', 
type = int, 
nargs = '+', 
required = True, 
help = 'values of statistics you want to calculate. Example: 90 for calculation of N90 and L90.')

parser.add_argument('-o', 
'--output', 
metavar = 'output_filename',
type = str, 
dest = 'output', 
required = True, 
help = 'output filename that contains the desired statistics without the CSV extension (.csv).')

args = parser.parse_args()

input_file = args.filename
numbers_list = list(args.n)
output_file = args.output

file = str(input_file)

##Extraction of fragments' IDs and their sequences.
def parse_fasta(file_path):
    HEADERS = []
    SEQUENCES = []
    for record in SeqIO.parse(file_path, 'fasta'):
        HEADERS.append(record.description)
        SEQUENCES.append(str(record.seq))
    
    return HEADERS, SEQUENCES

contigs, sequences = parse_fasta(file)

##Calculation of each fragment size.
contigs_size = [len(s) for s in sequences]

genome_size = sum(contigs_size)

contigs_dictionary = dict(zip(contigs,contigs_size))

sorted_contigs = dict(sorted(contigs_dictionary.items(), key=lambda x: x[1], reverse=True))

cumulative_contigs = np.cumsum(list(sorted_contigs.values()))

## Calculation of N and L statistics.
Ns, Ls, N_stats, L_stats = [], [], [], []

for number in numbers_list:
    
    if number in list(range(5,105,5)):

        ratio = genome_size*(number/100)
        contig_position = (np.min(np.where(cumulative_contigs >= ratio)))
        L = contig_position+1
        N = list(sorted_contigs.values())[contig_position]
        
        Ns.append("N{}".format(number))
        Ls.append("L{}".format(number))
        N_stats.append(N)
        L_stats.append(L)


    else:
        pass

names = ["Total length"] + Ns + Ls
stats = [genome_size] + N_stats + L_stats

statistics = pd.DataFrame(list(zip(names,stats)), columns = ["Assembly", "draft"])

statistics.to_csv("{}.csv".format(output_file), sep=",", index = False)
