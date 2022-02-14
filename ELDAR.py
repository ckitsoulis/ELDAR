from Bio import SeqIO
from Bio.SeqUtils import GC
import plotly.express as px
import pandas as pd
import numpy as np
import argparse

def main():

    usage = "A tool to compute basic and additional metrics for genome assembly, based on fragments' length."
    tool_name = "ELDAR"
    footer = "Who:\n Christos Kitsoulis (https://github.com/ckitsoulis) \n \nWhere: \n Genomics and Bioinformatics Group at IMBBC, HCMR"

    parser = argparse.ArgumentParser(description = usage, prog = tool_name, epilog = footer, formatter_class = argparse.RawDescriptionHelpFormatter,)

    parser.add_argument("-i", "--input", metavar = "path", dest = "filename", required = True, help = "absolute path to FASTA file that contains contigs or scaffolds.")
    parser.add_argument("-n", metavar = "int", type = int, nargs = "+", help = "values for additional metrics you want to calculate (e.g. 90 for computing of N90 and L90).")
    parser.add_argument("-p", metavar = "yes OR no", default = "no", type = str, help = "optional argument for computing the percentage of each base. yes or no (default)")
    parser.add_argument("-s", "--sorted", dest = "sort", nargs = "?", const = "", default = None, help = "mode to return a sorted version of FASTA file, in descending order based on length.")
    parser.add_argument("-o", "--output", metavar = "output_filename", type = str, dest = "output", required = True, help = "output file name, without the CSV extension.")
    
    args = parser.parse_args()

    input_file = args.filename
    percentages_arg = args.p
    sorted_file = args.sort
    output_file = args.output

    file = str(input_file)

    ## Extraction of fragments' IDs and their sequences.
    def parse_fasta(file_path):
        HEADERS = []
        SEQUENCES = []
        for record in SeqIO.parse(file_path, 'fasta'):
            HEADERS.append(record.description)
            SEQUENCES.append(record.seq)
        
        return HEADERS, SEQUENCES

    contigs, sequences = parse_fasta(file)

    def _sorted_file_(headers, sequences):
        sorted_ = sorted(zip(headers, sequences), key = lambda x: len(x[1]), reverse = True)
        
        if (sorted_file is not None):
            records_ = [SeqIO.SeqRecord(j,i,"","") for i,j in sorted_]
            SeqIO.write(records_, "sorted.fasta", "fasta")
        
        sorted_headers, sorted_seqs = [], []
        
        for H, S in sorted_:
            sorted_headers.append(H)
            sorted_seqs.append(str(S))
        return sorted_headers, sorted_seqs
    
    SORTED_HEADERS, SORTED_SEQUENCES = _sorted_file_(contigs, sequences)

    ## Calculation of #contigs, assembly size, longest contig size, average contig length.
    contigs_number = len(SORTED_SEQUENCES)
    stats = [contigs_number]

    assembly_size = sum([len(l) for l in SORTED_SEQUENCES])
    stats.append(assembly_size)

    longest_contig = len(SORTED_SEQUENCES[0])
    stats.append(longest_contig)

    average_contig_length = round(assembly_size/contigs_number, ndigits=3)
    stats.append(average_contig_length)

    cumulative_contigs = np.cumsum([len(l) for l in SORTED_SEQUENCES])

    temp = ''.join(SORTED_SEQUENCES)

    ## Calculation of N50-75 and L50-75
    position50 = np.min(np.where(cumulative_contigs>= (assembly_size*0.5)))
    N50 = len(SORTED_SEQUENCES[position50])
    L50 = position50+1

    position75 = np.min(np.where(cumulative_contigs>= (assembly_size*0.75)))
    N75 = len(SORTED_SEQUENCES[position75])
    L75 = position75+1

    stats.extend([N50, N75, L50, L75])

    ## Calculation of GC content
    GC_content = GC(temp)
    stats.append(round(GC_content, ndigits = 3))

    names = ["Number of contigs", "Total length", "Longest contig", "Average contig length", "N50", "N75", "L50", "L75", "GC(%)"]

    ## Calculation of additional N and L statistics.
    if (args.n is not None):
        numbers_list = list(args.n)
        
        Ns, Ls, N_stats, L_stats = [], [], [], []

        for number in numbers_list:
            
            if number in list(range(5,105,5)):
                if (number != 50) or (number != 75):
                    ratio = assembly_size*(number/100)
                    contig_position = (np.min(np.where(cumulative_contigs >= ratio)))
                    L = contig_position+1
                    N = len(SORTED_SEQUENCES[contig_position])
                
                    Ns.append("N{}".format(number))
                    Ls.append("L{}".format(number))
                    N_stats.append(N)
                    L_stats.append(L)
                
                else:
                    pass

            else:
                pass

        names = names + Ns + Ls
        stats = stats + N_stats + L_stats

    ## Calculation of bases' percentage
    if percentages_arg == "yes":
        bases = ["A", "T", "G", "C", "N"]
        percentages = []

        for base in bases:
            percentages.append(round(temp.count(base)/assembly_size, ndigits=3))
        
        names.extend(bases)
        stats.extend(percentages)

    ## Produce interactive plot for cumulative length 
    fig = px.line(x=[i+1 for i in range(len(cumulative_contigs))], y=cumulative_contigs, 
    labels = {'x': "N-th contig (in descending order)", "y": "cumulative length (bp)"}, 
    title = "GEMQT: Genome assembly cumulative length plot")
    fig.write_html("report.html")

    statistics = pd.DataFrame(list(zip(names,stats)), columns = ["Assembly", " "])

    statistics.to_csv("{}.csv".format(output_file), sep=",", index = False)

if __name__ == "__main__":
    main()
