#!/usr/bin/python3
#-*- coding : utf-8 -*-

############### IMPORT MODULES ############### 

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os,re,sys,argparse





     # OPTION LIST:
# -h or --help: Displays help information.
# -i or --input: Specifies the path to the input SAM file (.sam).
# -o or --output: Specifies the path to the output file 
# -cR or --count-reads: Counts total reads, mapped reads, unmapped reads,duplicated reads, with a filtering option on the mapping quality.
# -rC or --reads-per-chrom: Counts the number of reads per chromosome.
# -rMQ or --reads-per-mapq: Counts the number of reads for each MAPQ score.
# -cRF or --count-reads-by-flag: Counts the number of reads for each FLAG value.
# -sR or --save-results: Saves the results and associated plots to an HTML file.
# -fS or  --filterSam: Filter SAM file by mapping quality and create a new filterd sam file.

# SYNOPSIS:
# SamReader.py -h or --help               # Displays help.
# SamReader.py -i <file>                  # Check if sam file is valid
# SamReader.py -i <file> -cR              # Prints read statistics (total, mapped, unmapped, duplicated).
# SamReader.py -i <file> -cR  -m <mappingquality>  #Prints read statistics (total, mapped, unmapped, duplicated,filtered).
# SamReader.py -i <file> -rC              # Prints the number of reads per chromosome.
# SamReader.py -i <file> -rMQ             # Prints the number of reads per MAPQ score.
# SamReader.py -i <file> -cRF             # Prints the number of reads per FLAG value.
# SamReader.py -i <file> -sR              # Saves results (including plots) to an HTML file.
# SamReader.py -i <file> -cR -sR          # Combines read statistics and saves results to HTML.
# SamReader.py -i <file> -fS -o <outpufile> -m <mappingQuality>  # Filter SAM file by mapping quality and create a filterd sam file.



############### FUNCTIONS TO :

## 1/ Check, 



def fileVerifier(filePath):
    """
    Check if the file is a valid SAM file.
    
    :param filePath: path of the file to check
    :return: True if the file is valid, else False
    """
    # Check if it's a regular file
    if not os.path.isfile(filePath):
        print(f"Error: {filePath} is not valid.")
        return False

    # Check if the file has a .sam extension
    if not filePath.endswith('.sam'):
        print(f"Error: {filePath} doesn't have a '.sam' extension.")
        return False

    # Check if the file is empty
    if os.path.getsize(filePath) == 0:
        print(f"Error: file '{filePath}' is empty.")
        return False

    print(f"File '{filePath}' is valid and not empty.")

    # Validate the number of columns in the first 3 non-header lines
    with open(filePath, 'r') as file:
        line_count = 0
        for line in file:
            # Ignore headers
            if line.startswith('@'):
                continue

            # Count the number of columns separated by tabulation
            numColumns = len(line.strip().split('\t'))

            if numColumns < 11:
                print(f"Error: line '{line.strip()}' only has {numColumns} columns.")
                return False

            # Stop after 3 lines
            line_count += 1
            if line_count == 3:
                break

    print(f"File '{filePath}' has the expected number of columns.")
    return True


# 2/ Read

## flag dictionnary
flags = {
    0: "aligned forward",
    1: "multiple segments",
    2: "properly aligned",
    4: "unmapped",
    8: "mate unmapped",
    16: "aligned reverse",
    32: "mate aligned reverse",
    64: "first segment",
    128: "second segment",
    256: "not primary",
    512: "fails quality checks",
    1024: "duplicate",
    2048: "supplementary"
}




## Count reads and filter


def countReads(file_path, minQ=0):
    """
    Analyze the SAM file to count different types of reads and filter based on MAPQ score,
    using the FLAG dictionary to display detailed information.

    :param file_path: Path to the SAM file.
    :param minQ: Minimum MAPQ score for filtering reads. Defaults to 0.

    """
    filtered_reads = 0
    total_reads = 0
    unmapped_reads = 0
    duplicated_reads = 0
    mapped_reads = 0
    flag_details = {}  # Dictionary to count occurrences of each flag description

    with open(filePath, 'r') as file:
        for line in file:
            # Ignore header lines
            if line.startswith('@'):
                continue

            total_reads += 1
            fields = line.strip().split('\t')
            flag = int(fields[1])

            # Count occurrences of each FLAG description
            for key, description in flags.items():
                if flag & key:
                    flag_details[description] = flag_details.get(description, 0) + 1

            # Check flags for specific categories
            if flag & 4:  # Unmapped
                unmapped_reads += 1
            elif flag & 1024:  # Duplicated reads
                duplicated_reads += 1
            else:
                mapped_reads += 1  # Mapped reads
                mapq = int(fields[4])  # MAPQ score is in the 5th column

                # Apply the MAPQ filter
                if mapq >= minQ:
                    filtered_reads += 1

    # Print the results
    print(f"\n--- Read Statistics ---")
    
    print(f"Total number of unmapped reads: {unmapped_reads}")
    print(f"Total number of duplicated reads: {duplicated_reads}")
    print(f"Total number of mapped reads: {mapped_reads}")
    print(f"Total number of filtered reads (MAPQ >= {minQ}): {filtered_reads}")
    print(f"\n")
    for description, count in flag_details.items():
        print(f"{description}: {count}")



def readPerChrom(filePath):
    """
    Count the number of reads mapped to each chromosome.
    
    :param filePath: path to the SAM file
    :return: dictionary with counts of reads per chromosome
    """
    chromosomeCounts = {}

    with open(filePath, 'r') as file:
        for line in file:
            # Ignore headers
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            chromosome = fields[2]

            flag = int(fields[1])
            if flag & 4 == 0:  # if bit 4 is not set, the read is mapped
                chromosomeCounts[chromosome] = chromosomeCounts.get(chromosome, 0) + 1

    print("\n--- Reads per Chromosome ---")
    for chrom, count in chromosomeCounts.items():
        print(f"{chrom}: {count}")


def readPerMAPQ(filePath):
    """
    Count the number of reads for each MAPQ score.
    
    :param filePath: path to the SAM file
    :return: dictionary with counts of reads per MAPQ score
    """
    MappingQ_counts = {}

    with open(filePath, 'r') as file:
        for line in file:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            mapq = fields[4]

            flag = int(fields[1])
            if flag & 4 == 0:  # if bit 4 is not set, the read is mapped
                MappingQ_counts[mapq] = MappingQ_counts.get(mapq, 0) + 1

    print("\n--- Reads per MAPQ ---")
    for mapq, count in MappingQ_counts.items():
        print(f"MAPQ {mapq}: {count}")



def countReadsByFlags(filePath):
    """
    Count the number of reads for each flag value in a SAM file.
    
    Args:
    filePath (str): Path to the .sam file.
    
    Returns:
    dict: A dictionary with flag values as keys and counts as values.
    """
    flagCounts = {}

    with open(filePath, 'r') as file:
        for line in file:
            if line.startswith('@'):  # Skip header lines in SAM file
                continue
            
            columns = line.split('\t')
            flag = int(columns[1])  # The flag is the second column
            
            # Count occurrences of each flag
            if flag in flagCounts:
                flagCounts[flag] += 1
            else:
                flagCounts[flag] = 1
    
    return flagCounts



def plotFlagCounts(flagCounts):
    """
    Plot the count of reads for each flag value.
    
    Args:
    flagCounts (dict): A dictionary with flag values as keys and counts as values.
    """
    # Convert the dictionary to a pandas DataFrame for easy plotting
    df = pd.DataFrame(list(flagCounts.items()), columns=['Flag', 'Count'])
    
    # Plotting the data
    plt.figure(figsize=(10, 6))
    plt.bar(df['Flag'].astype(str), df['Count'], color='skyblue')
    plt.xlabel('Flag Value')
    plt.ylabel('Read Count')
    plt.title('Number of Reads for Each Flag Value')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

# Example usage
# filePath = 'C:/Users/aless/Downloads/mapping/src/mapping.sam'
# flagCounts = countReadsByFlags(filePath)
# plotFlagCounts(flagCounts)


###################################Data vizualizing


from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import pandas as pd

def saveResults(flagCounts, plotFileName='plot.png', htmlFileName='results.html'):
    """
    Save the flag counts and plot to an HTML file.
    
    Args:
    flagCounts (dict): Dictionary of flag counts.
    plotFileName (str): The name of the image file for the plot.
    htmlFileName (str): The name of the HTML file to save.
    """
    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(flagCounts.keys(), flagCounts.values(), color='skyblue')
    ax.set_xlabel('Flag Value')
    ax.set_ylabel('Read Count')
    ax.set_title('Number of Reads for Each Flag Value')
    
    # Save the plot as an image
    canvas = FigureCanvas(fig)
    canvas.print_figure(plotFileName)

    # Generate HTML content
    html_content = f"""
    <html>
        <head>
            <title>Read Analysis Results</title>
        </head>
        <body>
            <h1>Read Analysis</h1>
            <h2>Flag Counts</h2>
            <table border="1">
                <tr>
                    <th>Flag Value</th>
                    <th>Read Count</th>
                </tr>
    """
    for flag, count in flagCounts.items():
        html_content += f"""
                <tr>
                    <td>{flag}</td>
                    <td>{count}</td>
                </tr>
        """
    
    html_content += f"""
            </table>
            <h2>Graph of Read Counts by Flag</h2>
            <img src="{plotFileName}" alt="Flag Counts Graph">
        </body>
    </html>
    """
    
    # Save HTML file
    with open(htmlFileName, 'w') as html_file:
        html_file.write(html_content)
    print(f"Results saved to {htmlFileName}")



def executeFlagStats():
   countReadsByFlags(filePath)
   plotFlagCounts(flagCounts)



##Filter sam file by mapping quality


def filterSam(file_path, output_file, minQ=30):
    """
    Filters the SAM file by MAPQ score and writes the filtered reads to a new file.

    :param file_path: Path to the input SAM file.
    :param output_file: Path to the output SAM file.
    :param minQ: Minimum MAPQ score to retain a read (default is 30).
    """
    with open(file_path, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Skip header lines 
            if line.startswith('@'):
                outfile.write(line)
                continue
            
            # Split the line into columns
            fields = line.strip().split('\t')
            mapq = int(fields[4])  # MAPQ score is in the 5th column (index 4)

            # Filter based on MAPQ score
            if mapq >= minQ:
                outfile.write(line)


##Filter sam file by only keeping mapped reads

def mappedRead(file_path, output_file):
    """
    Filters the SAM file and keeping only the mapped reads  to a new file.

    :param file_path: Path to the input SAM file.
    :param output_file: Path to the output SAM file.
   
    """
    with open(file_path, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Skip header lines 
            if line.startswith('@'):
                outfile.write(line)
                continue
            
            # Split the line into columns
            fields = line.strip().split('\t')
            flag = int(fields[1])  # 

            # Filter based on MAPQ score
            if flag & 4 == 0:  # if bit 4 is not set, the read is mapped
                outfile.write(line)
    




## 3/ Store,

# Main script
import argparse

def main():
    parser = argparse.ArgumentParser(description="Analyze a SAM file and provide various statistics.")
    
    # Arguments for file input and output
    parser.add_argument('-i', '--input', required=True, help="Path to the SAM file.")
    parser.add_argument('-o', '--output', required=False, help="Path to the output SAM file.")
    
    # Arguments for read counting and analysis
    parser.add_argument('-cR', '--count-reads', action='store_true', help="Count reads (total, mapped, etc.).")
    parser.add_argument('-rC', '--reads-per-chrom', action='store_true', help="Count reads per chromosome.")
    parser.add_argument('-rMQ', '--reads-per-mapq', action='store_true', help="Count reads based on MAPQ scores.")
    parser.add_argument('-cRF', '--count-reads-by-flag', action='store_true', help="Count reads based on FLAG values.")
    
    # Arguments for saving results and generating output
    parser.add_argument('-sR', '--saveResults', action='store_true', help="Save results and graphs to an HTML file.")
    
    # Arguments for filtering and saving new SAM files
    parser.add_argument('-fS', '--filterSam', action='store_true', help="Filter SAM file by mapping quality")
    parser.add_argument('-mR', '--mappedRead', action='store_true', help="Filter SAM file by keeping mapped reads")
    
    # Argument for MAPQ score filtering
    parser.add_argument('-m', '--minQ', type=int, default=0, help="Minimum MAPQ score for filtering reads.")
    
    # Argument for executing flag stats (calls executeFlagStats)
    parser.add_argument('-eFS', '--execute-flag-stats', action='store_true', help="Execute Flag Stats analysis.")

    args = parser.parse_args()

    # Display help if no arguments are provided
    if not any(vars(args).values()):
        parser.print_help()
        return

    # Verify the input file (this function should be defined elsewhere in your code)
    if not fileVerifier(args.input):
        return

    # Execute the appropriate functions based on the arguments
    if args.count_reads:
        countReads(args.input, args.minQ)
    if args.reads_per_chrom:
        readPerChrom(args.input)
    if args.reads_per_mapq:
        readPerMAPQ(args.input)
    if args.count_reads_by_flag:
        flag_counts = countReadsByFlags(args.input)
        if args.saveResults:
            saveResults(flag_counts)
    if args.filterSam:
        filterSam(args.input, args.output, args.minQ)
    if args.mappedRead:
        mappedRead(args.input, args.output)
    if args.execute_flag_stats:
        executeFlagStats()

if __name__ == "__main__":
    main()


def executeFlagStats(filePath):
    """
    Exécuter l'analyse des statistiques sur les flags du fichier SAM.
    
    Args:
    filePath (str): Chemin du fichier SAM à analyser.
    
    Returns:
    dict: Un dictionnaire des comptes de flags.
    """
    # Chamar a função de contagem de flags
    flagCounts = countReadsByFlags(filePath)
    # Chamar a função de plot (se necessário)
    plotFlagCounts(flagCounts)
    
    return flagCounts

## 4/ Analyse 



# #### Convert the flag into binary ####
#  def flagBinary(flag) :

#     flagB = bin(int(flag)) # Transform the integer into a binary.
#     flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
#     flagB = list(flagB) 
#     if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
#         add = 12 - len(flagB) # We compute the difference between the maximal flag size (12) and the length of the binary flag.
#         for t in range(add):
#             flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
#     return flagB


# #### Analyze the unmapped reads (not paired) ####
# def unmapped(sam_line):
    
#     unmapped_count = 0
#     with open ("only_unmapped.fasta", "a+") as unmapped_fasta, open("summary_unmapped.txt", "w") as summary_file:
#         for line in sam_line:
#             col_line = line.split("\t")
#             flag = flagBinary(col_line[1])

#             if int(flag[-3]) == 1:
#                 unmapped_count += 1
#                 unmapped_fasta.write(toStringOutput(line))

#         summary_file.write("Total unmapped reads: " + str(unmapped_count) + "\n") 
#         return unmapped_count

# #### Analyze the partially mapped reads ####
# def partiallyMapped(sam_line):
    
#     partially_mapped_count = 0

#     with open ("only_partially_mapped.fasta", "a+") as partillay_mapped_fasta, open("summary_partially_mapped.txt", "w") as summary_file:
#         for line in sam_line:
#             col_line = line.split("\t")
#             flag = flagBinary(col_line[1]) # We compute the same 

#             if int(flag[-2]) == 1: 
#                 if col_line[5] != "100M":
#                     partially_mapped_count += 1
#                     partillay_mapped_fasta.write(toStringOutput(line))

#         summary_file.write("Total partially mapped reads: " + str(partially_mapped_count) + "\n") 
#         return partially_mapped_count


# ### Analyse the CIGAR = regular expression that summarise each read alignment ###
# def readCigar(cigar): 
   
#     ext = re.findall('\w',cigar) # split cigar 
#     key=[] 
#     value=[]    
#     val=""

#     for i in range(0,len(ext)): # For each numeric values or alpha numeric
#         if (ext[i] == 'M' or ext[i] == 'I' or ext[i] == 'D' or ext[i] == 'S' or ext[i] == 'H' or ext[i] == "N" or ext[i] == 'P'or ext[i] == 'X'or ext[i] == '=') :
#             key.append(ext[i])
#             value.append(val)
#             val = ""
#         else :
#             val = "" + val + ext[i]  # Else concatenate in order of arrival
    
#     dico = {}
#     n = 0
#     for k in key:   # Dictionnary contruction in range size lists              
#         if k not in dico.keys():    # for each key, insert int value
#             dico[k] = int(value[n])   # if key not exist, create and add value
#             n += 1
#         else:
#             dico[k] += int(value[n])  # inf key exist add value
#             n += 1
#     return dico

# ### Analyse the CIGAR = regular expression that summarise each read alignment ###
# def percentMutation(dico):
        
#     totalValue = 0 # Total number of mutations
#     for v in dico :
#         totalValue += dico[v]

#     mutList = ['M','I','D','S','H','N','P','X','=']
#     res = ""
#     for mut in mutList : # Calculated percent of mutation if mut present in the dictionnary, else, percent of mut = 0
#         if mut in dico.keys() :
#             res += (str(round((dico[mut] * 100) / totalValue, 2)) + ";")
#         else :
#             res += ("0.00" + ";")
#     return res

# def globalPercentCigar():
#     """
#       Global representation of cigar distribution.
#     """
    
#     with open ("outpuTable_cigar.txt","r") as outpuTable, open("Final_Cigar_table.txt", "w") as FinalCigar:
#         nbReads, M, I, D, S, H, N, P, X, Egal = [0 for n in range(10)]

#         for line in outpuTable :
#             mutValues = line.split(";")
#             nbReads += 2
#             M += float(mutValues[2])+float(mutValues[12])
#             I += float(mutValues[3])+float(mutValues[13])
#             D += float(mutValues[4])+float(mutValues[14])
#             S += float(mutValues[5])+float(mutValues[15])
#             H += float(mutValues[6])+float(mutValues[16])
#             N += float(mutValues[7])+float(mutValues[17])
#             P += float(mutValues[8])+float(mutValues[18])
#             X += float(mutValues[9])+float(mutValues[19])
#             Egal += float(mutValues[10])+float(mutValues[20])

#         FinalCigar.write("Global cigar mutation observed :"+"\n"
#                         +"Alignlent Match : "+str(round(M/nbReads,2))+"\n"
#                         +"Insertion : "+str(round(I/nbReads,2))+"\n"
#                         +"Deletion : "+str(round(D/nbReads,2))+"\n"
#                         +"Skipped region : "+str(round(S/nbReads,2))+"\n"
#                         +"Soft Clipping : "+str(round(H/nbReads,2))+"\n"
#                         +"Hard Clipping : "+str(round(N/nbReads,2))+"\n"
#                         +"Padding : "+str(round(P/nbReads,2))+"\n"
#                         +"Sequence Match : "+str(round(Egal/nbReads,2))+"\n"
#                         +"Sequence Mismatch : "+str(round(X/nbReads,2))+"\n")


 
# ### Summarise the results ####

# def Summary(fileName):
    
   




# #### Main function ####

# def main(argv):
    

# ############## LAUNCH THE SCRIPT ###############

# if __name__ == "__main__":
#     main(sys.argv[1:])


