#!/usr/bin/python3
#-*- coding : utf-8 -*-


# __authors__ = ("GONCALVES RIBEIRO Alessandra, IBRAHIM AMOUKOU Najat")
# __contact__ = ("alessandra.goncalves-ribeiro@etu.umontpellier.fr, najat.ibrahim-amoukou@etu.umontpellier.fr")
# __version__ = "0.0.1"
# __date__ = "10/15/2024"
# __licence__ ="This program is free software: you can redistribute it and/or modify
#         it under the terms of the GNU General Public License as published by
#         the Free Software Foundation, either version 3 of the License, or
#         (at your option) any later version.
#         This program is distributed in the hope that it will be useful,
#         but WITHOUT ANY WARRANTY; without even the implied warranty of
#         MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#         GNU General Public License for more details.
#         You should have received a copy of the GNU General Public License
#         along with this program. If not, see <https://www.gnu.org/licenses/>."


     
    ### OPTION LIST:
    
        ##-h or --help : help information
        ##-i or --input: input file (.sam)
        ##-o or --output: output name files (.txt)

    #Synopsis:
        ##SamReader.py -h or --help # launch the help.
        ##SamReader.py -i or --input <file> # Launch SamReader to analyze a samtools file (.sam) and print the result in the terminal
        ##SamReader.py -i or --input <file> -o or --output <name> # Launch SamReader to analyze a samtools file (.sam) and print the result in the file called <name>
  


############### IMPORT MODULES ###############

import os
import re
import sys

# File path
file_path = '/home/najat/mapping/src/mapping.sam'

# --- Functions ---

def fileVerifier(file_path):
    """
    Check if the file is a valid SAM file.
    
    :param file_path: path of the file to check
    :return: True if the file is valid, else False
    """
    # Check if it's a regular file
    if not os.path.isfile(file_path):
        print(f"Error: {file_path} is not valid.")
        return False

    # Check if the file has a .sam extension
    if not file_path.endswith('.sam'):
        print(f"Error: {file_path} doesn't have a '.sam' extension.")
        return False

    # Check if the file is empty
    if os.path.getsize(file_path) == 0:
        print(f"Error: file '{file_path}' is empty.")
        return False

    print(f"File '{file_path}' is valid and not empty.")

    # Validate the number of columns in the first 3 non-header lines
    with open(file_path, 'r') as file:
        line_count = 0
        for line in file:
            # Ignore headers
            if line.startswith('@'):
                continue

            # Count the number of columns separated by tabulation
            num_columns = len(line.strip().split('\t'))

            if num_columns < 11:
                print(f"Error: line '{line.strip()}' only has {num_columns} columns.")
                return False

            # Stop after 3 lines
            line_count += 1
            if line_count == 3:
                break

    print(f"File '{file_path}' has the expected number of columns.")
    return True


def countReads(file_path):
    """
    Analyze the SAM file to count different types of reads.
    
    :param file_path: path to the SAM file
    """
    total_reads = 0
    unmapped_reads = 0
    duplicated_reads = 0
    mapped_reads = 0

    with open(file_path, 'r') as file:
        for line in file:
            # Ignore headers
            if line.startswith('@'):
                continue

            total_reads += 1
            fields = line.strip().split('\t')
            flag = int(fields[1])

            # Check flags
            if flag & 4:  # unmapped
                unmapped_reads += 1
            elif flag & 1024:  # duplicated reads
                duplicated_reads += 1
            else:
                mapped_reads += 1  # mapped reads

    print(f"\n--- Read Statistics ---")
    print(f"Total number of reads: {total_reads}")
    print(f"Total number of unmapped reads: {unmapped_reads}")
    print(f"Total number of duplicated reads: {duplicated_reads}")
    print(f"Total number of mapped reads: {mapped_reads}")


def readPerChrom(file_path):
    """
    Count the number of reads mapped to each chromosome.
    
    :param file_path: path to the SAM file
    :return: dictionary with counts of reads per chromosome
    """
    chromosome_counts = {}

    with open(file_path, 'r') as file:
        for line in file:
            # Ignore headers
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            chromosome = fields[2]

            flag = int(fields[1])
            if flag & 4 == 0:  # if bit 4 is not set, the read is mapped
                chromosome_counts[chromosome] = chromosome_counts.get(chromosome, 0) + 1

    print("\n--- Reads per Chromosome ---")
    for chrom, count in chromosome_counts.items():
        print(f"{chrom}: {count}")


def readPerMAPQ(file_path):
    """
    Count the number of reads for each MAPQ score.
    
    :param file_path: path to the SAM file
    :return: dictionary with counts of reads per MAPQ score
    """
    MappingQ_counts = {}

    with open(file_path, 'r') as file:
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





# --- Main Execution ---

def main():
    # Verify the file
    if fileVerifier(file_path):
        # Count reads
        countReads(file_path)

        # Count reads per chromosome
        readPerChrom(file_path)

        # Count reads per MAPQ
        readPerMAPQ(file_path)


if __name__ == "__main__":
    main()



## 3/ Store,







# ## 4/ Analyse 

# #### Convert the flag into binary ####
# def flagBinary(flag) :

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


 
# #### Summarise the results ####

# def Summary(fileName):
    
   

# # #### Main function ####

#  def main(argv):
    

# # ############### LAUNCH THE SCRIPT ###############
#     if __name__ == "__main__":
#     main(sys.argv[1:])