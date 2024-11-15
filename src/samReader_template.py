#!/usr/bin/python3
#-*- coding : utf-8 -*-


__authors__ = ("GONCALVES RIBEIRO Alessandra, IBRAHIM AMOUKOU Najat")
__contact__ = ("alessandra.goncalves-ribeiro@etu.umontpellier.fr, najat.ibrahim-amoukou@etu.umontpellier.fr")
__version__ = "0.0.1"
__date__ = "10/15/2024"
__licence__ ="This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.
        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU General Public License for more details.
        You should have received a copy of the GNU General Public License
        along with this program. If not, see <https://www.gnu.org/licenses/>."


     
    ### OPTION LIST:
    
        ##-h or --help : help information
        ##-i or --input: input file (.sam)
        ##-o or --output: output name files (.txt)

    #Synopsis:
        ##SamReader.py -h or --help # launch the help.
        ##SamReader.py -i or --input <file> # Launch SamReader to analyze a samtools file (.sam) and print the result in the terminal
        ##SamReader.py -i or --input <file> -o or --output <name> # Launch SamReader to analyze a samtools file (.sam) and print the result in the file called <name>
  


############### IMPORT MODULES ###############

import os, sys, re ....


############### FUNCTIONS TO :

## 1/ Check, 


import os

def fileVerifier(file_path):
    """
    Vérifie si le fichier est un fichier SAM valide.
    
    :param file_path: Chemin du fichier à vérifier
    :return: True si le fichier est valide, False sinon
    """
    # Vérification du type de fichier
    if not os.path.isfile(file_path):
        print(f"Erreur : {file_path} n'est pas un fichier régulier.")
        return False

    if not file_path.endswith('.sam'):
        print(f"Erreur : {file_path} n'a pas l'extension '.sam'.")
        return False

    # Vérification si le fichier est vide
    if os.path.getsize(file_path) == 0:
        print(f"Erreur : Le fichier '{file_path}' est vide.")
        return False

    print(f"Le fichier '{file_path}' n'est pas vide.")
    
    # Vérification du nombre de colonnes dans les 3 premières lignes non-en-tête
    with open(file_path, 'r') as file:
        line_count = 0
        for line in file:
            # Ignorer les lignes d'en-tête
            if line.startswith('@'):
                continue

            # Compte le nombre de colonnes (séparées par des tabulations)
            num_columns = len(line.strip().split('\t'))

            if num_columns < 11:
                print(f"Erreur : La ligne '{line.strip()}' contient seulement {num_columns} colonnes.")
                return False

            # Compte jusqu'à 3 lignes
            line_count += 1
            if line_count == 3:
                break

    print(f"Le fichier '{file_path}' respecte le nombre minimum de colonnes.")
    return True

## 2/ Read, 

## Dictionnaire des flags
flags = {
    0: "read aligned to the reference in a forward strand",
    1: "template having multiple segments in sequencing",
    2: "each segment properly aligned according to the aligner",
    4: "read unmapped",
    8: "mate unmapped",
    16: "read aligned to the reverse strand",
    32: "mate aligned to the reverse strand",
    64: "the first segment in the template",
    128: "the second segment in the template",
    256: "not primary alignment",
    512: "read fails platform/vendor quality checks",
    1024: "PCR or optical duplicate",
    2048: "supplementary alignment"
}

def countReads(file_path):
    """
    Analyse un fichier SAM pour compter différents types de lectures.
    
    :param file_path: Chemin du fichier SAM à analyser
    """
    if not file_verifier(file_path):
        return

    total_reads = 0
    unmapped_reads = 0
    duplicated_reads = 0
    mapped_reads = 0

    with open(file_path, 'r') as file:
        for line in file:
            # Ignorer les lignes d'en-tête
            if line.startswith('@'):
                continue
            
            total_reads += 1
            fields = line.strip().split('\t')
            flag = int(fields[1])

            # Vérifier les flags
            if flag & 4:  # Non mappé
                unmapped_reads += 1
            elif flag & 1024:  # Duplicata
                duplicated_reads += 1
            else:
                mapped_reads += 1

    print(f"Le total de reads est: {total_reads}")
    print(f"Le total de reads non mappés est: {unmapped_reads}")
    print(f"Le total de reads en duplicité est: {duplicated_reads}")
    print(f"Le total de reads mappés est: {mapped_reads}")


## read per chromosom

def readPerChrom(file_path):
    if not fileVerifier(file_path):
    return

    # Initialiser un dictionnaire pour stocker le nombre de reads mappés par chromosome
chromosome_counts = {}

# Ouvrir le fichier SAM en mode lecture
with open(file_path, 'r') as file:
    for line in file:
        # Ignorer les lignes d'en-tête (commencent par '@')
        if line.startswith('@'):
            continue

        # Diviser la ligne en champs en utilisant une tabulation comme séparateur
        fields = line.strip().split('\t')
        
        # Extraire le FLAG (champ 2) et le convertir en entier
        flag = int(fields[1])
        
        # Extraire le chromosome (champ 3, aussi appelé RNAME)
        chromosome = fields[2]

        # Vérifier si le read est mappé :
        # FLAG 4 indique que le read est "unmapped", donc si le bit 4 n'est pas défini (flag & 4 == 0), le read est mappé
        if flag & 4 == 0:  
            # Si le chromosome est déjà dans le dictionnaire, incrémenter le compteur pour ce chromosome
            if chromosome in chromosome_counts:
                chromosome_counts[chromosome] += 1
            # Sinon, initialiser le compteur à 1 pour ce chromosome
            else:
                chromosome_counts[chromosome] = 1

# Retourner le dictionnaire contenant le nombre de reads mappés par chromosome
return chromosome_counts




## 3/ Store,







## 4/ Analyse 

#### Convert the flag into binary ####
def flagBinary(flag) :

    flagB = bin(int(flag)) # Transform the integer into a binary.
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB) 
    if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
        add = 12 - len(flagB) # We compute the difference between the maximal flag size (12) and the length of the binary flag.
        for t in range(add):
            flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
    return flagB


#### Analyze the unmapped reads (not paired) ####
def unmapped(sam_line):
    
    unmapped_count = 0
    with open ("only_unmapped.fasta", "a+") as unmapped_fasta, open("summary_unmapped.txt", "w") as summary_file:
        for line in sam_line:
            col_line = line.split("\t")
            flag = flagBinary(col_line[1])

            if int(flag[-3]) == 1:
                unmapped_count += 1
                unmapped_fasta.write(toStringOutput(line))

        summary_file.write("Total unmapped reads: " + str(unmapped_count) + "\n") 
        return unmapped_count

#### Analyze the partially mapped reads ####
def partiallyMapped(sam_line):
    
    partially_mapped_count = 0

    with open ("only_partially_mapped.fasta", "a+") as partillay_mapped_fasta, open("summary_partially_mapped.txt", "w") as summary_file:
        for line in sam_line:
            col_line = line.split("\t")
            flag = flagBinary(col_line[1]) # We compute the same 

            if int(flag[-2]) == 1: 
                if col_line[5] != "100M":
                    partially_mapped_count += 1
                    partillay_mapped_fasta.write(toStringOutput(line))

        summary_file.write("Total partially mapped reads: " + str(partially_mapped_count) + "\n") 
        return partially_mapped_count


### Analyse the CIGAR = regular expression that summarise each read alignment ###
def readCigar(cigar): 
   
    ext = re.findall('\w',cigar) # split cigar 
    key=[] 
    value=[]    
    val=""

    for i in range(0,len(ext)): # For each numeric values or alpha numeric
        if (ext[i] == 'M' or ext[i] == 'I' or ext[i] == 'D' or ext[i] == 'S' or ext[i] == 'H' or ext[i] == "N" or ext[i] == 'P'or ext[i] == 'X'or ext[i] == '=') :
            key.append(ext[i])
            value.append(val)
            val = ""
        else :
            val = "" + val + ext[i]  # Else concatenate in order of arrival
    
    dico = {}
    n = 0
    for k in key:   # Dictionnary contruction in range size lists              
        if k not in dico.keys():    # for each key, insert int value
            dico[k] = int(value[n])   # if key not exist, create and add value
            n += 1
        else:
            dico[k] += int(value[n])  # inf key exist add value
            n += 1
    return dico

### Analyse the CIGAR = regular expression that summarise each read alignment ###
def percentMutation(dico):
        
    totalValue = 0 # Total number of mutations
    for v in dico :
        totalValue += dico[v]

    mutList = ['M','I','D','S','H','N','P','X','=']
    res = ""
    for mut in mutList : # Calculated percent of mutation if mut present in the dictionnary, else, percent of mut = 0
        if mut in dico.keys() :
            res += (str(round((dico[mut] * 100) / totalValue, 2)) + ";")
        else :
            res += ("0.00" + ";")
    return res

def globalPercentCigar():
    """
      Global representation of cigar distribution.
    """
    
    with open ("outpuTable_cigar.txt","r") as outpuTable, open("Final_Cigar_table.txt", "w") as FinalCigar:
        nbReads, M, I, D, S, H, N, P, X, Egal = [0 for n in range(10)]

        for line in outpuTable :
            mutValues = line.split(";")
            nbReads += 2
            M += float(mutValues[2])+float(mutValues[12])
            I += float(mutValues[3])+float(mutValues[13])
            D += float(mutValues[4])+float(mutValues[14])
            S += float(mutValues[5])+float(mutValues[15])
            H += float(mutValues[6])+float(mutValues[16])
            N += float(mutValues[7])+float(mutValues[17])
            P += float(mutValues[8])+float(mutValues[18])
            X += float(mutValues[9])+float(mutValues[19])
            Egal += float(mutValues[10])+float(mutValues[20])

        FinalCigar.write("Global cigar mutation observed :"+"\n"
                        +"Alignlent Match : "+str(round(M/nbReads,2))+"\n"
                        +"Insertion : "+str(round(I/nbReads,2))+"\n"
                        +"Deletion : "+str(round(D/nbReads,2))+"\n"
                        +"Skipped region : "+str(round(S/nbReads,2))+"\n"
                        +"Soft Clipping : "+str(round(H/nbReads,2))+"\n"
                        +"Hard Clipping : "+str(round(N/nbReads,2))+"\n"
                        +"Padding : "+str(round(P/nbReads,2))+"\n"
                        +"Sequence Match : "+str(round(Egal/nbReads,2))+"\n"
                        +"Sequence Mismatch : "+str(round(X/nbReads,2))+"\n")


 
#### Summarise the results ####

def Summary(fileName):
    
   

#### Main function ####

def main(argv):
    

############### LAUNCH THE SCRIPT ###############

if __name__ == "__main__":
    main(sys.argv[1:])