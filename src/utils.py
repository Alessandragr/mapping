from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os,re,sys,argparse

from flags import flags



############################# FUNCTIONS TO :

############################# CHECK 

def fileVerifier(filePath):
    """
    Check if the file is a valid SAM file.
    
    :param filePath: path of the file to check
    :return: True if the file is valid, else False
    """
    print(f"Starting verification for the file: {filePath}")

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
        lineCount = 0
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
            lineCount += 1
            if lineCount == 3:
                break

    print(f"File '{filePath}' has the expected number of columns.")
    print(f"The file: '{filePath}' can be used for next step.")
    return True






############################# READ STAT

def countReads(filePath, minQ=0):
    """
    Analyze the SAM file to count different types of reads and filter based on MAPQ score,
    using the FLAG dictionary to display detailed information.

    :param filePath: Path to the SAM file.
    :param minQ: Minimum MAPQ score for filtering reads. Defaults to 0.

    """
    filteredReads = 0
    totalReads = 0
    unmapedReads = 0
    duplicatedReads = 0
    mappedReads = 0
    flagDetails = {}  # Dictionary to count occurrences of each flag description

    with open(filePath, 'r') as file:
        for line in file:
            # Ignore header lines
            if line.startswith('@'):
                continue

            totalReads += 1
            fields = line.strip().split('\t')
            flag = int(fields[1])

            # Count occurrences of each FLAG description
            for key, description in flags.items():
                if flag & key:
                    flagDetails[description] = flagDetails.get(description, 0) + 1
                    

            # Check flags for specific categories
            if flag & 4 or flag & 8:  # Unmapped
                unmapedReads += 1
            elif flag & 1024:  # Duplicated reads
                duplicatedReads += 1
            else:
                mappedReads += 1  # Mapped reads
                mapq = int(fields[4])  # MAPQ score is in the 5th column

                # Apply the MAPQ filter
                if mapq >= minQ:
                    filteredReads += 1

    # Print the results
    print(f"\n--- Read Statistics ---")
    
    print(f"Total number of unmapped reads: {unmapedReads}")
    print(f"Total number of duplicated reads: {duplicatedReads}")
    print(f"Total number of mapped reads: {mappedReads}")
    print(f"Total number of filtered reads (MAPQ >= {minQ}): {filteredReads}")
    print(f"\n")
    for description, count in flagDetails.items():
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
            # Only count reads where both the read and its mate are mapped
            if flag & 4 == 0 and flag & 8 == 0:  # Exclude unmapped reads and mates
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
    mappingQCount = {}

    with open(filePath, 'r') as file:
        for line in file:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            mapq = int(fields[4])

            flag = int(fields[1])
            if flag & 4 == 0 and flag & 8 == 0:  # if bit 4 and 8 are not set, the read is mapped
                mappingQCount[mapq] = mappingQCount.get(mapq, 0) + 1

    sortedMapingQCount = sorted(mappingQCount.items())
    
    print("\n--- Number of Reads ---")
    for mapq, count in sortedMapingQCount:
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


############################# MAPING
def parseSam(filePath):
    """
    Extract sequences from a SAM file.

    :param filePath: Path to the SAM file.
    :return: List of sequences (each as a tuple with read ID and sequence).
    """
    print("Passou no parseSam") 
    sequences = []
    with open(filePath, 'r') as file:
        for line in file:
            if line.startswith('@'):  # Ignore header lines
                continue
            parts = line.strip().split('\t')
            if len(parts) > 9:  # Verifica se há pelo menos 10 colunas
                readId = parts[0]   # Coluna 1: ID da leitura
                sequence = parts[9]  # Coluna 10: Sequência
                sequences.append((readId, sequence))
    print(f"Loaded {len(sequences)} sequences from the file.")
    return sequences




def smithWaterman(seqOne: str, seqTwo: str, matchScore: int = 2, mismatchPenalty: int = -2, gapPenalty: int = -3):
    m, n = len(seqOne), len(seqTwo)
    scoreMatrix = np.zeros((m + 1, n + 1), dtype=int)
    maxScore = 0
    maxPos = None

    # Preenchendo a matriz de pontuação com penalidades ajustadas
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = scoreMatrix[i - 1][j - 1] + (matchScore if seqOne[i - 1] == seqTwo[j - 1] else mismatchPenalty)
            delete = scoreMatrix[i - 1][j] + gapPenalty
            insert = scoreMatrix[i][j - 1] + gapPenalty
            scoreMatrix[i][j] = max(0, match, delete, insert)

            if scoreMatrix[i][j] > maxScore:
                maxScore = scoreMatrix[i][j]
                maxPos = (i, j)

    # Traceback para construir o alinhamento
    alignedSeqOne = []
    alignedSeqTwo = []
    i, j = maxPos

    while scoreMatrix[i][j] != 0:
        if scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + (matchScore if seqOne[i - 1] == seqTwo[j - 1] else mismatchPenalty):
            alignedSeqOne.append(seqOne[i - 1])
            alignedSeqTwo.append(seqTwo[j - 1])
            i -= 1
            j -= 1
        elif scoreMatrix[i][j] == scoreMatrix[i - 1][j] + gapPenalty:
            alignedSeqOne.append(seqOne[i - 1])
            alignedSeqTwo.append('-')
            i -= 1
        else:
            alignedSeqOne.append('-')
            alignedSeqTwo.append(seqTwo[j - 1])
            j -= 1

    return ''.join(reversed(alignedSeqOne)), ''.join(reversed(alignedSeqTwo)), maxScore




def compareSequences(referenceSequences, querySequences, smithWatermanFunc, scoreThreshold=50):
    """
    Compare query sequences with reference sequences using Smith-Waterman algorithm.

    :param referenceSequences: List of sequences from the reference SAM file.
    :param querySequences: List of sequences from the query SAM file.
    :param smithWatermanFunc: Function that implements Smith-Waterman algorithm.
    :param scoreThreshold: Minimum score threshold for considering an alignment valid.
    :return: List of matches with alignment details.
    """
    results = []

    for queryId, querySeq in querySequences:
        for refId, refSeq in referenceSequences:
            # Realiza o alinhamento usando Smith-Waterman
            alignedRef, alignedQuery, score = smithWatermanFunc(refSeq, querySeq)
            
            # Adiciona o resultado se a pontuação for maior que o limiar
            if score >= scoreThreshold:  # Verifica se o alinhamento é relevante
                results.append({
                    "queryId": queryId,
                    "refId": refId,
                    "score": score,
                    "alignedRef": alignedRef,
                    "alignedQuery": alignedQuery
                })
    
    return results



def processSequences(args, smithWaterman):
    """
    Process sequences by comparing reference and query sequences, and print the alignment results.

    :param args: Parsed command-line arguments.
    """
    print("Passou no processSequences") 
    if not args.reference or not args.input or not args.outputFile:
        print("Error: Reference, query, and output file paths are required.")
        return
    
    # Carregar as sequências dos arquivos SAM, ignorando os cabeçalhos
    referenceSequences = parseSam(args.reference)
    querySequences = parseSam(args.input)

    if not referenceSequences:
        print(f"Error: Failed to load sequences from reference file '{args.reference}'.")
        return

    if not querySequences:
        print(f"Error: Failed to load sequences from query file '{args.input}'.")
        return

    # Realizar o alinhamento e salvar os resultados
    with open(args.outputFile, 'w') as outputFile:
        matches = compareSequences(referenceSequences, querySequences, smithWaterman)
        for match in matches:
            print(f"Passou aqui {args.outputFile}")
            outputFile.write(f"Query {match['queryId']} aligned with Reference {match['refId']}\n")
            outputFile.write(f"Score: {match['score']}\n")
            outputFile.write(f"Aligned Reference: {match['alignedRef']}\n")
            outputFile.write(f"Aligned Query: {match['alignedQuery']}\n\n")

            print(f"Passou aqui também {args.outputFile}") 
    print(f"Alignment results saved to {args.outputFile}")


############################# DATA VIZUALIZING

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
# Example usage
# save_results_to_html(flagCounts)


def executeFlagStats(filePath):
    """
    Execute the analysis of flag statistics on the SAM file.
    
    Args:
    filePath (str): Path to the SAM file to analyze.
    
    Returns:
    dict: A dictionary of flag counts.
    """
    # Call the function to count flags
    flagCounts = countReadsByFlags(filePath)
    # Call the plot function (if necessary)
    plotFlagCounts(flagCounts)
    
    return flagCounts





############################# DATA STORE

def filterSam(filePath, outputFile, minQ=30):
    """
    Filters the SAM file by MAPQ score and writes the filtered reads to a new file.

    :param filePath: Path to the input SAM file.
    :param outputFile: Path to the output SAM file.
    :param minQ: Minimum MAPQ score to retain a read (default is 30).
    """
    with open(filePath, 'r') as infile, open(outputFile, 'w') as outfile:
        for line in infile:
            # Skip header lines 
            if line.startswith('@'):
                outfile.write(line)
                continue
            
            # Split the line into columns
            fields = line.strip().split('\t')
            if len(fields) > 4 and fields[4].isdigit():
                mapq = int(fields[4])  # MAPQ score is in the 5th column (index 4)

            # Filter based on MAPQ score
                if mapq >= minQ:
                    outfile.write(line)


def mappedRead(filePath, outputFile):
    """
    Filters the SAM file and keeping only the mapped reads  to a new file.

    :param filePath: Path to the input SAM file.
    :param outputFile: Path to the output SAM file.
   
    """
    with open(filePath, 'r') as infile, open(outputFile, 'w') as outfile:
        for line in infile:
            # Skip header lines 
            if line.startswith('@'):
                outfile.write(line)
                continue
            
            # Split the line into columns
            fields = line.strip().split('\t')
            flag = int(fields[1])  # 

            # Filter based on MAPQ score
            if flag & 4 == 0 and flag & 8 == 0:
                outfile.write(line)
    

def mappedPrimaryReads(filePath, outputFile):
    """
    Filters the SAM file and keeps only the primary and mapped reads to a new file.

    :param filePath: Path to the input SAM file.
    :param outputFile: Path to the output SAM file.
   
    """
    # Flags to be excluded:
    excludedFlags = [
        'read unmapped',
        'mate unmapped',
        'not primary alignment',
        'read fails platform/vendor quality checks',
        'PCR or optical duplicate',
        'supplementary alignment'
    ]

    with open(filePath, 'r') as infile, open(outputFile, 'w') as outfile:
        for line in infile:
            # Skip header lines 
            if line.startswith('@'):
                outfile.write(line)
                continue
            
            # Split the line into columns
            fields = line.strip().split('\t')
            flag = int(fields[1])  # 

            # Filter based on MAPQ score
            if not any(flag & flags[excludedFlags] for excludedFlags in excludedFlags):
                outfile.write(line)
    



