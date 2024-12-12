
############### IMPORT MODULES ############### 
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os,re,sys,argparse
from collections import defaultdict

from flags import flags


############################# FUNCTIONS TO :

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

            lineCount += 1
            

    print(f"File '{filePath}' has the expected number of columns.")
    print(f"The file: '{filePath}' can be used for next step.")
    return True


# 2/ Read



############################################ READ

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
    flagDetails = {} # Dictionary to count occurrences of each flag description

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

    unmapedReads = (unmapedReads / totalReads) * 100
    duplicatedReads = (duplicatedReads / totalReads) * 100
    mappedReads = (mappedReads / totalReads) * 100
    filteredReads = (filteredReads / totalReads) * 100

    readstatistics = {
        "Unmapped Reads (%)": f"{unmapedReads:.2f}%",
        "Duplicated Reads (%)": f"{duplicatedReads:.2f}%",
        "Mapped Reads (%)": f"{mappedReads:.2f}%"
    }

    # Print the results
    print(f"\n--- Read Statistics ---")
    
    print(f" Percentage of unmapped reads: {unmapedReads:.2f}")
    print(f"Percentage of duplicated reads: {duplicatedReads:.2f}")
    print(f"Percentage of mapped reads: {mappedReads:.2f}")
    print(f"Percentage of filtered reads (MAPQ >= {minQ}): {filteredReads:.2f}")
    print(f"\n")
    for description, count in flagDetails.items():
        print(f"{description}: {count}")
    print(f"\n")

    return readstatistics, flagDetails, totalReads





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

    # Improved readability for printed output
    print("\n--- Reads per Mapping Quality ---")
    print(f"{'Mapping quality ':<10}{'Number of reads':>10}")
    print("-" * 22)
    for mapq, count in sorted(mappingQCount.items()):
        print(f"{mapq:<10}{count:>10}")
    
    return mappingQCount



def countReadsByFlags(filePath):
    """
    Count the number of reads for each flag value in a SAM file.
    
    Args:
    filePath (str): Path to the .sam file.
    
    Returns:
    dict: A dictionary with flag values as keys and counts as values.
    """
    flagCounts = {}

    try:
        with open(filePath, 'r') as file:
            for line in file:
                if line.startswith('@'):  # Skip header lines in SAM file
                    continue
                
                columns = line.split('\t')
                if len(columns) < 11:  # SAM file should have at least 11 columns
                    continue
                
                try:
                    flag = int(columns[1])  # The flag is the second column
                except ValueError:
                    print(f"Warning: Invalid flag value in line: {line.strip()}")
                    continue  # Skip lines where flag is not an integer
                
                # Count occurrences of each flag
                if flag in flagCounts:
                    flagCounts[flag] += 1
                else:
                    flagCounts[flag] = 1
    
    except FileNotFoundError:
        print(f"Error: File not found: {filePath}")
        return {}
    except Exception as e:
        print(f"An error occurred: {e}")
        return {}

    return flagCounts  # Return the dictionary




###Filter Sam file

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
            if flag & 4 == 0:  # if bit 4 is not set, the read is mapped
                outfile.write(line)
    




################################## 3 plots



def plotReadsPerMAPQ(mappingQCount, outputFile="reads_per_mapq.png"):
    """
    Plot the distribution of reads per MAPQ score as a bar chart with intervals of 20 
    and save it as a PNG file.

    Args:
        mappingQCount (dict): A dictionary with MAPQ scores as keys and counts as values.
        outputFile (str): The filename to save the plot.
    """
    try:
        # Group MAPQ scores into intervals of 20
        interval_counts = defaultdict(int)
        for mapq, count in mappingQCount.items():
            interval = (mapq // 20) * 20  # Calculate the lower bound of the interval
            interval_counts[interval] += count
        
        # Sort intervals for plotting
        sorted_intervals = sorted(interval_counts.keys())
        counts = [interval_counts[interval] for interval in sorted_intervals]

        # Create the bar chart
        plt.figure(figsize=(12, 6))
        plt.bar(
            [f"{interval}-{interval + 19}" for interval in sorted_intervals],
            counts,
            color="skyblue",
            edgecolor="black"
        )

        # Add details
        plt.xlabel("MAPQ Score Intervals", fontsize=14, fontweight="bold")
        plt.ylabel("Number of Reads", fontsize=14, fontweight="bold")
        plt.title("Distribution of Reads by MAPQ Score Intervals", fontsize=14, fontweight="bold")
        plt.xticks(rotation=45, fontsize=8)
        plt.yticks(fontsize=8)

        # Ensure the directory exists and save the plot
        os.makedirs(os.path.dirname(outputFile), exist_ok=True)
        plt.savefig(outputFile, format="png", dpi=600)
        print(f"Plot saved at: {outputFile}")
    except Exception as e:
        print(f"Error saving the plot: {e}")
    finally:
        plt.close()

def plotFlagCounts(flagCounts, output_file="flag_counts.png"):
    """
    Plot the flag counts as a bar chart, annotate with percentages, and save the plot as a PNG file.
    
    :param flagCounts: Dictionary containing flag values and their counts.
    :param output_file: Path to save the plot image file.
    :return: DataFrame containing flag counts and percentages for verification.
    """
    # Convert flagCounts dictionary to a DataFrame
    df = pd.DataFrame(list(flagCounts.items()), columns=['Flag', 'Count'])
    
    # Calculate the total number of reads
    total_reads = df['Count'].sum()
    
    # Add a percentage column
    df['Percentage'] = (df['Count'] / total_reads) * 100
    
    # Plot the flag counts as a bar chart
    plt.figure(figsize=(10, 6))
    plt.bar(df['Flag'].astype(str), df['Count'], color='skyblue')
    
    # Annotate percentages on the bars
    for i, (count, percentage) in enumerate(zip(df['Count'], df['Percentage'])):
        plt.text(i, count, f"{percentage:.1f}%", ha='center', va='bottom', fontsize=8)
    
    # Customize plot
    plt.xlabel('Flag Value')
    plt.ylabel('Read Count')
    plt.title('Number of Reads for Each Flag Value')
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    # Save the plot to a file
    plt.savefig(output_file, format='png', dpi=300)
    print(f"Plot successfully saved to {output_file}")
    
    
    
    # Return the updated DataFrame for verification if needed
    return df


def plotReadsPercentage(readstatistics, outputFile="read_distribution.png"):
    """
    Plot the read distribution as percentages in a pie chart and save it as a PNG file.

    Args:
        readstatistics (dict): A dictionary containing read percentages, potentially nested.
        outputFile (str): The filename to save the plot.
    """
    try:
        # Si le dictionnaire est imbriqué, récupérer la première valeur (supposée être le sous-dictionnaire)
        if any(isinstance(value, dict) for value in readstatistics.values()):
            readstatistics = list(readstatistics.values())[0]  # Récupérer le sous-dictionnaire

        # Validation des données
        if not isinstance(readstatistics, dict):
            raise ValueError("Input readstatistics must be a dictionary of percentages.")

        labels = list(readstatistics.keys())
        sizes = []
        for value in readstatistics.values():
            if isinstance(value, str) and value.endswith('%'):
                sizes.append(float(value.strip('%')))
            elif isinstance(value, (int, float)):
                sizes.append(float(value))
            else:
                raise ValueError(f"Invalid value format: {value}. Expected percentage strings or numeric values.")

        # Préparation des couleurs et de l'explosion
        colors = ['lightblue', 'orange', 'lightgreen', 'red', 'purple'][:len(labels)]
        explode = [0.1 if i == 0 else 0.05 for i in range(len(labels))]

        # Création du graphique
        plt.figure(figsize=(10, 8))
        wedges, texts, autotexts = plt.pie(
            sizes,
            labels=labels,
            autopct='%1.2f%%',
            colors=colors,
            startangle=90,
            explode=explode,
            wedgeprops={'edgecolor': 'black', 'linewidth': 1.2}
        )

        # Titre
        plt.title('Read Distribution (%)', fontsize=16, fontweight='bold')
        plt.axis('equal')  # Maintenir un cercle parfait

        # Vérification du répertoire et sauvegarde
        os.makedirs(os.path.dirname(outputFile), exist_ok=True)
        plt.savefig(outputFile, format='png', dpi=300)
        print(f"Plot saved at: {outputFile}")
    except ValueError as ve:
        print(f"ValueError during plotting: {ve}")
    except Exception as e:
        print(f"Error saving the plot: {e}")
    finally:
        plt.close()

def summaryTable(filePath, output_file="flag_summary.txt"):
    """
    Create a summary table based on flag distribution from a SAM file and save it to a file.
    
    Args:
    filePath (str): Path to the .sam file.
    output_file (str): Path to save the generated text file.
    """
    # Predefined flag descriptions for known flag values (SAM flag values)
    flagDescriptions = {
         0: 'Unmapped, unpaired read',
        1: 'Read is paired',
        2: 'Both reads mapped in proper pair',
        4: 'Read unmapped',
        8: 'Mate unmapped',
        16: 'Read reverse strand',
        32: 'Mate reverse strand',
        64: 'First in pair',
        128: 'Second in pair',
        256: 'Not primary alignment',
        512: 'Failed quality check',
        1024: 'Duplicate read',
        2048: 'Supplementary alignment',
        77: 'Read unmapped',
        99: 'Paired read',
        145: 'Second in pair',
        147: 'Properly paired',
        163: 'First in pair',
        97: 'Read reverse strand',
        141: 'Mate unmapped',
        83: 'Mate reverse strand',
        
    }
    
    # Get the flag counts from the SAM file
    flagCounts = countReadsByFlags(filePath)

    # Prepare the table headers
    table_content = f"""
    Flag Distribution Summary
    {'-'*80}
    {'Flag':<20} {'Description':<30} {'Read Count':<15} {'Percentage (%)'}
    {'-'*80}
    """

    # Add rows for each flag in flagCounts
    totalReads = sum(flagCounts.values())  # Total number of reads (sum of all flag counts)
    for flag, count in flagCounts.items():
        description = flagDescriptions.get(flag, 'Unknown')  # Get description for the flag
        percentage = (count / totalReads) * 100  # Calculate percentage
        table_content += f"{flag:<20} {description:<30} {count:<15} {percentage:.2f}%\n"

    # Save the plain text table to the output file
    with open(output_file, 'w') as file:
        file.write(table_content)
    
    print(f"Flag summary table saved to {output_file}")

      

# Function to execute flag stats and generate plots
def executePlots(filePath, minQ=0):
    """
    Execute the analysis of flag statistics and read statistics on the SAM file,
    generating a bar chart, pie chart, and MAPQ distribution, and generating a flag summary file.
    
    Args:
        filePath (str): Path to the SAM file to analyze.
        minQ (int): Minimum MAPQ score for filtering reads. Defaults to 0.
    
    Returns:
        tuple: Paths to the generated plot images, read statistics, and summary flag table file.
    """
    # Count read statistics and generate percentages
    readstatistics = countReads(filePath, minQ=minQ)
    
    # Convert readstatistics to a dictionary if it's in a tuple form
    if isinstance(readstatistics, tuple):
        readstatistics = {
            "Mapped Reads (%)": readstatistics[0],
            "Unmapped Reads (%)": readstatistics[1],
            # Add more metrics from the tuple if necessary
        }
    
    # Count occurrences of each FLAG
    flagCounts = countReadsByFlags(filePath)
    totalReads = sum(flagCounts.values())  # Total reads based on flag counts
    
    # Generate the MAPQ counts
    mappingQCount = readPerMAPQ(filePath)
    
    # Define file paths for plots
    flag_plot_path = '/home/najat/mapping/src/flag_counts_plot.png'
    pie_plot_path = '/home/najat/mapping/src/read_distribution.png'
    mappingq_plot_path = '/home/najat/mapping/src/reads_per_mapq.png'
    summary_flag_table_path = '/home/najat/mapping/src/flag_summary.txt'

    # Generate and save the bar chart for flag counts
    df_flag = plotFlagCounts(flagCounts)  # Returns DataFrame with percentages
    
    # Save flag counts plot
    plt.figure(figsize=(10, 6))
    plt.bar(df_flag['Flag'].astype(str), df_flag['Count'], color='skyblue')
    for i, (count, percentage) in enumerate(zip(df_flag['Count'], df_flag['Percentage'])):
        plt.text(i, count, f"{percentage:.1f}%", ha='center', va='bottom', fontsize=8)
    plt.xlabel('Flag Value')
    plt.ylabel('Read Count')
    plt.title('Number of Reads for Each Flag Value (with Percentages)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(flag_plot_path)
    plt.close()
    
    # Generate the pie chart for read distribution
    plotReadsPercentage(readstatistics, pie_plot_path)
    
    # Generate the bar chart for MAPQ distribution
    plotReadsPerMAPQ(mappingQCount, mappingq_plot_path)
    
    # Call the summaryTable function to generate and save the flag summary table
    summaryTable(filePath, output_file=summary_flag_table_path)
    
    # Return all plot paths, statistics, and summary flag table path
    return flag_plot_path, pie_plot_path, mappingq_plot_path, readstatistics, mappingQCount, summary_flag_table_path

def saveResults(plot_paths, flagDetails, summary_flag_table_path, html_output_path="/home/najat/mapping/src/Result.html"):
    """
    Create an HTML file with multiple plots embedded in it and a flag distribution summary from an external text file.
    
    :param plot_paths: List of plot image file paths.
    :param flagDetails: Dictionary containing counts for each flag.
    :param summary_flag_table_path: Path to the summary flag table file (e.g., summaryflagtable.txt).
    :param html_output_path: Path to save the generated HTML file.
    """
    # Read the summary flag table from the provided file
    try:
        with open(summary_flag_table_path, 'r') as f:
            summary_flag_table_content = f.read()
    except FileNotFoundError:
        print(f"Error: The file {summary_flag_table_path} was not found.")
        return

    # Generate the HTML content with improved styling
    html_content = f"""
    <html>
    <head>
        <title>Sequencing Data Analysis Report</title>
        <style>
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                color: #333;
                background-color: #f9f9f9;
                margin: 0;
                padding: 0;
            }}
            h1, h2 {{
                color: #0056b3;
                text-align: center;
            }}
            .container {{
                width: 80%;
                margin: 0 auto;
                padding: 20px;
            }}
            .section {{
                background-color: #fff;
                padding: 20px;
                margin-bottom: 20px;
                border-radius: 8px;
                box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.1);
            }}
            .section-title {{
                font-size: 24px;
                font-weight: bold;
                margin-bottom: 15px;
                color: #333;
            }}
            <pre {{
                font-family: Courier, monospace;
                background-color: #f1f1f1;
                padding: 20px;
                border-radius: 8px;
                white-space: pre;  /* Cela permet de conserver la structure de ligne du texte sans ajouter d'espaces supplémentaires */
                word-wrap: break-word;  /* Permet de couper les mots longs sans dépasser la largeur de l'écran */
                overflow-x: auto;
                margin-top: 0;  /* Enlever l'espace supérieur */
                }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin-top: 20px;
            }}
            table, th, td {{
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }}
            th {{
                background-color: #0056b3;
                color: #fff;
            }}
            tr:nth-child(even) {{
                background-color: #f2f2f2;
            }}
            img {{
                width: 48%;
                margin: 10px;
                border-radius: 8px;
                box-shadow: 0px 2px 5px rgba(0, 0, 0, 0.1);
            }}
            .plot-container {{
                display: flex;
                justify-content: space-between;
                flex-wrap: wrap;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>Sequencing Data Analysis</h1>
            
            <!-- Flag Distribution Summary Section -->
            <div class="section">
                <div class="section-title">Flag Distribution Summary</div>
                <pre>{summary_flag_table_content}</pre>
            </div>

            <!-- Plots Section -->
            <div class="section">
                <div class="section-title">Plots</div>
                <div class="plot-container">
    """
    
    # Add each plot image
    for plot_path in plot_paths:
        html_content += f'<img src="{plot_path}" alt="Plot Image">'

    # Finalize HTML
    html_content += """
                </div>
            </div>
        </div>
    </body>
    </html>
    """

    # Write the HTML content to a file
    with open(html_output_path, 'w') as f:
        f.write(html_content)
    
    print(f"HTML report saved to {html_output_path}")
    
    
    
    
    
    ###################Mapping 
    
    
def parseSam(filePath):
        """
        Extract sequences from a SAM file.

        :param filePath: Path to the SAM file.
        :return: List of sequences (each as a tuple with read ID and sequence).
        """
        print("Entered parseSam") 
        sequences = []
        with open(filePath, 'r') as file:
            for line in file:
                if line.startswith('@'):  # Skip header lines
                    continue
                parts = line.strip().split('\t')
                if len(parts) > 9:  # Ensure there are at least 10 columns
                    readId = parts[0]   # Column 1: Read ID
                    sequence = parts[9]  # Column 10: Sequence
                    sequences.append((readId, sequence))
        print(f"Loaded {len(sequences)} sequences from the file.")
        return sequences



def smithWaterman(sequences: list,  reference: str, matchScore: int = 2, mismatchPenalty: int = -2, gapPenalty: int = -3):
    """
    Perform the Smith-Waterman alignment algorithm.

    :param sequences: List of tuples containing read ID and query sequences.
    :param reference: Reference sequence to align against.
    :param outputFilePath: Path to save the alignment results.
    :param matchScore: Score for matching characters.
    :param mismatchPenalty: Penalty for mismatched characters.
    :param gapPenalty: Penalty for gaps.
    """
  
    m, n = len(sequences), len(reference)
    scoreMatrix = np.zeros((m + 1, n + 1), dtype=int)
    maxScore = 0
    maxPos = None

    # Preenchendo a matriz de pontuação com penalidades ajustadas
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = scoreMatrix[i - 1][j - 1] + (matchScore if sequences[i - 1] == reference[j - 1] else mismatchPenalty)
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
        if scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + (matchScore if sequences[i - 1] == reference[j - 1] else mismatchPenalty):
            alignedSeqOne.append(sequences[i - 1])
            alignedSeqTwo.append(reference[j - 1])
            i -= 1
            j -= 1
        elif scoreMatrix[i][j] == scoreMatrix[i - 1][j] + gapPenalty:
            alignedSeqOne.append(sequences[i - 1])
            alignedSeqTwo.append('-')
            i -= 1
        else:
            alignedSeqOne.append('-')
            alignedSeqTwo.append(reference[j - 1])
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

