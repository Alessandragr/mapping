#!/usr/bin/python3
#-*- coding : utf-8 -*-

############### IMPORT MODULES ############### 
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os,re,sys,argparse
from collections import defaultdict

from flags import flags





# OPTION LIST:
# -h or --help: Displays help information.
# -i or --input: Specifies the path to the input SAM file (.sam).
# -o or --output: Specifies the path to the output file 
# -cR or --countReads: Counts total reads, mapped reads, unmapped reads,duplicated reads, with a filtering option on the mapping quality.
# -rC or --readPerChrom: Counts the number of reads per chromosome.
# -rMQ or --readPerMAPQ: Counts the number of reads for each MAPQ score.
# -cRF or --countReadsByFlags: Counts the number of reads for each FLAG value.
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
#-mr -eFS


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


def smithWaterman(seqOne: str, seqTwo: str, matchScore: int = 2, mismatchPenalty: int = -2, gapPenalty: int = -3):
    """
    Perform the Smith-Waterman alignment algorithm.

    :param seqOne: First sequence (query).
    :param seqTwo: Second sequence (reference).
    :param matchScore: Score for matching characters.
    :param mismatchPenalty: Penalty for mismatched characters.
    :param gapPenalty: Penalty for gaps.
    :return: Tuple containing the aligned sequences and the maximum alignment score.
    """
    m, n = len(seqOne), len(seqTwo)
    scoreMatrix = np.zeros((m + 1, n + 1), dtype=int)
    maxScore = 0
    maxPos = None

    # Fill the scoring matrix with adjusted penalties
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = scoreMatrix[i - 1][j - 1] + (matchScore if seqOne[i - 1] == seqTwo[j - 1] else mismatchPenalty)
            delete = scoreMatrix[i - 1][j] + gapPenalty
            insert = scoreMatrix[i][j - 1] + gapPenalty
            scoreMatrix[i][j] = max(0, match, delete, insert)

            if scoreMatrix[i][j] > maxScore:
                maxScore = scoreMatrix[i][j]
                maxPos = (i, j)

    # Traceback to construct the alignment
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




    


# Main script
def main():
    parser = argparse.ArgumentParser(description="Analyze a SAM file and provide various statistics.")

    # Arguments for file input and output
    parser.add_argument('-i', '--input', required=True, help="Path to the SAM file.")
    parser.add_argument('-o', '--outputFile', required=False, help="Path to the output SAM file.")

    # Arguments for read counting and analysis
    parser.add_argument('-cR', '--countReads', action='store_true', help="Count reads (total, mapped, etc.).")
    parser.add_argument('-rC', '--readPerChrom', action='store_true', help="Count reads per chromosome.")
    parser.add_argument('-rMQ', '--readPerMAPQ', action='store_true', help="Count reads based on MAPQ scores.")
    parser.add_argument('-cRF', '--countReadsByFlags', action='store_true', help="Count reads based on FLAG values.")

    # Arguments for saving results and generating output
    parser.add_argument('-sR', '--saveResults', action='store_true', help="Save results and graphs to an HTML file.")

    # Arguments for filtering and saving new SAM files
    parser.add_argument('-fS', '--filterSam', action='store_true', help="Filter SAM file by mapping quality")
    parser.add_argument('-mR', '--mappedRead', action='store_true', help="Filter SAM file by keeping mapped reads")
    parser.add_argument('-r', '--reference', help="Reference sequence for alignment.")
    parser.add_argument('-q', '--query', help="Query sequence for alignment.")

    # Argument for MAPQ score filtering
    parser.add_argument('-m', '--minQ', type=int, default=0, help="Minimum MAPQ score for filtering reads.")

    # Argument for executing flag stats (calls executePlots)
    parser.add_argument('-eP', '--executePlots', action='store_true', help="Execute Flag Stats analysis.")

    # Arguments for Smith-Waterman alignment
    parser.add_argument('-sw', '--smithWaterman', action='store_true', help="Perform Smith-Waterman alignment.")

    args = parser.parse_args()

    # Display help if no arguments are provided
    if not any(vars(args).values()):
        parser.print_help()
        return

    # Verify the input file (this function should be defined elsewhere in your code)
    if not fileVerifier(args.input):
        print(f"Error: The input file {args.input} does not exist or cannot be opened.")
        return

    # Execute the appropriate functions based on the arguments
    if args.smithWaterman:
        if not args.reference or not args.query:
            print("Error: Both reference and query sequences are required for Smith-Waterman alignment.")
            return

        alignment_result = smithWaterman(args.reference, args.query)
        print("Smith-Waterman Alignment Result:")
        print(f"Aligned Sequence 1: {alignment_result[0]}")
        print(f"Aligned Sequence 2: {alignment_result[1]}")
        print(f"Alignment Score: {alignment_result[2]}")

    # Existing code for read counting, flag stats, and result saving
    if args.countReads:
        readstatistics = countReads(args.input, args.minQ)

    if args.readPerChrom:
        readPerChrom(args.input)

    if args.readPerMAPQ:
        readPerMAPQ(args.input)

    if args.countReadsByFlags:
        flagCounts = countReadsByFlags(args.input)
        flagDetails = {key: count for key, count in flagCounts.items()}

    if args.filterSam:
        filterSam(args.input, args.outputFile, args.minQ)

    if args.mappedRead:
        mappedRead(args.input, args.outputFile)

    if args.executePlots:
        # Execute flag stats and retrieve results
        flag_plot_path, pie_plot_path, mappingq_plot_path, readstatistics, mappingQCount, summary_flag_table_path = executePlots(args.input, args.minQ)

        # Collect the generated plot paths
        plot_paths = [flag_plot_path, pie_plot_path, mappingq_plot_path]

        # Update flag counts based on the SAM file
        flagCounts = countReadsByFlags(args.input)

        # Create flag details dictionary, ensuring it contains flag counts and descriptions
        flagDetails = {flag: {'count': count, 'description': flags.get(flag, 'Unknown flag')} for flag, count in flagCounts.items()}

        print(f"Flag Details: {flagDetails}")
        print(f"Read Statistics: {readstatistics}")
        print(f"Mapping Quality Counts: {mappingQCount}")

    # Save results to HTML if required
    if args.saveResults:
        if not plot_paths:
            print("Error: Plot paths are empty. Ensure executePlots is run before saving results.")
            return

        if not flagCounts:
            print("Error: Flag counts are empty. Ensure countReadsByFlags returns valid results.")
            return

        saveResults(
            plot_paths=plot_paths,
            flagDetails=flagDetails,
            summary_flag_table_path=summary_flag_table_path
        )

if __name__ == "__main__":
    main()


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



