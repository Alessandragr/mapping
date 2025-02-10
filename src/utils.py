
############### IMPORT MODULES ############### 
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import re
from collections import defaultdict
import datetime, shutil
import math
from tabulate import tabulate 

from flags import flags



############################# FUNCTIONS TO :

############################# CHECK 
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

    # Validate that the number of columns is at least 11
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




############################# READS STAT
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
                    if flag not in flagDetails:
                        flagDetails[flag] = {'count': 0, 'description': description}
                    flagDetails[flag]['count'] += 1

            # Check flags for specific categories
            if flag & 4 or flag & 8:  # Unmapped reads
                unmapedReads += 1
            elif flag & 1024:  # Duplicated reads
                duplicatedReads += 1
            else:
                mappedReads += 1  # Mapped reads
                mapq = int(fields[4])  # MAPQ score is in the 5th column

                # Apply the MAPQ filter
                if mapq >= minQ:
                    filteredReads += 1

    unmapedReadsPercentage = (unmapedReads / totalReads) * 100
    duplicatedReadsPercentage = (duplicatedReads / totalReads) * 100
    mappedReadsPercentage = (mappedReads / totalReads) * 100
    filteredReadsPercentage = (filteredReads / totalReads) * 100

    readstatistics = {
        "Unmapped Reads (%)": f"{unmapedReadsPercentage:.2f}%",
        "Duplicated Reads (%)": f"{duplicatedReadsPercentage:.2f}%",
        "Mapped Reads (%)": f"{mappedReadsPercentage:.2f}%"
    }

    # Stats showed on the terminal to the user
    print("\n-------------------------- Reads Statistics --------------------------\n")

    print(f"{'Description':<40} {'Percentage (%)':<20} {'Total':<18}")
    print(f"{'--' * 35}\n")


    print(f"{'Unmapped Reads':<40} {unmapedReadsPercentage:.2f}%  {unmapedReads:>18.0f}")
    print(f"{'Duplicated Reads':<40} {duplicatedReadsPercentage:.2f}%  {duplicatedReads:>18.0f}")
    print(f"{'Mapped Reads':<40} {mappedReadsPercentage:.2f}%  {mappedReads:>18.0f}")
    print(f"{'Filtered Reads (MAPQ >= {minQ})':<40} {filteredReadsPercentage:.2f}%  {filteredReads:>18.0f}")

    print(f"\n{'--' * 35}\n")

   
    print(f"{'Total Reads:':<40} {totalReads:>26}")
    print(f"{'--' * 35}\n\n\n")


    # Display Flag Details Table
    print(f"{'--' * 20}- Flag Details -{'--' * 20} \n")
    print(f"{'Flag':<10} {'Description':<55} {'Percentage (%)':<22} {'Total':<34}")
    print(f"{'--' * 48}\n")

    for key, details in sorted(flagDetails.items()):
        count = details['count']
        description = details['description']
        percentage = (count / totalReads) * 100 if totalReads > 0 else 0
        print(f"{key:<10} {description:<40} {percentage:>15.2f}%  {count:>18}")


    print(f"\n{'--' * 48}\n")


    return readstatistics, flagDetails, totalReads




def readPerChrom(filePath):
    """
    Count the number of reads mapped to each chromosome. Also show some stats on the terminal: mean and standard deviation
    
    :param filePath: path to the SAM file
    """
    chromosomeCounts = {}
    totalReads = 0
    chromosomeMAPQSum = {}  # Sum the values of MAPQ per chromosome
    chromosomeMAPQCount = {}  # Count the number of reads per chromosome
    chromosomeMAPQSquaresSum = {}  # Sum the squares of the values of MAPQ per chromosome

    with open(filePath, 'r') as file:
        for line in file:
            # Ignore headers
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            chromosome = fields[2]
            mapq = int(fields[4])  # MAPQ score

            flag = int(fields[1])
            # Only count reads where both the read and its mate are mapped
            if flag & 4 == 0 and flag & 8 == 0:  # Exclude unmapped reads and mates
                chromosomeCounts[chromosome] = chromosomeCounts.get(chromosome, 0) + 1
                chromosomeMAPQSum[chromosome] = chromosomeMAPQSum.get(chromosome, 0) + mapq
                chromosomeMAPQCount[chromosome] = chromosomeMAPQCount.get(chromosome, 0) + 1
                chromosomeMAPQSquaresSum[chromosome] = chromosomeMAPQSquaresSum.get(chromosome, 0) + mapq**2
                totalReads += 1
                
    # Stats showed on the terminal to the user:
    print(f"\n{'--' * 20}  Reads per Chromosome -{'--' * 20}\n")

    print(f"{'Chromosome':<20} {'Percentage (%)':<18} {'Total':<18} {'Mean MAPQ':<18} {'Standard deviation MAPQ':<18}")
    print(f"{'--' * 52}\n")

    for chrom in sorted(chromosomeCounts.keys()):
        count = chromosomeCounts[chrom]
        percentage = (count / totalReads) * 100
        avgMAPQ = chromosomeMAPQSum[chrom] / chromosomeMAPQCount[chrom] if chromosomeMAPQCount[chrom] > 0 else 0

        # Standard deviation
        variance = (chromosomeMAPQSquaresSum[chrom] / chromosomeMAPQCount[chrom]) - (avgMAPQ ** 2)
        stddev = math.sqrt(variance) if variance > 0 else 0
        print(f"{chrom:<6} {percentage:>18.2f}% {count:>15} {avgMAPQ:>20.2f} {stddev:>25.2f}")


    print(f"\n{'--' * 52}\n")
    print(f"{'Total Reads:':<17} {totalReads:>27}")
    print(f"{'--' * 52}\n")




def readPerMAPQ(filePath):
    """
    Count the number of reads for each MAPQ score.
    
    :param filePath: path to the SAM file
    :return: dictionary with counts of reads per MAPQ score
    """
    mappingQCount = {}
    totalReads = 0

    with open(filePath, 'r') as file:
        for line in file:
            # Ignore headers
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            mapq = int(fields[4])

            flag = int(fields[1])
            if flag & 4 == 0 and flag & 8 == 0:  # if bit 4 and 8 are not set, the read is mapped
                mappingQCount[mapq] = mappingQCount.get(mapq, 0) + 1
                totalReads += 1

    print("\n\n\n--------------------- Reads per Mapping Quality ----------------------\n")

    print(f"{'Mapping quality ':<30}{'Percentage (%)':<30} {'Total':<18}")

    print(f"\n{'--' * 35}\n")

    for mapq, count in sorted(mappingQCount.items()):
        percentage = (count / totalReads) * 100
        print(f"{mapq:<10} {percentage:>26.2f}% {count:>27}")

    print(f"\n{'--' * 35}\n")

    print(f"{'Total Reads:':<39} {totalReads:>26}")
    print(f"{'--' * 35}\n")

    return mappingQCount




def countReadsByFlags(filePath):
    """
    Count the number of reads for each flag value in a SAM file.
    
    :param filePath: Path to the .sam file.
    
    :returns: A dictionary with flag values as keys and counts as values.
    """
    flagCounts = {}

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
    print("The flag count was executed succefully!")
    return flagCounts  # Return the dictionary
    



############################# FILTER FILES
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
    print(" Filtered sam file was created succesfully")




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
            if flag & 4 == 0 and flag & 8 == 0:  # if bit 4 is not set, the read is mapped
                outfile.write(line)
    print(" Filtered sam file was created succesfully")            
    
    

    

def mappedPrimaryReads(filePath, outputFile):
    """
    Filters the SAM file and keeps only the primary and mapped reads to a new file.

    :param filePath: Path to the input SAM file.
    :param outputFile: Path to the output SAM file.
    """
    # Flags to be excluded (using numerical values from the `flags` dictionary):
    excludedFlagValues = [
        4,    # read unmapped
        8,    # mate unmapped
        256,  # not primary alignment
        512,  # read fails platform/vendor quality checks
        1024, # PCR or optical duplicate
        2048  # supplementary alignment
    ]

    try:
        with open(filePath, 'r') as infile, open(outputFile, 'w') as outfile:
            for line in infile:
                # Write header lines directly to the output
                if line.startswith('@'):
                    outfile.write(line)
                    continue
                
                # Process alignment lines
                fields = line.strip().split('\t')
                flag = int(fields[1])  # Extract flag field (column 2 in SAM format)

                # Exclude reads with any of the excluded flags
                if not any(flag & excludedFlag for excludedFlag in excludedFlagValues):
                    outfile.write(line)

        print("Filtered SAM file created successfully.")
    
    except FileNotFoundError:
        print(f"Error: The file {filePath} was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
     

    

############################# DATA VIZUALIZATION
def plotReadsPerMAPQ(mappingQCount, outputFile="reads_per_mapq.png"):
    """
    Plot the distribution of reads per MAPQ score as a bar chart with intervals of 20 
    and save it as a PNG file.

    :param mappingQCount: A dictionary with MAPQ scores as keys and counts as values.
    :param outputFile: The filename to save the plot.
    :return: The path to the saved plot.
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
        plt.xticks(rotation=45, fontsize=10)
        plt.yticks(fontsize=10)

        # Ensure the directory exists and save the plot
        os.makedirs(os.path.dirname(outputFile), exist_ok=True)
        plt.savefig(outputFile, format="png", dpi=600)
        print(f"Plot saved at: {outputFile}")
    except Exception as e:
        print(f"Error saving the plot: {e}")
    finally:
        plt.close()
    
    return outputFile





def plotFlagCounts(flagCounts, outputFile="flag_counts.png"):
    """
    Plot the flag counts as a bar chart, annotate with percentages, and save the plot as a PNG file.

    :param flagCounts: Dictionary containing flag values and their counts.
    :param outputFile: Path to save the plot image file.
    :return: DataFrame containing flag counts and percentages for verification.
    """
    try:
        # Convert flagCounts dictionary to a DataFrame
        dataFrame = pd.DataFrame(list(flagCounts.items()), columns=['Flag', 'Count'])

        # Calculate the total number of reads
        total_reads = dataFrame['Count'].sum()

        # Add a percentage column
        dataFrame['Percentage'] = (dataFrame['Count'] / total_reads) * 100

        # Plot the flag counts as a bar chart
        plt.figure(figsize=(10, 6))
        plt.bar(dataFrame['Flag'].astype(str), dataFrame['Count'], color='skyblue')

        # Annotate percentages on the bars
        for i, (count, percentage) in enumerate(zip(dataFrame['Count'], dataFrame['Percentage'])):
            plt.text(i, count, f"{percentage:.1f}%", ha='center', va='bottom', fontsize=8)

        # Customize plot
        plt.xlabel('Flag Value', fontsize=12)
        plt.ylabel('Read Count', fontsize=12)
        plt.title('Number of Reads for Each Flag Value', fontsize=14)
        plt.xticks(rotation=45, fontsize=10)
        plt.tight_layout()

        # Save the plot to a file
        plt.savefig(outputFile, format='png', dpi=300)
        print(f"Plot successfully saved to {outputFile}")
    except Exception as e:
        print(f"Error creating the flag counts plot: {e}")
    finally:
        plt.close()

    return dataFrame





def plotReadsPercentage(readstatistics, outputFile="read_distribution.png"):
    """
    Plot the read distribution as percentages in a pie chart and save it as a PNG file.

    :param readstatistics: A dictionary containing read percentages, potentially nested.
    :param outputFile: The filename to save the plot.
    """
    try:
        # Handle nested dictionaries by extracting the first level
        if any(isinstance(value, dict) for value in readstatistics.values()):
            readstatistics = list(readstatistics.values())[0]

        # Validate that the input is a dictionary
        if not isinstance(readstatistics, dict):
            raise ValueError("Input readstatistics must be a dictionary of percentages.")

        # Extract labels and sizes
        labels = list(readstatistics.keys())
        sizes = []
        for value in readstatistics.values():
            if isinstance(value, str) and value.endswith('%'):
                sizes.append(float(value.strip('%')))
            elif isinstance(value, (int, float)):
                sizes.append(float(value))
            else:
                raise ValueError(f"Invalid value format: {value}. Expected percentage strings or numeric values.")

        # Define colors and explode effect
        colors = ['lightblue', 'orange', 'lightgreen', 'red', 'purple'][:len(labels)]
        explode = [0.1 if i == 0 else 0.05 for i in range(len(labels))]

        # Create the pie chart
        plt.figure(figsize=(10, 8))
        plt.pie(
            sizes,
            labels=labels,
            autopct='%1.2f%%',
            colors=colors,
            startangle=90,
            explode=explode,
            wedgeprops={'edgecolor': 'black', 'linewidth': 1.2}
        )

        # Add title and save the plot
        plt.title('Read Distribution (%)', fontsize=16, fontweight='bold')
        plt.axis('equal')  # Ensure a perfect circle

        os.makedirs(os.path.dirname(outputFile), exist_ok=True)
        plt.savefig(outputFile, format='png', dpi=300)
        print(f"Plot saved at: {outputFile}")
    except ValueError as ve:
        print(f"ValueError during plotting: {ve}")
    except Exception as e:
        print(f"Error saving the plot: {e}")
    finally:
        plt.close()




def decodeFlag(flag):
    """
    Decode a SAM flag into its constituent components using bitwise operations.
    :param flag (int): The SAM flag to decode.
    :returns (str): A human-readable description of the flag.
    """
    descriptions = [description for bit, description in flags.items() if flag & bit]
    return ", ".join(descriptions) if descriptions else "Unknown"





def summaryTable(filePath, outputFile="flag_summary.txt"):
    """
    Create a summary table based on flag distribution from a SAM file and save it to a file.

    :param filePath: Path to the .sam file.
    :param outputFile: Path to save the generated text file.
    """
    try:
        # Get the flag counts from the SAM file
        flagCounts = countReadsByFlags(filePath)

        # Prepare table content
        totalReads = sum(flagCounts.values())  # Total number of reads
        table_data = []
        for flag, count in flagCounts.items():
            description = decodeFlag(flag)
            percentage = (count / totalReads) * 100 if totalReads > 0 else 0
            table_data.append([flag, description, count, f"{percentage:.2f}%"])

        # Format table using tabulate
        headers = ["Flag", "Description", "Read Count", "Percentage (%)"]
        table_content = tabulate(table_data, headers=headers, tablefmt="grid", numalign="right", stralign="left")

        # Save the table to the output file
        with open(outputFile, 'w') as file:
            file.write(table_content)

        print(f"Flag summary table saved to {outputFile}")
    except Exception as e:
        print(f"Error creating the summary table: {e}")





############################# EXECUTE FLAGS AND CREATE PLOTS
def executePlots(filePath, outputDir="result", outputFile="result", minQ=0):
    """
    Execute the analysis of flag statistics and read statistics on the SAM file,
    generating plots and a summary flag file.
    
    :param filePath (str): Path to the SAM file to analyze.
    :param outputDir (str): Directory where output files will be saved. Defaults to "result".
    :param outputFile (str): Prefix for output files. Defaults to "result".
    :param minQ (int): Minimum MAPQ score for filtering reads. Defaults to 0.
    
    :returns: dict: Paths to the generated plot images, read statistics, and summary flag table file.
    """
    try:
        # Resolve the output directory to an absolute path
        outputDir = os.path.abspath(outputDir)
        
        # Ensure the output directory exists (delete and recreate if it does)
        if os.path.exists(outputDir):
            shutil.rmtree(outputDir)
        os.makedirs(outputDir, exist_ok=True)

        print(f"Results will be stored in: {outputDir}")

        # Generate a timestamp for unique filenames
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')  
        outputFile = f"{outputFile}_{timestamp}"

        # Count read statistics
        readStatistics = countReads(filePath, minQ=minQ)
        if isinstance(readStatistics, tuple):  # Handle tuple output for compatibility
            readStatistics = {
                "Mapped Reads (%)": readStatistics[0],
                "Unmapped Reads (%)": readStatistics[1],
            }

        # Count occurrences of each FLAG
        flagCounts = countReadsByFlags(filePath)

        # Generate MAPQ counts
        mappingQCount = readPerMAPQ(filePath)

        # Define output file paths
        mappingq_plot_path = os.path.join(outputDir, "mapq_distribution.png")
        pie_plot_path = os.path.join(outputDir, "read_distribution.png")
        flag_plot_path = os.path.join(outputDir, "flag_counts.png")
        summary_flag_table_path = os.path.join(outputDir, "flag_summary.txt")

        # Generate plots and summary files
        try:
            plotFlagCounts(flagCounts, flag_plot_path)
            print(f"Flag counts plot saved at: {flag_plot_path}")
        except Exception as e:
            print(f"Error creating the flag counts plot: {e}")

        try:
            plotReadsPercentage(readStatistics, pie_plot_path)
            print(f"Read distribution plot saved at: {pie_plot_path}")
        except Exception as e:
            print(f"Error creating the read distribution plot: {e}")

        try:
            plotReadsPerMAPQ(mappingQCount, mappingq_plot_path)
            print(f"MAPQ distribution plot saved at: {mappingq_plot_path}")
        except Exception as e:
            print(f"Error creating the MAPQ distribution plot: {e}")

        try:
            summaryTable(filePath, summary_flag_table_path)
            print(f"Flag summary table saved to: {summary_flag_table_path}")
        except Exception as e:
            print(f"Error creating the flag summary table: {e}")

        # Return paths of generated files
        plot_paths = {
            'flag_plot': flag_plot_path,
            'pie_plot': pie_plot_path,
            'mappingq_plot': mappingq_plot_path,
            'readStatistics': readStatistics,
            'mappingQCount': mappingQCount,
            'summary_flag_table': summary_flag_table_path
        }

        return plot_paths

    except Exception as e:
        print(f"Error during analysis and plotting: {e}")
        raise




def saveResults(plot_paths, summary_flag_table_path, final_cigar_table_path="Final_Cigar_table.txt", html_output_path="Result.html"):
    """
    Create an HTML file with multiple plots embedded in it and a flag distribution summary from an external text file.
    
    :param plot_paths: Dictionary containing paths to plot image files and summary flag table.
    :param summary_flag_table_path: Path to the summary flag table text file.
    :param final_cigar_table_path: Path to the final CIGAR table text file (optional).
    :param html_output_path: Path to save the generated HTML file.
    """
    # Read the summary flag table content (which is a .txt file)
    try:
        with open(summary_flag_table_path, 'r') as f:
            summary_flag_table_content = f.read()
    except FileNotFoundError:
        print(f"Error: The file {summary_flag_table_path} was not found.")
        return
    except Exception as e:
        print(f"Error reading the file {summary_flag_table_path}: {e}")
        return

    # Read the final CIGAR table content (optional)
    final_cigar_table_content = ""
    if final_cigar_table_path and os.path.exists(final_cigar_table_path):
        try:
            with open(final_cigar_table_path, 'r') as f:
                final_cigar_table_content = f.read()
        except FileNotFoundError:
            print(f"Error: The file {final_cigar_table_path} was not found.")
        except Exception as e:
            print(f"Error reading the file {final_cigar_table_path}: {e}")

    # Generate the HTML content
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
            pre {{
                font-family: Courier, monospace;
                background-color: #f1f1f1;
                padding: 20px;
                border-radius: 8px;
                white-space: pre;
                word-wrap: break-word;
                overflow-x: auto;
                margin-top: 0;
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
                <pre>{summary_flag_table_content}</pre> <!-- Embed the content from the .txt file -->
            </div>
    """
    
    # Add the final CIGAR table if it exists
    if final_cigar_table_content:
        html_content += f"""
            <!-- Final CIGAR Table Section -->
            <div class="section">
                <div class="section-title">Final CIGAR Table</div>
                <pre>{final_cigar_table_content}</pre> <!-- Embed the CIGAR table content -->
            </div>
        """

    # Plots Section
    html_content += """
            <!-- Plots Section -->
            <div class="section">
                <div class="section-title">Plots</div>
                <div class="plot-container">
    """
    
    # Add each plot image from plot_paths
    for plot_name, plot_path in plot_paths.items():
        if plot_name != 'summary_flag_table':  # Skip the summary flag table file
            html_content += f'<img src="{plot_path}" alt="{plot_name} Plot">'

    # Finalize HTML
    html_content += """
                </div>
            </div>
        </div>
    </body>
    </html>
    """

    # Write the HTML content to a file
    try:
        with open(html_output_path, 'w') as f:
            f.write(html_content)
        print(f"HTML report saved to {html_output_path}")
    except Exception as e:
        print(f"Error writing the HTML report: {e}")

        
        

############################# CIGAR    
def toStringOutput(line):
    """Converts a SAM line to FASTA format. Simplified for this example."""
    cols = line.split("\t")
    return f">{cols[0]}\n{cols[9]}\n"  # Example: FASTA header and sequence.




def flagBinary(flag):
    """Converts a flag into a binary string of length 12."""
    flagB = bin(int(flag))[2:]  # Convert to binary and remove '0b'
    flagB = flagB.zfill(12)     # Pad with zeros to ensure 12 bits
    return flagB




def unmapped(sam_lines):
    """Analyze unmapped reads."""
    unmapped_count = 0
    with open("only_unmapped.fasta", "w") as unmapped_fasta, open("summary_unmapped.txt", "w") as summary_file:
        for line in sam_lines:
            col_line = line.split("\t")
            if len(col_line) < 11:
                continue  # Skip invalid lines
            flag = flagBinary(col_line[1])
            if flag[-3] == '1':  # Check if the read is unmapped
                unmapped_count += 1
                unmapped_fasta.write(toStringOutput(line))
        summary_file.write(f"Total unmapped reads: {unmapped_count}\n")
    return unmapped_count




def partiallyMapped(sam_lines):
    """Analyze partially mapped reads."""
    partially_mapped_count = 0
    with open("only_partially_mapped.fasta", "w") as partially_mapped_fasta, open("summary_partially_mapped.txt", "w") as summary_file:
        for line in sam_lines:
            col_line = line.split("\t")
            if len(col_line) < 11:
                continue  # Skip invalid lines
            flag = flagBinary(col_line[1])
            if flag[-2] == '1' and col_line[5] != "100M":  # Check if partially mapped
                partially_mapped_count += 1
                partially_mapped_fasta.write(toStringOutput(line))
        summary_file.write(f"Total partially mapped reads: {partially_mapped_count}\n")
    return partially_mapped_count




def processSAMFileAndCigar(filePath):
    """Parse CIGAR strings from a SAM file and output the parsed results to outputTable_cigar.txt."""
    # Open the output file to write the CIGAR operations counts
    with open("outputTable_cigar.txt", "w") as output_file:
        # Write the header for the output file
        output_file.write("M;I;D;S;H;N;P;X;Egal\n")

        try:
            with open(filePath, "r") as f:
                for line in f:
                    if line.startswith('@'):  # Skip header lines in SAM file
                        continue
                    columns = line.strip().split("\t")
                    if len(columns) < 11:  # Skip invalid lines (lines without a full SAM format)
                        continue

                    # Extract the CIGAR string from the 6th column (index 5)
                    cigar = columns[5]

                    if cigar != "*":  # Skip if the CIGAR string is "*"
                        # Parse the CIGAR string and get the counts of different operations
                        dico = readCigarOperations(cigar)

                        # Convert the parsed result into a semicolon-separated string and write to file
                        output_line = ";".join(str(dico.get(mut, 0)) for mut in ['M', 'I', 'D', 'S', 'H', 'N', 'P', 'X', '=']) + "\n"
                        output_file.write(output_line)

        except FileNotFoundError:
            print(f"Error: The file {filePath} was not found.")
        except Exception as e:
            print(f"Error reading the SAM file: {e}")

def readCigarOperations(cigar):
    """Parse the CIGAR string and return a dictionary with counts for each operation."""
    # Match CIGAR operations like '10M', '5I', '2D', etc.
    elements = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    dico = {mut: 0 for mut in ['M', 'I', 'D', 'S', 'H', 'N', 'P', 'X', '=']}  # Initialize dictionary

    for count, op in elements:
        dico[op] += int(count)

    return dico




def percentMutation(dico):
    """Calculate the percentage of each mutation type."""
    total = sum(dico.values())
    if total == 0:  # Avoid division by zero
        return "0.00;" * 9
    mutList = ['M', 'I', 'D', 'S', 'H', 'N', 'P', 'X', '=']
    return ";".join(f"{(dico.get(mut, 0) * 100 / total):.2f}" for mut in mutList) + ";"





def globalPercentCigar():
    """Global representation of CIGAR distribution."""
    try:
        # Open files with a context manager
        with open("outputTable_cigar.txt", "r") as outputTable, open("Final_Cigar_table.txt", "w") as FinalCigar:
            # Initialize counters for each operation
            nbReads, M, I, D, S, H, N, P, X, Egal = [0] * 10

            # Iterate over each line in the input file
            for line in outputTable:
                values = line.strip().split(";")
                
                # Skip header or invalid lines
                if values == ['M', 'I', 'D', 'S', 'H', 'N', 'P', 'X', 'Egal']:  # Assuming this is your header
                    continue

                if len(values) < 9:  # Skip invalid lines (those without 9 fields)
                    continue

                nbReads += 1

                # Process each value, catching conversion errors for individual fields
                try:
                    M += float(values[0]) if values[0] else 0
                    I += float(values[1]) if values[1] else 0
                    D += float(values[2]) if values[2] else 0
                    S += float(values[3]) if values[3] else 0
                    H += float(values[4]) if values[4] else 0
                    N += float(values[5]) if values[5] else 0
                    P += float(values[6]) if values[6] else 0
                    X += float(values[7]) if values[7] else 0
                    Egal += float(values[8]) if values[8] else 0
                except ValueError as e:
                    print(f"Warning: Invalid data in line '{line.strip()}': {e}")
                    continue  # Skip this line if there's an issue with the data

            # Check if there are any reads to process to avoid division by zero
            if nbReads == 0:
                FinalCigar.write("No reads to calculate global CIGAR distribution.\n")
                return
            
            # Calculate percentages for each operation
            M_percent = (M / nbReads) * 100
            I_percent = (I / nbReads) * 100
            D_percent = (D / nbReads) * 100
            S_percent = (S / nbReads) * 100
            H_percent = (H / nbReads) * 100
            N_percent = (N / nbReads) * 100
            P_percent = (P / nbReads) * 100
            X_percent = (X / nbReads) * 100
            Egal_percent = (Egal / nbReads) * 100



            # Calculate and write the global CIGAR statistics
            FinalCigar.write(f"Global cigar mutation observed:\n"
                             f"Alignment Match: {M/nbReads:.2f}\n"
                             f"Insertion: {I/nbReads:.2f}\n"
                             f"Deletion: {D/nbReads:.2f}\n"
                             f"Skipped Region: {S/nbReads:.2f}\n"
                             f"Soft Clipping: {H/nbReads:.2f}\n"
                             f"Hard Clipping: {N/nbReads:.2f}\n"
                             f"Padding: {P/nbReads:.2f}\n"
                             f"Sequence Match: {Egal/nbReads:.2f}\n"
                             f"Sequence Mismatch: {X/nbReads:.2f}\n")


            # Print the statistics to the terminal in a tabular format
            print("\n------------------------ CIGAR Distribution Stats ------------------------\n")
            print(f"{'Description':<40} {'Percentage (%)':<25} {'Total':<18}")
            print(f"{'--' * 35}\n")
                
            print(f"{'Alignment Match':<40} {M_percent:<20.2f} {M:>18.0f}")
            print(f"{'Insertion':<40} {I_percent:<20.2f} {I:>18.0f}")
            print(f"{'Deletion':<40} {D_percent:<20.2f} {D:>18.0f}")
            print(f"{'Skipped Region':<40} {S_percent:<20.2f} {S:>18.0f}")
            print(f"{'Soft Clipping':<40} {H_percent:<20.2f} {H:>18.0f}")
            print(f"{'Hard Clipping':<40} {N_percent:<20.2f} {N:>18.0f}")
            print(f"{'Padding':<40} {P_percent:<20.2f} {P:>18.0f}")
            print(f"{'Sequence Match':<40} {Egal_percent:<20.2f} {Egal:>18.0f}")
            print(f"{'Sequence Mismatch':<40} {X_percent:<20.2f} {X:>18.0f}")
                
            print(f"\n{'--' * 35}\n")
            print(f"{'Total Reads:':<40} {nbReads:>26}")
            print(f"{'--' * 35}\n\n\n")


            print("Final Cigar file created: Final_Cigar_table.txt")
            return "Final_Cigar_table.txt"




    except FileNotFoundError:
        print("Error: outputTable_cigar.txt not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


    ###################Mapping 
    
def parseSam(filePath):
    """
    Extract sequences from a SAM file.

    :param filePath: Path to the SAM file.
    :return: List of sequences (each as a tuple with read ID and sequence).
    """
    print("Entered parseSam") 
    sequences = []
    with open(filePath, 'r') as file:  # Open the SAM file for reading
        for line in file:  # Iterate over each line in the file
            if line.startswith('@'):  # Skip header lines (lines starting with '@')
                continue
            parts = line.strip().split('\t')  # Split the line by tabs to get the columns
            if len(parts) > 9:  # Ensure there are at least 10 columns in the line
                readId = parts[0]   # Column 1: Read ID
                sequence = parts[9]  # Column 10: Sequence
                sequences.append((readId, sequence))  # Store the read ID and sequence as a tuple
    print(f"Loaded {len(sequences)} sequences from the file.")  # Print the number of sequences loaded
    return sequences  # Return the list of sequences
    


def parseFastq(filePath):
    """
    Extract sequences from a FASTQ file.

    :param filePath: Path to the FASTQ file.
    :return: List of sequences (each as a tuple with read ID and sequence).
    """
    print("Entered parseFastq")
    sequences = []
    with open(filePath, 'r') as file:  # Open the FASTQ file for reading
        while True:
            try:
                # Read four lines at a time
                read_id = file.readline().strip()
                if not read_id:  # End of file
                    break
                sequence = file.readline().strip()
                file.readline()  # Skip the '+' line
                file.readline()  # Skip the quality line

                if read_id.startswith('@'):  # Validate FASTQ identifier
                    read_id = read_id[1:]  # Remove '@' for consistency
                    sequences.append((read_id, sequence))
            except Exception as e:
                print(f"Error processing FASTQ file: {e}")
                break

    print(f"Loaded {len(sequences)} sequences from the file.")  # Print the number of sequences loaded
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
  
    m, n = len(sequences), len(reference)  # Get the length of sequences and the reference
    scoreMatrix = np.zeros((m + 1, n + 1), dtype=int)  # Initialize the score matrix
    maxScore = 0  # Initialize the maximum score
    maxPos = None  # Initialize the position of the maximum score

    # Fill the score matrix with adjusted penalties based on the Smith-Waterman algorithm
    for i in range(1, m + 1):  # Iterate over the rows of the matrix (query sequences)
        for j in range(1, n + 1):  # Iterate over the columns of the matrix (reference sequence)
            # Calculate the match score, delete penalty, and insert penalty
            match = scoreMatrix[i - 1][j - 1] + (matchScore if sequences[i - 1] == reference[j - 1] else mismatchPenalty)
            delete = scoreMatrix[i - 1][j] + gapPenalty
            insert = scoreMatrix[i][j - 1] + gapPenalty
            scoreMatrix[i][j] = max(0, match, delete, insert)  # Take the maximum score (with zero threshold)

            if scoreMatrix[i][j] > maxScore:  # Update the maximum score and position
                maxScore = scoreMatrix[i][j]
                maxPos = (i, j)

    # Traceback to build the alignment
    alignedSeqOne = []  # Aligned sequence 1 (query)
    alignedSeqTwo = []  # Aligned sequence 2 (reference)
    i, j = maxPos  # Start the traceback from the maximum score position

    while scoreMatrix[i][j] != 0:  # Continue tracing back until we reach a score of 0
        if scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + (matchScore if sequences[i - 1] == reference[j - 1] else mismatchPenalty):
            alignedSeqOne.append(sequences[i - 1])  # Append the aligned query sequence
            alignedSeqTwo.append(reference[j - 1])  # Append the aligned reference sequence
            i -= 1  # Move to the previous row
            j -= 1  # Move to the previous column
        elif scoreMatrix[i][j] == scoreMatrix[i - 1][j] + gapPenalty:
            alignedSeqOne.append(sequences[i - 1])  # Append the query sequence with a gap
            alignedSeqTwo.append('-')  # Append a gap for the reference sequence
            i -= 1  # Move to the previous row
        else:
            alignedSeqOne.append('-')  # Append a gap for the query sequence
            alignedSeqTwo.append(reference[j - 1])  # Append the reference sequence with a gap
            j -= 1  # Move to the previous column

    return ''.join(reversed(alignedSeqOne)), ''.join(reversed(alignedSeqTwo)), maxScore  # Return the aligned sequences and the score




def compareSequences(reference, sequences, smithWatermanFunc, scoreThreshold=50):
    """
    Compare query sequences with reference sequences using Smith-Waterman algorithm.

    :param reference: List of sequences from the reference SAM file.
    :param sequences: List of sequences from the query SAM file.
    :param smithWatermanFunc: Function that implements Smith-Waterman algorithm.
    :param scoreThreshold: Minimum score threshold for considering an alignment valid.
    :return: List of matches with alignment details.
    """
    results = []  # Initialize an empty list to store the results

    for queryId, querySeq in sequences:  # Iterate over the query sequences
        for refId, refSeq in reference:  # Iterate over the reference sequences
            # Perform the alignment using Smith-Waterman
            alignedRef, alignedQuery, score = smithWatermanFunc(refSeq, querySeq)
            
            # Add the result if the score is higher than the threshold
            if score >= scoreThreshold:  # Check if the alignment score meets the threshold
                results.append({
                    "queryId": queryId,  # Store the query ID
                    "refId": refId,  # Store the reference ID
                    "score": score,  # Store the alignment score
                    "alignedRef": alignedRef,  # Store the aligned reference sequence
                    "alignedQuery": alignedQuery  # Store the aligned query sequence
                })
    
    return results  # Return the list of alignment results




def alignSequences(args, smithWaterman):
    """
    Process sequences by comparing reference and query sequences, and print the alignment results.

    :param args: Parsed command-line arguments.
    """
    # Check if the necessary arguments (reference, input, outputFile) are provided
    if not args.reference or not args.input or not args.outputFile:
        print("Error: Reference, query, and output file paths are required.")
        return  # If not, print an error and exit the function
    
    # Load sequences from the reference and input (query) SAM files, ignoring the headers
    reference = parseSam(args.reference)  # Load sequences from the reference file
    sequences = parseFastq(args.input)  # Load sequences from the query file

    # Check if the reference sequences were loaded successfully
    if not reference:
        print(f"Error: Failed to load sequences from reference file '{args.reference}'.")
        return  # If the reference sequences could not be loaded, print an error and exit

    # Check if the query sequences were loaded successfully
    if not sequences:
        print(f"Error: Failed to load sequences from query file '{args.input}'.")
        return  # If the query sequences could not be loaded, print an error and exit

    # Perform the alignment using the Smith-Waterman algorithm and save the results to the output file
    with open(args.outputFile, 'w') as outputFile:  # Open the output file for writing
        # Compare the reference and query sequences using the Smith-Waterman algorithm
        matches = compareSequences(reference, sequences, smithWaterman)
        # Iterate over the matches and write the alignment details to the output file
        for match in matches:
            # Write the alignment details for each match
            outputFile.write(f"Query {match['queryId']} aligned with Reference {match['refId']}\n")
            outputFile.write(f"Score: {match['score']}\n")
            outputFile.write(f"Aligned Reference: {match['alignedRef']}\n")
            outputFile.write(f"Aligned Query: {match['alignedQuery']}\n\n")

    # Print a message indicating that the alignment results have been saved
    print(f"Alignment results saved to {args.outputFile}")
