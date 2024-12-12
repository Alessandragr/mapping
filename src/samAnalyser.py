#!/usr/bin/python3
#-*- coding : utf-8 -*-

############### IMPORT MODULES ############### 

import os,re,sys,argparse
from collections import defaultdict

from flags import flags

from utils import fileVerifier, countReads, readPerChrom, readPerMAPQ, countReadsByFlags, parseSam, saveResults,executePlots, filterSam, mappedRead





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
    # Arguments for ParseSam
    parser.add_argument('-pS', '--parseSam', action='store_true', help="Parse sam file")
  
    args = parser.parse_args()

    # Display help if no arguments are provided
    if not any(vars(args).values()):
        parser.print_help()
        return

    # Verify the input file (this function should be defined elsewhere in your code)
    if not fileVerifier(args.input):
        print(f"Error: The input file {args.input} does not exist or cannot be opened.")
        return

    
    if args.parseSam:
        if not args.input:
            print("Error:sequence is required.")
            return
        parseSam(args.input)
        print('file parsed')
    
    

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




