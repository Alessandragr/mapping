#!/usr/bin/python3
#-*- coding : utf-8 -*-

############### IMPORT MODULES ############### 

import argparse
from collections import defaultdict

from flags import flags

from utils import fileVerifier, countReads, readPerChrom, readPerMAPQ, countReadsByFlags, saveResults, executePlots, filterSam, mappedRead, mappedPrimaryReads, globalPercentCigar, processSAMFileAndCigar




# OPTION LIST:
# -h or --help: Displays help information.
# -i or --input: Specifies the path to the input SAM file (.sam).
# -o or --outputFile: Specifies the path to the output file.
# -oD or --outputDir: Specifies the directory for output files.
# -cR or --countReads: Counts total reads, mapped reads, unmapped reads, duplicated reads, with a filtering option on the mapping quality.
# -rC or --readPerChrom: Counts the number of reads per chromosome.
# -rMQ or --readPerMAPQ: Counts the number of reads for each MAPQ score.
# -cRF or --countReadsByFlags: Counts the number of reads for each FLAG value.
# -sR or --saveResults: Saves the results and associated plots to an HTML file.
# -fS or --filterSam: Filters SAM file by mapping quality and creates a new filtered SAM file.
# -mR or --mappedRead: Filters SAM file to keep only mapped reads.
# -mPR or --mappedPrimaryReads: Filters SAM file to keep only primarily mapped reads.
# -m or --minQ: Sets the minimum MAPQ score for filtering reads.
# -eP or --executePlots: Executes Flag Stats analysis and generates plots.
# -gPC or --globalPercentCigar: Calculate and display global CIGAR mutation percentages.



# SYNOPSIS:
# SamReader.py -h or --help                                      # Displays help.
# SamReader.py -i <file>                                         # Check if the SAM file is valid.
# SamReader.py -i <file> -cR                                     # Prints read statistics (total, mapped, unmapped, duplicated).
# SamReader.py -i <file> -cR -m <mappingQuality>                 # Prints filtered read statistics (total, mapped, unmapped, duplicated).
# SamReader.py -i <file> -rC                                     # Prints the number of reads per chromosome.
# SamReader.py -i <file> -rMQ                                    # Prints the number of reads for each MAPQ score.
# SamReader.py -i <file> -cRF                                    # Prints the number of reads for each FLAG value.
# SamReader.py -i <file> -fS -o <outputFile> -m <mappingQuality> # Filters SAM file by mapping quality and creates a filtered SAM file.
# SamReader.py -i <file> -mR -o <outputFile>                     # Filters SAM file to keep only mapped reads.
# SamReader.py -i <file> -mPR -o <outputFile>                    # Filters SAM file to keep only primarily mapped reads.
# SamReader.py -i <file> -eP                                     # Executes Flag Stats analysis and generates plots.
# SamReader.py -i <file> -cRF -eP                                # Combines flag-based read counting with Flag Stats plot generation.
# SamReader.py -i <file> -eP -sR -oD <outputDir>                                # Combines plots creation with saving results to HTML.
# SamReader.py -i <file> -gPC                                    # Calculates and displays global CIGAR mutation percentages.
# SamReader.py -i <file> -gPC -eP -sR  -oD <outputDir>           # Combines plots creation with saving results to HTML and cigar analysis.



def main():
    """
    Main function to handle various SAM file analyses and sequence alignment tasks.
    """
    parser = argparse.ArgumentParser(description="Analyze a SAM file and provide various statistics.")

    # Arguments for file input and output
    parser.add_argument('-i', '--input', required=True, help="Path to the SAM file.")
    parser.add_argument('-o', '--outputFile', required=False, help="Path to the output SAM file.")
    parser.add_argument('-oD', '--outputDir', required=False, help="Path to the output SAM file.")

    # Arguments for read counting and analysis
    parser.add_argument('-cR', '--countReads', action='store_true', help="Count reads (total, mapped, etc.).")
    parser.add_argument('-rC', '--readPerChrom', action='store_true', help="Count reads per chromosome.")
    parser.add_argument('-rMQ', '--readPerMAPQ', action='store_true', help="Count reads based on MAPQ scores.")
    parser.add_argument('-cRF', '--countReadsByFlags', action='store_true', help="Count reads based on FLAG values.")

    

    # Arguments for filtering and saving new SAM files
    parser.add_argument('-fS', '--filterSam', action='store_true', help="Filter SAM file by mapping quality.")
    parser.add_argument('-mR', '--mappedRead', action='store_true', help="Filter SAM file by keeping mapped reads.")
    parser.add_argument('-mPR', '--mappedPrimaryReads', action='store_true', help="Filter SAM file by keeping primarly mapped reads.")
    

    # Argument for MAPQ score filtering
    parser.add_argument('-m', '--minQ', type=int, default=0, help="Minimum MAPQ score for filtering reads.")

    # Argument for executing plots (calls executePlots)
    parser.add_argument('-eP', '--executePlots', action='store_true', help="Execute Flag Stats analysis.")
    
    # Arguments for saving results and generating output
    parser.add_argument('-sR', '--saveResults', action='store_true', help="Save results and graphs to an HTML file.")
    
   
    # Arguments for CIGAR ANALYSIS
    parser.add_argument('-gPC','--globalPercentCigar', action='store_true', help="Calculate and display global CIGAR mutation percentages.")

    args = parser.parse_args()

    # Display help if no arguments are provided
    if not any(vars(args).values()):
        parser.print_help()
        return

    # Verify the input file (this function should be defined elsewhere in your code)
    if not fileVerifier(args.input):
        print(f"Error: The input file {args.input} does not exist or cannot be opened.")
        return

    plot_paths = []  # Initialize plot paths


    # Existing code for read counting, flag stats, and result saving
    if args.countReads:
        countReads(args.input, args.minQ)

    if args.readPerChrom:
        readPerChrom(args.input)

    if args.readPerMAPQ:
        readPerMAPQ(args.input)

    if args.countReadsByFlags:
        countReadsByFlags(args.input)
        
    if args.filterSam:
        if not args.outputFile:
            print("Error: Output file path (-o/--outputFile) is required for filtering.")
            return
        filterSam(args.input, args.outputFile, args.minQ)

    if args.mappedRead:
        if not args.outputFile:
            print("Error: Output file path (-o/--outputFile) is required for mapped read filtering.")
            return
        mappedRead(args.input, args.outputFile)
        
    if args.mappedPrimaryReads:
        mappedPrimaryReads(args.input, args.outputFile)


    if args.globalPercentCigar:
        processSAMFileAndCigar(args.input)
        globalPercentCigar()
    
    if args.executePlots:
        # Execute the function and retrieve results
        plot_paths = executePlots(
            filePath=args.input,
            outputDir="result",
            outputFile="analysis",
            minQ=args.minQ
        )


# Ensure the plot_paths is populated and then save results
    if args.saveResults:
        if not plot_paths:
            print("Error: Plot paths are empty. Ensure executePlots is run before saving results.")
            return

        # Safely get the summary flag table path
        summary_flag_table_path = plot_paths.get('summary_flag_table')

        # Get the path to the final CIGAR table
        final_cigar_table_path = globalPercentCigar()

        # Proceed only if the summary_flag_table_path exists
        if summary_flag_table_path:
            saveResults(
                plot_paths=plot_paths,  # Pass the dictionary containing all paths
                summary_flag_table_path=summary_flag_table_path,
            )
        else:
            print("Error: Summary flag table path is not available. Ensure executePlots ran correctly.")



if __name__ == "__main__":
    main()





