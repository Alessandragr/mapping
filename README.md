
# SAM File Analysis Tool

## Description
This is a Python-based command-line tool for analyzing SAM (Sequence Alignment/Map) files. The script processes SAM files to extract essential information, perform analysis, and visualize the data. It is useful for bioinformatics workflows involving read alignment analysis.

## Features
- **File Verification**: Ensures the input SAM file exists, is valid, and is non-empty.
- **Read Counting**: Computes statistics including:
  - Total reads
  - Mapped reads
  - Unmapped reads
  - Duplicated reads
  - Filtered reads based on MAPQ.
- **FLAG Interpretation**: Decodes SAM file FLAGS to explain the alignment properties.
- **Visualization**: Generates a html file showing the statistics of the sam file.
- **Filtering SAm file**:

## Requirements
To run this tool, you need:
- Python 3.6 or higher
- The following Python libraries:
  - `matplotlib`
  - `pandas`
  - `numpy`
  - `os`

### Install Required Libraries
Install the necessary libraries using `pip`:

pip install matplotlib pandas numpy os


### Usage



## OPTION LIST:
# -h or --help: Displays help information.
# -i or --input: Specifies the path to the input SAM file (.sam).
# -o or --output: Specifies the path to the output file 
# -cR or --countReads: Counts total reads, mapped reads, unmapped reads,duplicated reads, with a filtering option on the mapping quality.
# -rC or --readPerChrom: Counts the number of reads per chromosome.
# -rMQ or --readPerMAPQ: Counts the number of reads for each MAPQ score.
# -cRF or --countReadsByFlags: Counts the number of reads for each FLAG value.
# -sR or --save-results: Saves the results and associated plots to an HTML file.
# -fS or  --filterSam: Filter SAM file by mapping quality and create a new filterd sam file.
# -mr or  --mappedRead:  # Filter SAM file by mapped reaad  and create a filterd sam file
# -eFS  or  --executeFlagStats'         # Execute Flag Stats analysis and create plots
  

## SYNOPSIS:
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
# SamReader.py -i <file> -mr-o <outpufile> -  # Filter SAM file by mapped reaad  and create a filterd sam file.
# SamReader.py -i <file> --eFS            # Execute Flag Stats analysis and create plots
  








### Installing



### Executing program


## Help





## Authors

Contributors names and contact info

ex. Najat AMOUKOU  
ex. Alessandra Gonçalves

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments


