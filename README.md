
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
- Python 3.12.7 or higher
- The following Python libraries:
  - `matplotlib`
  - `pandas`
  - `numpy`
  - `os`

### Install Required Libraries
Install the necessary libraries using `pip`:

pip install matplotlib pandas numpy os


### Usage

## Usage

### Command-Line Options
| Option                  | Description                                                                                 |
|-------------------------|---------------------------------------------------------------------------------------------|
| `-h` or `--help`        | Displays help information.                                                                  |
| `-i` or `--input`       | Specifies the path to the input SAM file (.sam).                                            |
| `-o` or `--output`      | Specifies the path to the output file.                                                      |
| `-cR` or `--countReads` | Counts total, mapped, unmapped, and duplicated reads, with an option for MAPQ filtering.    |
| `-rC` or `--readPerChrom` | Counts the number of reads per chromosome.                                                |
| `-rMQ` or `--readPerMAPQ` | Counts the number of reads for each MAPQ score.                                            |
| `-cRF` or `--countReadsByFlags` | Counts the number of reads for each FLAG value.                                      |
| `-sR` or `--save-results` | Saves the results and associated plots to an HTML file.                                   |
| `-fS` or `--filterSam`  | Filters the SAM file by mapping quality and creates a new filtered SAM file.                |
| `-mr` or `--mappedRead` | Filters the SAM file by mapped reads and creates a filtered SAM file.                       |
| `-eFS` or `--executeFlagStats` | Executes FLAG analysis and generates plots.                                           |

### Command Synopsis
Here are examples of how to use the tool:

## SYNOPSIS:
 | Command                                                       | Description                                                                                     |
|---------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| `SamReader.py -h` or `--help`                                  | Displays help.                                                                                 |
| `SamReader.py -i <file>`                                       | Check if SAM file is valid.                                                                     |
| `SamReader.py -i <file> -cR`                                   | Prints read statistics (total, mapped, unmapped, duplicated).                                  |
| `SamReader.py -i <file> -cR -m <mappingquality>`               | Prints read statistics (total, mapped, unmapped, duplicated, filtered).                        |
| `SamReader.py -i <file> -rC`                                   | Prints the number of reads per chromosome.                                                      |
| `SamReader.py -i <file> -rMQ`                                  | Prints the number of reads per MAPQ score.                                                      |
| `SamReader.py -i <file> -cRF`                                  | Prints the number of reads per FLAG value.                                                      |
| `SamReader.py -i <file> -sR`                                   | Saves results (including plots) to an HTML file.                                                 |
| `SamReader.py -i <file> -cR -sR`                               | Combines read statistics and saves results to HTML.                                              |
| `SamReader.py -i <file> -fS -o <outputfile> -m <mappingQuality>`| Filter SAM file by mapping quality and create a filtered SAM file.                              |
| `SamReader.py -i <file> -mr -o <outputfile>`                   | Filter SAM file by mapped read and create a filtered SAM file.                                  |
| `SamReader.py -i <file> --eFS`                                 | Execute FLAG Stats analysis and create plots.                                                   |


#### Help and Validation
```bash
python3 SamReader.py -h
python3 SamReader.py --help
python3 SamReader.py -i path/to/input.sam
```
  

  



### Installing






## Authors

Contributors names and contact info

ex. Najat AMOUKOU  
ex. Alessandra 

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments


