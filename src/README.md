
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
    - Filterring by mapping quality
    - FIltering by only keeping mapped reads
    - Filtering by keeping primarly mapped reads
  

## Requirements
To run this tool, you need:
- Python 3.12.7 or higher
- The following Python libraries:
  - `matplotlib`
  - `pandas`
  - `numpy`
  - `os`
  - `argparse`

### Install Required Libraries
Install the necessary libraries using `pip`:

```bash
pip install matplotlib pandas numpy argparser
```

### Usage



### Command-Line Options
| Option                  | Description                                                                                 |
|-------------------------|---------------------------------------------------------------------------------------------|
| `-h` or `--help`        | Displays help information.                                                                  |
| `-i` or `--input`       | Specifies the path to the input SAM file (.sam).                                            |
| `-o` or `--output`      | Specifies the path to the output file.                                                      |
| `-m` or `--minQ`       | Minimum MAPQ score for filtering reads.                                          |
| `-r` or `--reference`       | Reference sequence for alignment.                                           |
| `-cR` or `--countReads` | Counts total, mapped, unmapped, and duplicated reads, with an option for MAPQ filtering.    |
| `-rC` or `--readPerChrom` | Counts the number of reads per chromosome.                                                |
| `-rMQ` or `--readPerMAPQ` | Counts the number of reads for each MAPQ score.                                            |
| `-cRF` or `--countReadsByFlags` | Counts the number of reads for each FLAG value.                                      |                             |
| `-fS` or `--filterSam`  | Filters the SAM file by mapping quality and creates a new filtered SAM file.                |
| `-mr` or `--mappedRead` | Filters the SAM file by mapped reads and creates a filtered SAM file.                       |
| `-mPR` or `--mappedPrimarlyReads` | Filters the SAM file by primarly  mapped reads and creates a filtered SAM file.                       |
| `-eP` or `--executePlots` | Executes FLAG analysis and generates plots.                                           |
| `-aS` or `--alignSequences ` | Aligns reference and query sequences using an alignment algorithm.                                            |
| `-sR`  or `--saveResults` |  Saves results (including plots) to an HTML file.                                        |

### Command Synopsis
Here are examples of how to use the tool:

## SYNOPSIS:
| Command                                                       | Description                                                                                     |
|---------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| `SamAnalyser.py -h` or `--help`                               | Displays help.                                                                                 |
| `SamAnalyser.py -i <file>`                                    | Check if the SAM file is valid.                                                                |
| `SamAnalyser.py -i <file> -cR`                                | Prints read statistics (total, mapped, unmapped, duplicated).                                  |
| `SamAnalyser.py -i <file> -cR -m <mappingQuality>`            | Prints filtered read statistics (total, mapped, unmapped, duplicated).                         |
| `SamAnalyser.py -i <file> -rC`                                | Prints the number of reads per chromosome.                                                     |
| `SamAnalyser.py -i <file> -rMQ`                               | Prints the number of reads for each MAPQ score.                                                |
| `SamAnalyser.py -i <file> -cRF`                               | Prints the number of reads for each FLAG value.                                                |
| `SamAnalyser.py -i <file> -fS -o <outputFile> -m <mappingQuality>` | Filters SAM file by mapping quality and creates a filtered SAM file.                           |
| `SamAnalyser.py -i <file> -mR -o <outputFile>`                | Filters SAM file to keep only mapped reads.                                                    |
| `SamAnalyser.py -i <file> -mPR -o <outputFile>`               | Filters SAM file to keep only primarily mapped reads.                                          |
| `SamAnalyser.py -i <file> -eP`                                | Executes Flag Stats analysis and generates plots.                                              |                            |
| `SamAnalyser.py -i <file> -eP -sR`                            | Combines read statistics plots with saving results to HTML.                             |
| `SamAnalyser.py -i <file> -r <reference> -o <outputFile> -aS` | Aligns reference and query sequences using an alignment algorithm.                             |


#### Help and Validation
```bash
python3 SamAnalyser.py -h
python3 SamAnalyser.py --help
python3 SamAnalyser.py -i path/to/input.sam
```
  

  



### Installing
To use this software:

1. Clone the repository from GitHub to your computer.
2. Install the required dependencies listed above.
3. Run the script as per the provided commands.







## Authors

Contributors names and contact info

ex. Najat AMOUKOU  
ex. Alessandra RIBEIRO

## Version History

* 0.1
    * Initial Release

## License

Université de Montpellier

Copyright (c) [2024] [GONCALVES RIBEIRO Alessandra, IBRAHIM AMOUKOU Najat]

Authors contact: [alessandra.goncalves-ribeiro@etu.umontpellier.fr, najat.ibrahim-amoukou@etu.umontpellier.fr]

Version: 0.0.1

Date: 17 December 2024

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 0.0.1 of the License, or (at your option) any later version. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.


## Acknowledgments


