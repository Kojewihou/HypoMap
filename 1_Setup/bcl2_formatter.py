"""
This script processes FASTQ files in a given directory. It filters out files matching a specific pattern,
groups them by sample name, and then processes each sample to identify and rename key files based on sequence 
characteristics. The primary goal is to prepare the raw data for integration model setup by renaming the files 
appropriately and identifying redundant sequences.

Main functionalities include:
1. Filtering out BCL2-related FASTQ files.
2. Grouping and processing samples based on their filenames.
3. Analyzing sequences within each file to determine characteristics such as the longest sequence length and T-richness.
4. Renaming files based on the analysis, marking redundant files accordingly.

Modules:
- argparse: For command-line argument parsing.
- glob: To find all files matching a specific pattern.
- gzip: For reading compressed FASTQ files.
- multiprocessing: To parallelize the processing of multiple samples.
- re: Regular expressions for pattern matching.
- os: For file manipulation.

Global variables:
- DATADIR: The directory containing FASTQ files.
- BCL2_PATTERN: A regex pattern to filter out BCL2-related FASTQ files.
"""

import argparse
from glob import glob
import gzip
import multiprocessing
import re
import os

global DATADIR
BCL2_PATTERN = r'_S[0-9]*_L[0-9]*_[R,I][0-9]_[0-9]*\.fastq\.gz'

def main(args: argparse.Namespace) -> None:
    """
    Main function that sets up the directory for processing and starts the multiprocessing pool to process samples.

    Args:
        args (argparse.Namespace): Parsed command-line arguments containing the directory to be processed.
    """
    global DATADIR
    DATADIR = args.dir
    
    files_glob = glob(f'{DATADIR}/*.fastq.gz')
    filtered_files_glob = [f for f in files_glob if not re.search(BCL2_PATTERN, f)]
    
    sample_glob = set(get_sample_name(f) for f in filtered_files_glob)
    
    with multiprocessing.Pool() as pool:
        pool.map(process_sample, sample_glob)

def get_sample_name(file_path: str) -> str:
    """
    Extracts the sample name from a given file path.

    Args:
        file_path (str): Path to the FASTQ file.

    Returns:
        str: The extracted sample name.
    """
    base = os.path.basename(file_path)
    filename_without_extension = base.split('.')[0]
    
    if filename_without_extension.count('_') == 2:
        sample = filename_without_extension.split('_')[:-1]
        return "_".join(sample)
    else:
        sample = filename_without_extension
        return sample
            
def process_sample(sample: str) -> None:
    """
    Processes a sample by analyzing its associated FASTQ files and renaming them based on the analysis.

    Args:
        sample (str): The sample name to process.
    """
    global DATADIR
    files = glob(os.path.join(DATADIR, f'{sample}*.fastq.gz'))

    filelist = [file_analysis(f) for f in files]
    
    sorted_filelist = sorted(filelist, key=lambda f: (f.line_length, f.t_score), reverse=True)
    
    for i, seq_file in enumerate(sorted_filelist):
        if i == 0:
            os.rename(seq_file.filename, os.path.join(DATADIR, f'{sample}_S1_L001_R2_001.fastq.gz'))
        elif i == 1:
            os.rename(seq_file.filename, os.path.join(DATADIR, f'{sample}_S1_L001_R1_001.fastq.gz'))
        else:
            os.rename(seq_file.filename, seq_file.filename + '.redundant')

class file_analysis:
    """
    A class for analyzing sequences in a FASTQ file to determine the longest sequence length and T-richness (t_score).
    
    Attributes:
        filename (str): The name of the FASTQ file being analyzed.
        line_length (int): The length of the longest sequence in the file.
        t_score (float): A score representing the T-richness of sequences in the file.
    """
    def __init__(self, filename: str, max_lines: int = 36) -> None:
        """
        Initializes the file_analysis class by reading sequences from the given FASTQ file and calculating
        the line length and t_score.

        Args:
            filename (str): The name of the FASTQ file to analyze.
            max_lines (int, optional): The maximum number of sequences to analyze. Defaults to 36.
        """
        self.filename = filename
        
        sequences = []
        with gzip.open(filename, 'rt') as fp:
            for i, line in enumerate(fp):
                if i % 4 == 1:
                    sequences.append(line.strip())
                if i > max_lines:
                    break
        
        self.line_length = max(len(seq) for seq in sequences)
        
        total_sequence = "".join(sequences)
        t_matches = re.findall(r'T+', total_sequence.upper())
        
        if t_matches:
            self.t_score = 1 / sum(len(match) ** 2 for match in t_matches)
        else:
            self.t_score = 1

    def __repr__(self) -> str:
        """
        Returns a string representation of the file_analysis object, including the filename, t_score, and line_length.

        Returns:
            str: A string representation of the file_analysis object.
        """
        return f"{self.filename}\n\tt_score = {self.t_score}\n\tline_length = {self.line_length}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Setup raw data for integration model setup')

    parser.add_argument('dir', metavar='[ DATA DIR ]')

    args = parser.parse_args()
    main(args)
