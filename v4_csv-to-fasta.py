
# Import necessary Python modules
import os  # For file and path operations
import sys  # For system-specific parameters and functions
import math  # For mathematical operations
import re  # For regular expressions
import csv  # For reading and writing CSV files
import argparse  # For parsing command-line arguments
from typing import List, Dict, Tuple  # For type hinting
from datetime import datetime  # For generating timestamps

# Define a class to represent a binder molecule
class Binder:
    def __init__(self, csv_line: str, orig_line_index: int):
        # Initialize a Binder object from a CSV line
        # csv_line: a string containing comma-separated values
        # orig_line_index: the original index of this line in the CSV file
        
        # Split the CSV line into individual fields
        fields = csv_line.split(',')
        
        # Create a unique label for this binder
        # Combine the first field, the original line index, and the last field (stripped of whitespace)
        self.label = f"{fields[0]}_{orig_line_index}_{fields[6].strip()}"
        
        # Store the sequence of the binder
        self.sequence = fields[1]
        
        # Convert numerical fields to float and store them
        self.ipae = float(fields[2])  # IPAE score
        self.rmsd = float(fields[3])  # RMSD value
        self.plddt = float(fields[4])  # pLDDT score
        self.iptm = float(fields[5])  # IPTM score
        
        # Store the original line index for reference
        self.orig_line_index = orig_line_index

    def max_alanine_stretch_length(self) -> int:
        # Find the longest stretch of consecutive alanines in the sequence
        # Use regular expression to find all stretches of 'A' (Alanine)
        ala_stretches = re.findall('A+', self.sequence)
        
        # Return the length of the longest stretch
        # If no stretches found, return 0 (this is what the 'default=0' does)
        return max((len(stretch) for stretch in ala_stretches), default=0)

# Function to parse command-line arguments
def parse_arguments():
    # Create an ArgumentParser object to handle command-line arguments
    parser = argparse.ArgumentParser(description='Process CSV file to FASTA and ranked CSV')
    
    # Add arguments to the parser
    # Each add_argument call defines a new command-line option
    parser.add_argument('-i', '--input_csv', required=True, help='Path to input CSV file')
    parser.add_argument('-f', '--output_fasta', help='Path to output FASTA file')
    parser.add_argument('-o', '--output_csv', help='Path to output ranked CSV file')
    parser.add_argument('-n', '--write_top_n', type=int, help='Number of top sequences to write to FASTA')
    parser.add_argument('-p', '--min_plddt', type=float, default=0.8, help='Minimum pLDDT threshold (default: 0.8)')
    parser.add_argument('-s', '--similarity_cutoff', type=int, default=18, help='Sequence similarity cutoff (default: 18)')
    parser.add_argument('-w', '--ipae_window_size', type=float, default=1.0, help='IPAE window size (default: 1.0)')
    parser.add_argument('-a', '--max_ala_stretch', type=int, default=7, help='Maximum allowed alanine stretch (default: 7)')

    # Parse the arguments
    args = parser.parse_args()

    # Validate the similarity_cutoff argument
    if args.similarity_cutoff < 1:
        print('ERROR: similarity_cutoff must be >= 1')
        sys.exit(1)  # Exit the program if the validation fails

    # Generate default output filenames if not provided by the user
    if not args.output_fasta:
        # If no FASTA output file is specified, create one based on the input CSV filename
        args.output_fasta = os.path.splitext(args.input_csv)[0] + '_ranked.fasta'
    if not args.output_csv:
        # If no CSV output file is specified, create one based on the input CSV filename
        args.output_csv = os.path.splitext(args.input_csv)[0] + '_ranked.csv'

    # Check if output files already exist and append timestamp if they do
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if os.path.exists(args.output_fasta):
        base, ext = os.path.splitext(args.output_fasta)
        args.output_fasta = f"{base}_{timestamp}{ext}"
    if os.path.exists(args.output_csv):
        base, ext = os.path.splitext(args.output_csv)
        args.output_csv = f"{base}_{timestamp}{ext}"

    return args

# Function to read CSV file
def read_csv_file(file_path: str) -> List[str]:
    # Open the file in read mode
    with open(file_path, 'r') as csvf:
        # Read all lines from the file and return them as a list
        return csvf.readlines()

# Function to sort binder dictionary by a specific attribute
def sort_binder_dict(bdict: Dict[str, Binder], key: str) -> Dict[str, Binder]:
    # Sort the dictionary based on the specified attribute (key) of the Binder objects
    # The 'sorted' function returns a list of tuples, which we convert back to a dictionary
    return dict(sorted(bdict.items(), key=lambda item: getattr(item[1], key)))

# Function to subset binder dictionary based on IPAE range
def subset_binder_dict_by_ipae(bdict: Dict[str, Binder], min_ipae: float, max_ipae: float) -> Dict[str, Binder]:
    # Create a new dictionary containing only the binders whose IPAE falls within the specified range
    return {lbl: binder for lbl, binder in bdict.items() 
            if min_ipae <= binder.ipae < max_ipae}

# Function to resort binder dictionary by RMSD within IPAE windows
def resort_binder_dict_by_rmsd(bdict: Dict[str, Binder], ipae_window_size: float) -> Dict[str, Binder]:
    # Check if the input dictionary is empty
    if not bdict:
        return {}  # Return an empty dictionary if bdict is empty
    
    # Find the minimum and maximum IPAE values in the dictionary
    min_ipae = min(binder.ipae for binder in bdict.values())
    max_ipae = max(binder.ipae for binder in bdict.values())
    
    # Calculate the number of windows needed
    num_windows = math.ceil((max_ipae - min_ipae) / ipae_window_size)
    
    resorted_dict = {}
    # Iterate through each window
    for i in range(num_windows):
        # Calculate the IPAE range for this window
        window_min_ipae = min_ipae + (i * ipae_window_size)
        window_max_ipae = window_min_ipae + ipae_window_size
        
        # Get the subset of binders that fall within this IPAE window
        window_dict = subset_binder_dict_by_ipae(bdict, window_min_ipae, window_max_ipae)
        
        # Sort this subset by RMSD and add it to the resorted dictionary
        resorted_dict.update(sort_binder_dict(window_dict, 'rmsd'))
    
    return resorted_dict

# Function to calculate sequence difference
def seq_difference(seq1: str, seq2: str) -> int:
    # Calculate the number of positions where the two sequences differ
    # This uses a generator expression with the 'sum' function
    return sum(a != b for a, b in zip(seq1, seq2))

# Function to find minimum sequence difference against a set of sequences
def min_seq_difference_vs_set(seq_to_check: str, set_dict: Dict[str, Binder]) -> int:
    # If the set_dict is empty, return the length of the sequence to check
    # This ensures that if there are no sequences to compare against, the difference is maximized
    if not set_dict:
        return len(seq_to_check)
    
    # Find the minimum difference between seq_to_check and all sequences in set_dict
    return min(seq_difference(seq_to_check, binder.sequence) for binder in set_dict.values())

# Function to process binders and separate them into sets based on criteria
def process_binders(binders: List[Binder], min_plddt: float, similarity_cutoff: int, ipae_window_size: float) -> Tuple[Dict[str, Binder], Dict[str, Binder], Dict[str, Binder], int]:
    # Create initial dictionary of all binders
    binder_dict = {binder.label: binder for binder in binders}
    
    # Separate binders based on pLDDT threshold
    set1 = {lbl: binder for lbl, binder in binder_dict.items() if binder.plddt >= min_plddt}
    set2 = {lbl: binder for lbl, binder in binder_dict.items() if binder.plddt < min_plddt}
    
    # Sort sets by IPAE
    set1 = sort_binder_dict(set1, 'ipae')
    set2 = sort_binder_dict(set2, 'ipae')
    
    # Resort set1 by RMSD if ipae_window_size > 0
    if ipae_window_size > 0:
        set1 = resort_binder_dict_by_rmsd(set1, ipae_window_size)
    
    # Process set1 to remove similar sequences
    set1_checked = {}
    set3 = {}
    num_identical = 0
    
    # Iterate through set1 and categorize binders based on their similarity to already processed binders
    for lbl, binder in set1.items():
        min_difference = min_seq_difference_vs_set(binder.sequence, set1_checked)
        if min_difference >= similarity_cutoff:
            # If the binder is sufficiently different, add it to set1_checked
            set1_checked[lbl] = binder
        elif 0 < min_difference < similarity_cutoff:
            # If the binder is similar but not identical, add it to set3
            set3[lbl] = binder
        else:
            # If the binder is identical to an already processed binder, increment the counter
            num_identical += 1
    
    return set1_checked, set2, set3, num_identical

# Function to write binders to FASTA file
def write_fasta(file_path: str, binder_dict: Dict[str, Binder], write_top_n: int = None):
    # Open the file in write mode
    with open(file_path, 'w') as ffile:
        # Iterate through the binder dictionary
        for i, (lbl, binder) in enumerate(binder_dict.items()):
            # If write_top_n is specified and we've reached that number, stop writing
            if write_top_n is not None and i >= write_top_n:
                break
            # Write the FASTA header line
            ffile.write(f'>{binder.label}_{binder.ipae}_{binder.rmsd}_{binder.plddt}_{binder.iptm}\n')
            # Write the sequence line
            ffile.write(f'{binder.sequence}\n')

# Function to write ranked binders to CSV file
def write_ranked_csv(file_path: str, binder_dict: Dict[str, Binder]):
    # Open the file in write mode, with newline='' to handle line endings correctly
    with open(file_path, 'w', newline='') as csvfile:
        # Create a CSV writer object
        writer = csv.writer(csvfile)
        # Write the header row
        writer.writerow(['Rank', 'Label', 'Sequence', 'IPAE', 'RMSD', 'pLDDT', 'IPTM'])
        # Write data rows
        for rank, (lbl, binder) in enumerate(binder_dict.items(), start=1):
            writer.writerow([rank, binder.label, binder.sequence, binder.ipae, binder.rmsd, binder.plddt, binder.iptm])

# Main function to orchestrate the entire process
def main():
    # Parse command-line arguments
    args = parse_arguments()
    
    # Read CSV file
    csv_lines = read_csv_file(args.input_csv)
    
    # Process binders
    binders = []
    num_with_excessive_ala_stretches = 0
    
    # Iterate through CSV lines (skipping the header) and create Binder objects
    for i, line in enumerate(csv_lines[1:], start=1):
        binder = Binder(line, i)
        if binder.max_alanine_stretch_length() > args.max_ala_stretch:
            # Count binders with excessive alanine stretches
            num_with_excessive_ala_stretches += 1
        else:
            # Add valid binders to the list
            binders.append(binder)
    
    # Separate binders into sets based on various criteria
    set1, set2, set3, num_identical = process_binders(binders, args.min_plddt, args.similarity_cutoff, args.ipae_window_size)
    
    # Write output files
    write_fasta(args.output_fasta, set1, args.write_top_n)
    write_ranked_csv(args.output_csv, set1)
    
    # Print summary statistics
    print(f'\n{args.input_csv} contained {len(binders) + num_with_excessive_ala_stretches} binders.')
    print(f'{num_with_excessive_ala_stretches} had Ala stretches of > {args.max_ala_stretch}.')
    print(f'{len(set2)} of those remaining had plddt < {args.min_plddt}.')
    print(f'Among those with plddt >= {args.min_plddt} and acceptable Ala stretches, {num_identical} were redundant, and {len(set3)} were excluded because they had < {args.similarity_cutoff} differences from higher ranking seqs.')
    print(f'Final set after removing excessive Ala stretches, low plddt, redundant, and similar seqs: {len(set1)} seqs left.\n')
    print(f'Ranked CSV file has been written to: {args.output_csv}')
    print(f'FASTA file has been written to: {args.output_fasta}')

if __name__ == "__main__":
    # This block only executes if the script is run directly (not imported as a module)
    
    # Usage instructions
    usage = '''
    python v4_csv-to-fasta.py -i INPUT_CSV -f OUTPUT_FASTA -o OUTPUT_CSV [-n WRITE_TOP_N] [-p MIN_PLDDT] [-s SIMILARITY_CUTOFF] [-w IPAE_WINDOW_SIZE] [-a MAX_ALA_STRETCH]

    Options:
    -i, --input_csv: Path to input CSV file (required)
    -f, --output_fasta: Path to output FASTA file (optional, default: INPUT_CSV_ranked.fasta)
    -o, --output_csv: Path to output ranked CSV file (optional, default: INPUT_CSV_ranked.csv)
    -n, --write_top_n: Number of top sequences to write to FASTA (optional, if not provided, all sequences will be written in ranked order)
    -p, --min_plddt: Minimum pLDDT threshold (default: 70.0)
    -s, --similarity_cutoff: Sequence similarity cutoff (default: 18)
    -w, --ipae_window_size: IPAE window size (default: 1.0, if 0, seqs will be sorted only by IPAE)
    -a, --max_ala_stretch: Maximum allowed alanine stretch (default: 7)

    Seqs with < similarity_cutoff aa differences from seqs ranked above them will be removed.
    Seqs with alanine stretches longer than max_ala_stretch will be removed.
    A ranked CSV file will be output in addition to the FASTA file.
    If output files already exist, new files with timestamps will be generated.
    '''
    # Print the usage instructions
    print(usage)
    # Call the main function to start the program
    main()
