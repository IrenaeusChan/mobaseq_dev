import regex as re
import os, sys
import pandas as pd
import gzip
import glob
from pathlib import Path 
import mobaseq.utils.logger as log
import numpy as np
from concurrent.futures import ThreadPoolExecutor

'''
This file contains various tools and functions that are used in the processing of the raw MOBAseq files.
'''

def ensure_abs_path(path):
    path = Path(path) # Ensure path is an absolute path
    if not path.is_absolute(): path = path.resolve()
    path.mkdir(parents=True, exist_ok=True) # Ensure the output directory exists
    return path

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    return ''.join([complement.get(base, 'N') for base in seq][::-1])

def get_sample_name(file, file_type, debug):
    dispatch = {
        'merge_read_csv' : get_sample_name_merge_read_csv,
        'barcode_clean_txt' : get_sample_name_barcode_clean_txt,
        'spike_count' : get_sample_name_spike_count
    }
    function = dispatch[file_type]
    return function(file, debug)

def get_sample_name_merge_read_csv(file, debug):
    log.logit(f"Processing {os.path.basename(file)}.", color="green")
    sample_name = os.path.basename(file).split('_MergeReadOut.csv')[0]
    log.logit(f"Finished processing {os.path.basename(file)}.")
    return sample_name

def get_sample_name_barcode_clean_txt(file, debug):
    log.logit(f"Processing {os.path.basename(file)}.", color="green")
    sample_name = os.path.basename(file).split('_BarcodeClean.txt')[0]
    log.logit(f"Finished processing {os.path.basename(file)}.")
    return sample_name

def get_sample_name_spike_count(file, debug):
    log.logit(f"Processing {os.path.basename(file)}.", color="green")
    sample_name = os.path.basename(file).split('_SpikeInCounts.txt')[0]
    log.logit(f"Finished processing {os.path.basename(file)}.")
    return sample_name

def process_fastq_files(read_one, read_two, debug):
    log.logit(f"Processing {os.path.basename(read_one)} and {os.path.basename(read_two)}.", color="green")
    if debug: log.logit(f"Determining sample name.", color="yellow")
    pattern = r"(.*)_R[12]\.(fastq|fq)(\.gz)?$" # Pattern to match the sample name, may not be Universal. Assumes anything before _R1 or _R2 is the sample name
    match1 = re.match(pattern, os.path.basename(read_one))
    match2 = re.match(pattern, os.path.basename(read_two))
    if match1.group(1) == match2.group(1) and match1 != None and match2 != None:
        sample_name = match1.group(1)
    else:
        err_msg = f"ERROR: FASTQ files names do not match. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    log.logit(f"Finished processing FASTQ files.")
    return sample_name

def convert_sgID_format(sgID):
    '''Convert sgID format, for example 'sgSafe30-MMU' to 'Safe30mmu-sg1'.'''
    match = re.match(r"sg([A-Za-z0-9]+)-?(mmu|MMU)?$", sgID)
    if match != None:
        gene = match.group(1)  # Extract the gene part
        species = match.group(2).lower() if match.group(2) else 'mmu'  # Default species to 'mmu'
        return f"{gene}{species}-sg1"
    return sgID

def split_sgID_name(sgID):
    '''Split sgID into Gene, Species, and sgRNA parts.'''
    if sgID.startswith("Dummy"): # Special case for "Dummy"
        return "Dummy", "mmu", sgID.split('-')[-1]
    match = re.match(r"([A-Za-z0-9]+?)(mmu|hsa)?-(sg\d+)", sgID) # General regex pattern for parsing the sgID
    if match:
        gene = match.group(1)  # Gene name
        species = match.group(2) if match.group(2) else ''  # Species code
        sgRNA = match.group(3)  # sgRNA part (e.g., sg1)
        return gene, species, sgRNA
    return sgID, '', ''  # Default return if no match

def revert_sgID_format(sgID):
    '''Revert sgID format, for example 'Safe30mmu-sg1' to 'sgSafe30-MMU'.'''
    if sgID == 'Dummymmu-sg1':
        return 'Dummy'  # Revert to original
    return sgID  # Return unchanged for others

def process_sgid_file(sgid_file, debug):
    log.logit(f"Processing {os.path.basename(sgid_file)}.", color="green")
    log.logit(f"Reading in sgID file.")
    sgID_df = pd.read_excel(sgid_file, header=None, names=['sgID', 'gRNA_seq']) # Load the sgID file into a DataFrame
    log.logit(f"Loaded sgID file: {sgid_file}")
    sgIDs = {'sgDummy': 'AGAGACGCTCGAGCGTCTCT'} # Create the sgID Dictionary
    sgIDs.update(sgID_df.set_index('sgID')['gRNA_seq'].to_dict()) # Populate the dictionary with values from the DataFrame
    sgID_dummy_df = pd.DataFrame([['sgDummy', 'AGAGACGCTCGAGCGTCTCT']], columns=['sgID', 'gRNA_seq']) # Create a DataFrame with a single row for 'sgDummy'
    sgID_df = pd.concat([sgID_dummy_df, sgID_df]).reset_index(drop=True) # Concatenate the dummy DataFrame with the original DataFrame    
    log.logit(f"Converting sgID format.")
    sgID_df['sgID'] = sgID_df['sgID'].apply(convert_sgID_format)
    log.logit(f"Splitting sgID names.")
    sgID_df[['sgID', 'Species', 'sgRNA']] = sgID_df['sgID'].apply(lambda x: pd.Series(split_sgID_name(x)))
    #log.logit(f"Reverting sgID format.")
    # df_sgID['sgID'] = df_sgID['sgID'].apply(revert_sgID_format)
    log.logit(f"Finished processing sgID Excel file.")
    return sgIDs

def process_spike_ins(spike_ins, debug):
    log.logit(f"Processing {os.path.basename(spike_ins)}.", color="green")
    log.logit(f"Reading in spike-in file.")
    spike_in_df = pd.read_csv(spike_ins)
    # Check if all sequences are 17 NT in length
    if not all(len(seq) == 17 for seq in spike_in_df['sequence']):
        err_msg = f"ERROR: Not all spike-in sequences are 17 nucleotides in length. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    log.logit(f"Finished processing spike-in information file.")
    return spike_in_df

def get_info_about_sgID_file(sgIDs, debug):
    seq_len = [len(seq) for seq in sgIDs.values()] # Calculate the minimum and maximum sequence lengths
    min_length = min(seq_len)
    max_length = max(seq_len)
    all_start_with_G = all(seq.startswith('G') for seq in sgIDs.values())
    log.logit(f"All sequences start with 'G': {all_start_with_G}")
    log.logit(f"Minimum sequence length: {min_length}")
    log.logit(f"Maximum sequence length: {max_length}")
    return min_length, max_length, all_start_with_G

def read_sequences_from_fastq(file_path, debug):
    log.logit(f"Reading sequences from {os.path.basename(file_path)}.", color="green")
    sequences = []
    # Determining file size to determine if it's better to use pandas or just read line by line
    file_size = os.path.getsize(file_path)

    # Open the file, using gzip if the file is compressed
    if file_path.endswith('.gz'): open_func = gzip.open
    else: open_func = open

    if file_size > 1 * 1024 * 1024 * 1024:  # If file size is greater than 1 GB
        if debug: log.logit(f"File size is greater than 1 GB. Using pandas to read the FASTQ file.", color="yellow")
        # Read the FASTQ file using pandas
        with open_func(file_path, 'rt') as f:
            lines = f.read().splitlines()
        # Extract sequences (every 4th line starting from the second line)
        sequences = [lines[i] for i in range(1, len(lines), 4)]
        if debug: log.logit(f"Finished reading sequences from {os.path.basename(file_path)}.", color="yellow")
    else:
        if debug: log.logit(f"File size is less than 1 GB. Reading the FASTQ file line by line.", color="yellow")
        # Read the FASTQ file line by line
        with open_func(file_path, 'rt') as f:
            line_number = 0
            for line in f:
                line_number += 1
                # Sequence lines are every 4th line starting from the second line
                if line_number % 4 == 2:
                    sequences.append(line.strip())
        if debug: log.logit(f"Finished reading sequences from {os.path.basename(file_path)}.", color="yellow")
    # NumPy Array is faster and better for vectorized operations
    #if debug: log.logit(f"Converting sequences to NumPy array.", color="yellow")
    #sequences = np.array(sequences)
    #sequences = pd.Series(sequences).values
    #if debug: log.logit(f"Finished converting sequences to NumPy array.", color="yellow")
    return sequences

# TODO: Attemping to introduce multithreading to read the sequences from the FASTQ file
def read_sequences_from_fastq_read_two(file_path):
    sequences = read_sequences_from_fastq(file_path)
    with ThreadPoolExecutor() as executor:
        rev_comps = list(executor.map(reverse_complement, sequences))
    return rev_comps

def get_list_of_files(input_dir, what_file, debug):
    dispatch = {
        'fastq' : get_list_of_fastqs,
        'merge_reads_csv' : get_list_of_merge_reads_csv,
        'barcode_clean_txt' : get_list_of_barcode_clean_txt,
        'cell_count' : get_list_of_cell_count
    }
    function = dispatch[what_file]
    return function(input_dir, debug)

def get_list_of_fastqs(input_dir, debug):
    log.logit(f"Getting list of FASTQ files from {input_dir}.", color="green")
    # Get a list of all FASTQ files in the input directory
    fastq_files = glob.glob(os.path.join(input_dir, "*.fastq")) + \
                    glob.glob(os.path.join(input_dir, "*.fq")) + \
                    glob.glob(os.path.join(input_dir, "*.fastq.gz")) + \
                    glob.glob(os.path.join(input_dir, "*.fq.gz"))
    # Check if the list is empty
    if len(fastq_files) == 0:
        err_msg = f"ERROR: No FASTQ files found in {input_dir}. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    log.logit(f"Found {len(fastq_files)} FASTQ files in {input_dir}.")
    # Need to pair the FASTQ files
    fastq_files.sort()
    paired_files = []
    # We are going to assume the same FASTQ files will have the same name except for the R1 and R2
    # therefore we will pair them together by sorting them and organizing them in pairs
    for i in range(0, len(fastq_files), 2):
        sample_name = process_fastq_files(fastq_files[i], fastq_files[i+1], debug)
        paired_files.append((fastq_files[i], fastq_files[i+1], sample_name))
    log.logit(f"Paired {len(paired_files)} FASTQ files.")
    return paired_files

def get_list_of_merge_reads_csv(input_dir, debug):
    log.logit(f"Getting list of MergeReadOut files from {input_dir}.", color="green")
    merge_read_csv_files = glob.glob(os.path.join(input_dir, "*_MergeReadOut.csv"))
    # Check if the list is empty
    if len(merge_read_csv_files) == 0:
        err_msg = f"ERROR: No MergeReadOut files found in {input_dir}. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    log.logit(f"Found {len(merge_read_csv_files)} MergeReadOut files in {input_dir}.")
    merge_read_csv_files.sort()
    files_and_sample_names = []
    for merge_read_csv in merge_read_csv_files:
        files_and_sample_names.append((merge_read_csv, get_sample_name(merge_read_csv, "merge_read_csv", debug)))
    return files_and_sample_names

def get_list_of_barcode_clean_txt(input_dir, debug):
    log.logit(f"Getting list of BarcodeClean files from {input_dir}.", color="green")
    barcode_clean_txt_files = glob.glob(os.path.join(input_dir, "*_BarcodeClean.txt"))
    # Check if the list is empty
    if len(barcode_clean_txt_files) == 0:
        err_msg = f"ERROR: No BarcodeClean files found in {input_dir}. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    log.logit(f"Found {len(barcode_clean_txt_files)} BarcodeClean files in {input_dir}.")
    barcode_clean_txt_files.sort()
    files_and_sample_names = []
    for barcode_clean_txt in barcode_clean_txt_files:
        files_and_sample_names.append((barcode_clean_txt, get_sample_name(barcode_clean_txt, "barcode_clean_txt", debug)))
    return files_and_sample_names

def get_list_of_cell_count(input_dir, debug):
    log.logit(f"Getting list of UnfilteredSampleInfo, FilteredSampleInfo, and SpikeInCount files from {input_dir}.", color="green")
    unfiltered_sample_info_files = glob.glob(os.path.join(input_dir, "*_UnfilteredSampleInfo.txt"))
    filtered_sample_info_files = glob.glob(os.path.join(input_dir, "*_FilteredSampleInfo.txt"))
    spike_count_files = glob.glob(os.path.join(input_dir, "*_SpikeInCounts.txt"))
    # Check if the list is empty
    if len(unfiltered_sample_info_files) == 0:
        err_msg = f"ERROR: No UnfilteredSampleInfo files found in {input_dir}. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    log.logit(f"Found {len(unfiltered_sample_info_files)} UnfilteredSampleInfo files in {input_dir}.")
    if len(filtered_sample_info_files) == 0:
        err_msg = f"ERROR: No FilteredSampleInfo files found in {input_dir}. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    log.logit(f"Found {len(filtered_sample_info_files)} FilteredSampleInfo files in {input_dir}.")
    if len(spike_count_files) == 0:
        err_msg = f"ERROR: No SpikeInCount files found in {input_dir}. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    log.logit(f"Found {len(spike_count_files)} SpikeInCount files in {input_dir}.")
    unfiltered_sample_info_files.sort()
    filtered_sample_info_files.sort()
    spike_count_files.sort()
    files_and_sample_names = []
    for unfiltered, filtered, spike_count in zip(unfiltered_sample_info_files, filtered_sample_info_files, spike_count_files):
        files_and_sample_names.append(((unfiltered, filtered, spike_count), get_sample_name(spike_count, "spike_count", debug)))
    return files_and_sample_names

def process_input_files(input_files_dir, debug):
    # Look for sgid_file, spike_ins, and library_info
    sgid_file = glob.glob(os.path.join(input_files_dir, "*sgID*.xls*"))
    spike_ins = glob.glob(os.path.join(input_files_dir, "B*key*.csv"))
    library_info = glob.glob(os.path.join(input_files_dir, "CountingTable.xls*"))
    # Check if the list is empty
    if len(sgid_file) == 0:
        err_msg = f"ERROR: No sgID file found in {input_files_dir}. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    if len(spike_ins) == 0:
        err_msg = f"ERROR: No spike-in file found in {input_files_dir}. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    if len(library_info) == 0:
        err_msg = f"ERROR: No library information file found in {input_files_dir}. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    sgid_file = sgid_file[0]
    spike_ins = spike_ins[0]
    library_info = library_info[0]
    log.logit(f"Found sgID file: {sgid_file}")
    log.logit(f"Found spike-in file: {spike_ins}")
    log.logit(f"Found library information file: {library_info}")
    return sgid_file, spike_ins, library_info