import regex as re
import gzip, os, glob
import pandas as pd
from pathlib import Path 
import matplotlib.pyplot as plt
import csv
import math
import subprocess

from itertools import zip_longest
import mobaseq.utils.logger as log
import mobaseq.process.tools as tools

# sgID matching function to tolerate one mismatch
def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip_longest(s1, s2))

def getsgID(sgIDDict, R1sgID="GTAAGGAG"):
    for gene, SGIDSeq in sgIDDict.items():
        if hamming_distance(R1sgID, SGIDSeq) <= 1:
            return(gene+","+str(hamming_distance(R1sgID, SGIDSeq)))
    return(R1sgID+",NA")

def countreads(read_one, read_two, sample_name, sgID_dict, min_length, max_length, all_start_with_G, out_dir, debug):
    '''
    This function reads the merged reads, merge reads together and get high quality ones
    
    File1: Forward read file 
    return the a dictionary of barcodes
    '''
    # 1. Read in both forward and reverse fastq files, get the second line of sequence
    # the sequencing quality is so high that there is no need to consider the phred score if the forward and reverse reads matches.
    # Open the two fastq files
    log.logit(f"Reading {read_one} and {read_two}")
    f1 = gzip.open(read_one,'rt')
    f2 = gzip.open(read_two,'rt')
    log.logit(f"Done reading.")

    sgIDBCdict = {}
    
    line_seq1 = 1
    
    while (line_seq1):
        # read in the files line by line
        # remove the first and third lines
        # store the second line into line_seq, this is the sequence
        # store the fourth line into line, this is the quality score
        
        # for R1 file
        line1 = f1.readline().rstrip() # skip the first line
        line_seq1 = f1.readline().rstrip() # line for sequencing quality
        line1 = f1.readline().rstrip() # skip the third line
        line_qua1 = f1.readline().rstrip() # get the sequencing quality
        
		# for R2 file
        line2 = f2.readline().rstrip()
        line_seq2 = f2.readline().rstrip()
        line2 = f2.readline().rstrip()
        line_qua2 = f2.readline().rstrip()
        
        line_seq2_complement=tools.revcom(line_seq2)
        
        # DEBUG: Prints full read sequence 
        #print(f"R1 sequence: {line_seq1}")
        #print(f"R2 sequence (complement): {line_seq2_complement}")
        
        # find the locator for R1 and R2
        # regex = re.compile('GTT' + (17N random BC) + 'ATGG' + (sgRNA/sgID) + 'GTTTAA') Moba 500
        if min_length > 8:
            # Use regex for guide RNAs longer than 8 if all sequences start with G
            if all_start_with_G:
                regexR1 = re.compile(f'.*GTT(.{{17}})ATG(.{{{min_length},{max_length}}})GTT(TAA){{e<2}}')  # 15 nt for alignment purpose
                regexR2 = re.compile(f'.*GTT(.{{17}})ATG(.{{{min_length},{max_length}}})GTT(TAA){{e<3}}')  # Looser criteria for R2
            # Regex for when gRNA sequences do not all start with G
            else:
                regexR1 = re.compile(f'.*GTT(.{{17}})ATGG(.{{{min_length},{max_length}}})GTT(TAA){{e<2}}')  # 15 nt for alignment purpose
                regexR2 = re.compile(f'.*GTT(.{{17}})ATGG(.{{{min_length},{max_length}}})GTT(TAA){{e<3}}')  # Looser criteria for R2
        else:
            # Original 3D Moba regex = re.compile(‘GA’ + (sgID) + (random_barcode) + ‘ATGCCCAAGAAG’)\n”,
                regexR1 = re.compile(r'GA(.{8})(GC.....TA.....GC.....TA.....GC)(ATGCCCA){e<2}')  # Allowing 2 errors in the ATGCCCA
                regexR2 = re.compile(r'A(.{8})(GC.....TA.....GC.....TA.....GC)A(TGCCCA){e<3}')  # Allowing 3 errors in the TGCCCA

        if regexR1.search(line_seq1) and regexR2.search(line_seq2_complement):
            # if we can see patterns in both reads
            k1 = regexR1.search(line_seq1)  # align R1
            if min_length > 8:
                R1sgID = k1.group(2)  # For longer guide RNAs, R1sgID is from group 2
                R1BC = k1.group(1)     # R1BC is from group 1
            else:
                R1sgID = k1.group(1)  # For shorter guide RNAs, R1sgID is from group 1
                R1BC = k1.group(2)    # R1BC is from group 2

            k2 = regexR2.search(line_seq2_complement)  # align R2
            if min_length > 8:
                R2sgID = k2.group(2)  # For longer guide RNAs, R2sgID is from group 2
                R2BC = k2.group(1)     # R2BC is from group 1
            else:
                R2sgID = k2.group(1)  # For shorter guide RNAs, R2sgID is from group 1
                R2BC = k2.group(2)    # R2BC is from group 2
            # combine reads only if barcode are correct from both ends
            # use the sequence of R1sgID as the sgID
            sgID = getsgID(sgID_dict, R1sgID=R1sgID)
            
            # if exists barcode and sgID combination
            if R1BC == R2BC and ('N' not in R1BC):
                myKey = sgID+","+R1BC
                if myKey in sgIDBCdict:
                    sgIDBCdict[myKey] += 1
                else:
                    sgIDBCdict[myKey] = 1
    f1.close()
    f2.close()

    # Write the dictionary to a file
    log.logit(f"Writing output to {sample_name}_MergeReadOut.csv.")
    with open(out_dir + "/" + sample_name + "_MergeReadOut.csv", "w") as f:
        for key in sorted(sgIDBCdict, key=sgIDBCdict.get, reverse=True):
            f.write(",".join(map(str, [key, sgIDBCdict[key]])) + "\n")
    log.logit(f"Finished writing output to {sample_name}_MergeReadOut.csv.")

def clean_barcodes(input_dir, out_dir, debug):
    merge_reads_csv_files = tools.get_list_of_files(input_dir, "merge_reads_csv", debug)
    # We need to make the Mobaseq_Sample_Info.txt
    df = pd.DataFrame([
        (csv, name) for csv, name in merge_reads_csv_files
    ], columns=['MergeReadOut', 'Sample'])
    df.to_csv(out_dir + "/Mobaseq_Sample_Info.txt", index=False, sep = "\t")
    current_dir = Path(__file__).parent
    rscript_path = current_dir / 'filter_fake_bc.r'
    if not rscript_path.exists():
        raise FileNotFoundError(f"R script not found at {rscript_path}")
        
    log.logit("Running R script to filter fake barcodes...")
    cmd = ['Rscript', str(rscript_path)] + [str(out_dir), str(out_dir + "/Mobaseq_Sample_Info.txt")]
    try:
        # Run R script
        result = subprocess.run(
            cmd,
            cwd=current_dir,
            check=True,
            capture_output=True,
            text=True
        )
        log.logit("R script completed successfully.")
        # Remove the temporary file
        os.remove(out_dir + "/Mobaseq_Sample_Info.txt")
        return result.stdout
    except subprocess.CalledProcessError as e:
        err_msg = f"Error running R script: {e}"
        err_msg += f"R Error Output: {e.stderr}"
        raise RuntimeError(err_msg)

# Summarizes information generated from the countreads function
def summarize(input_dir, sgIDDict, out_dir, debug):
    merge_reads_csv_files = glob.glob(os.path.join(input_dir, "*_MergeReadOut.csv"))
    merge_reads_csv_files.sort()

    # Arrays to store total reads, unmapped reads, and column names
    total_reads_array = []
    unmapped_reads_array = []
    names_array = []

    # Dictionary to store total counts and number of rows for each sgID
    sgID_total_reads = {}
    sgID_row_counts = {}

    # Process each MergeReadOut file
    for idx, merge_readout_file in enumerate(merge_reads_csv_files):
        merge_readout_file = Path(merge_readout_file)
        # Check if the file is empty (file size is 0)
        if Path(merge_readout_file).stat().st_size == 0:
            log.logit(f"File '{merge_readout_file.name}' is empty, skipping.")
            total_reads_array.append(0)
            unmapped_reads_array.append(0)
            names_array.append(merge_readout_file.name.replace('_MergeReadOut.csv', ''))
            continue

        df = pd.read_csv(merge_readout_file, sep=",", header=None)

        # Check if DataFrame is empty (e.g., after reading a malformed file)
        if df.empty:
            print(f"File '{merge_readout_file.name}' contains no valid data, skipping.")
            total_reads_array.append(0)
            unmapped_reads_array.append(0)
            names_array.append(merge_readout_file.name.replace('_MergeReadOut.csv', ''))
            continue

        # Ensure columns exist in the DataFrame
        df.columns = ['sgID', 'dist', 'BC', 'reads']
        
        # Count the total reads and reads that are unmapped (not filtering by dist)
        total_reads = df['reads'].sum()
        unmapped_reads = df.loc[~df['sgID'].isin(sgIDDict.keys()), 'reads'].sum()
        
        # Store the total number of reads in the arrays
        total_reads_array.append(total_reads)
        unmapped_reads_array.append(unmapped_reads)
        names_array.append(merge_readout_file.name.replace('_MergeReadOut.csv', ''))
        
        # Filter rows where dist = 0 for sgID calculations
        df_dist_zero = df[df['dist'] == 0]
        
        # Update the total reads and row counts for each sgID
        for sgID, group in df_dist_zero.groupby('sgID'):
            if sgID not in sgID_total_reads:
                sgID_total_reads[sgID] = [0] * len(merge_reads_csv_files)
                sgID_row_counts[sgID] = [0] * len(merge_reads_csv_files)
            
            sgID_total_reads[sgID][idx] = group['reads'].sum()
            sgID_row_counts[sgID][idx] = len(group)
        
        # Ensure all sgIDs have an entry for the current file
        for sgID in sgID_total_reads:
            if len(sgID_total_reads[sgID]) <= idx:
                sgID_total_reads[sgID].append(0)
                sgID_row_counts[sgID].append(0)

    # Convert lists to DataFrames
    total_reads_df = pd.DataFrame(sgID_total_reads, index=names_array).T
    row_counts_df = pd.DataFrame(sgID_row_counts, index=names_array).T

    # Calculate relative reads and counts by dividing by the column total
    relative_total_reads_df = total_reads_df.div(total_reads_df.sum(axis=0), axis=1)
    relative_row_counts_df = row_counts_df.div(row_counts_df.sum(axis=0), axis=1)

    # Create a list of DataFrames to save
    data_frames = [
        ('total_reads_per_sgid', total_reads_df),
        ('barcodes_per_sgid', row_counts_df),
        ('relative_reads_per_sgid', relative_total_reads_df),
        ('relative_barcodes_per_sgid', relative_row_counts_df)
    ]

    # Save each DataFrame to a CSV file
    for name, df in data_frames:
        csv_filename = os.path.join(out_dir, f'{name}.csv')
        df.to_csv(csv_filename, index_label='sgID')
        log.logit(f"CSV file '{csv_filename}' has been created successfully.")

def filter_fake_bc(input_dir, out_dir, debug):
    pass

# Function to extract tissue and number
def sort_key(key):
    # Define tissue sorting order
    tissue_order = {"Liver": 1, "Lung": 2, "Brain": 3, "Blood": 4, "WB": 5, "BM": 6, "preTran": 7}
    match = re.match(r"(\d+)-(\w+)", key)  # Extract number and tissue type
    if match:
        num = int(match.group(1))  # Extract numeric part
        tissue = match.group(2)  # Extract tissue
        tissue_rank = tissue_order.get(tissue, 99)  # Assign tissue rank, default to 99 if not found
        return (tissue_rank, num)  # Sort by tissue, then by number
    return (99, 99999)  # Default sorting for unexpected values

def plot_mapped_and_unmapped(input_dir, sgIDDict, out_dir, debug):
    merge_reads_csv_files = glob.glob(os.path.join(input_dir, "*_MergeReadOut.csv"))
    merge_reads_csv_files.sort()

    # Arrays to store total reads, unmapped reads, and column names
    total_reads_array = []
    unmapped_reads_array = []
    names_array = []

    # Process each MergeReadOut file
    for idx, merge_readout_file in enumerate(merge_reads_csv_files):
        merge_readout_file = Path(merge_readout_file)
        # Check if the file is empty (file size is 0)
        if Path(merge_readout_file).stat().st_size == 0:
            log.logit(f"File '{merge_readout_file.name}' is empty, skipping.")
            total_reads_array.append(0)
            unmapped_reads_array.append(0)
            names_array.append(merge_readout_file.name.replace('_MergeReadOut.csv', ''))
            continue

        df = pd.read_csv(merge_readout_file, sep=",", header=None)

        # Check if DataFrame is empty (e.g., after reading a malformed file)
        if df.empty:
            print(f"File '{merge_readout_file.name}' contains no valid data, skipping.")
            total_reads_array.append(0)
            unmapped_reads_array.append(0)
            names_array.append(merge_readout_file.name.replace('_MergeReadOut.csv', ''))
            continue

        # Ensure columns exist in the DataFrame
        df.columns = ['sgID', 'dist', 'BC', 'reads']
        
        # Count the total reads and reads that are unmapped (not filtering by dist)
        total_reads = df['reads'].sum()
        unmapped_reads = df.loc[~df['sgID'].isin(sgIDDict.keys()), 'reads'].sum()
        
        # Store the total number of reads in the arrays
        total_reads_array.append(total_reads)
        unmapped_reads_array.append(unmapped_reads)
        names_array.append(merge_readout_file.name.replace('_MergeReadOut.csv', ''))
    
    # Import data and format as percentages
    sgIDs = names_array

    # Initialize mapped and unmapped percentage lists
    mapped_percentages = []
    unmapped_percentages = []

    # Process unmapped reads and total reads arrays
    for unmapped, total in zip(unmapped_reads_array, total_reads_array):
        if total == 0:  # Handle empty or missing data
            mapped_percentages.append(0)
            unmapped_percentages.append(0)
        else:
            unmapped_fraction = unmapped / total
            unmapped_percentage = unmapped_fraction * 100
            mapped_percentage = 100 - unmapped_percentage
            mapped_percentages.append(mapped_percentage)
            unmapped_percentages.append(unmapped_percentage)
            
    # Define the file name for the CSV
    csv_filename = os.path.join(out_dir, "mapped_percentages.csv")

    # Combine categories and percentages into a list of tuples
    #data = list(zip(sgIDs, mapped_percentages, unmapped_percentages))

    # Combine categories and percentages into a list of tuples
    data = list(zip(sgIDs, mapped_percentages, unmapped_percentages))

    # Sort the data by tissue type first, then by number
    sorted_data = sorted(data, key=lambda x: sort_key(x[0]))

    # Unzip sorted data back into separate lists
    sgIDs, mapped_percentages, unmapped_percentages = zip(*sorted_data)
    # Write data to CSV file
    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header
        writer.writerow(["sgID", "Mapped Percentage", "Unmapped Percentage"])
        
        # Write data rows
        for row in data:
            writer.writerow(row)
            
    # Create a bar plot
    fig, ax = plt.subplots(figsize=(22, 16))

    # Plot the percentage of mapped 
    bars1 = ax.bar(sgIDs, mapped_percentages, label='Mapped', color='navy')

    # Plot the percentage of unmapped on top of the first one
    bars2 = ax.bar(sgIDs, unmapped_percentages, label='Unmapped', color='lightblue', bottom=mapped_percentages)

    # Outline bars with mapped percentage under 30% in red
    for bar, mapped in zip(bars1, mapped_percentages):
        if mapped < 30:
            bar.set_edgecolor('red')
            bar.set_linewidth(2)
            
    # Display the percentage values on top of the bars
    for bar in bars1:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval+0.1, f'{int(yval)}%', ha='center', va='bottom', color="black")

    for bar1, bar2 in zip(bars1, bars2):
        yval = bar2.get_height()
        ax.text(bar2.get_x() + bar2.get_width()/2, (yval+bar1.get_height()+0.1), f'{int(yval)}%', ha='center', va='bottom', color="black")

    # Set labels and title
    ax.set_ylabel('Percentage')
    ax.set_title('Percentage of unmapped and mapped reads for each sample')
    ax.legend()

    # Reorder legend items to have "Mapped" at the bottom
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])

    # Show the plot
    fig.autofmt_xdate()
    plt.tight_layout()
    
    fig.savefig(os.path.join(out_dir, 'mapped_percentages.pdf'), format='pdf')

# Function to filter DataFrame and extract values greater than 0.1
def filter_significant_entries(df, value_threshold):
    significant_entries = []
    for sgID in df.index:
        for sample in df.columns:
            value = df.at[sgID, sample]
            if value > value_threshold:
                significant_entries.append({'Sample': sample, 'sgID': sgID, 'Value': value})
    return significant_entries

def plot_per_sgid(reads_csv, barcodes_csv, rel_reads_csv, rel_barcodes_csv, sgIDs, out_dir, no_dummy, debug):
    log.logit(f"Plotting the relative distribution of reads and barcodes per sgID by sample.")
    # Assuming sgIDDict is a dictionary containing sgIDs, define the dynamic threshold
    dynamic_threshold = 1/(2*math.sqrt(len(sgIDs)))

    log.logit(f"Reading in CSV Files...")
    # Read CSV files into pandas DataFrames
    df_total_reads = pd.read_csv(reads_csv).set_index('sgID').transpose()
    df_unique_barcodes = pd.read_csv(barcodes_csv).set_index('sgID').transpose()
    df_relative_reads = pd.read_csv(rel_reads_csv).set_index('sgID').transpose()
    df_relative_barcodes = pd.read_csv(rel_barcodes_csv).set_index('sgID').transpose()

    if no_dummy:
        log.logit(f"Filtering out sgDummy from the DataFrames...")
        # Pre-filter out rows with sgID == 'sgDummy'
        cols_to_exclude = ['sgDummy', 'SpikeIn']

        df_total_reads = df_total_reads.loc[:, ~df_total_reads.columns.isin(cols_to_exclude)]
        df_unique_barcodes = df_unique_barcodes.loc[:, ~df_unique_barcodes.columns.isin(cols_to_exclude)]
        df_relative_reads = df_relative_reads.loc[:, ~df_relative_reads.columns.isin(cols_to_exclude)]
        df_relative_barcodes = df_relative_barcodes.loc[:, ~df_relative_barcodes.columns.isin(cols_to_exclude)]

    # Rename preTran1 to pretransplantation
    df_total_reads.rename(columns={'preTran1': 'pretransplantation'}, inplace=True)
    df_unique_barcodes.rename(columns={'preTran1': 'pretransplantation'}, inplace=True)
    df_relative_reads.rename(columns={'preTran1': 'pretransplantation'}, inplace=True)
    df_relative_barcodes.rename(columns={'preTran1': 'pretransplantation'}, inplace=True)

    # Sort the rows based on the index values (tissue names)
    df_total_reads = df_total_reads.sort_index(key=lambda idx: [sort_key(x) for x in idx])
    df_unique_barcodes = df_unique_barcodes.sort_index(key=lambda idx: [sort_key(x) for x in idx])
    df_relative_reads = df_relative_reads.sort_index(key=lambda idx: [sort_key(x) for x in idx])
    df_relative_barcodes = df_relative_barcodes.sort_index(key=lambda idx: [sort_key(x) for x in idx])

    if no_dummy:
        log.logit(f"Recalculating relative reads and barcodes after filtering out sgDummy...")
        # Recalculate total reads after filtering out sgDummy
        df_total_reads_new = df_total_reads.sum(axis=1)

        # Recalculate relative reads based on the new total reads
        df_relative_reads = df_total_reads.div(df_total_reads_new, axis=0)

        # Recalculate relative barcodes based on the new total reads (optional depending on your case)
        df_relative_barcodes_new_total = df_unique_barcodes.sum(axis=1)
        df_relative_barcodes = df_unique_barcodes.div(df_relative_barcodes_new_total, axis=0)

        # Sort columns in ascending order for the legend
        sorted_columns = sorted(df_total_reads.columns)

        # Sort columns in descending order for plotting
        df_total_reads = df_total_reads.reindex(sorted(df_total_reads.columns, reverse=True), axis=1)
        df_unique_barcodes = df_unique_barcodes.reindex(sorted(df_unique_barcodes.columns, reverse=True), axis=1)
        df_relative_reads = df_relative_reads.reindex(sorted(df_relative_reads.columns, reverse=True), axis=1)
        df_relative_barcodes = df_relative_barcodes.reindex(sorted(df_relative_barcodes.columns, reverse=True), axis=1)

    if no_dummy:
        df_relative_reads.to_csv(os.path.join(out_dir, 'TEST.csv'), index_label='Sample')
        significant_relative_reads = filter_significant_entries(df_relative_reads, 0.1)
        significant_relative_barcodes = filter_significant_entries(df_relative_barcodes, 0.1)
        all_significant_entries = significant_relative_reads + significant_relative_barcodes
        df_significant = pd.DataFrame(all_significant_entries)
        csv_filename_significant = os.path.join(out_dir, 'significant_relative_values_filtered_noDummy.csv')
    else:
        significant_relative_reads = filter_significant_entries(df_relative_reads, 0.02)
        significant_relative_barcodes = filter_significant_entries(df_relative_barcodes, 0.02)
        all_significant_entries = significant_relative_reads + significant_relative_barcodes
        df_significant = pd.DataFrame(all_significant_entries)
        csv_filename_significant = os.path.join(out_dir, 'significant_relative_values_filtered.csv')
    # Write the DataFrame to CSV
    df_significant.to_csv(csv_filename_significant, index=False)

    # Plotting all figures side by side
    fig, axs = plt.subplots(2, 2, figsize=(24, 16))  # 2 rows, 2 columns
    colors = plt.cm.tab10.colors[:len(df_unique_barcodes.columns)]  # Using Tab10 colormap

    # Initialize bottom position for stacking
    bottom1 = bottom2 = bottom3 = bottom4 = None

    # Plot stacked bars for total reads per sgID
    for sgID in df_total_reads.columns:
        if bottom1 is None:
            axs[0, 0].bar(df_total_reads.index, df_total_reads[sgID], label=sgID)
            bottom1 = df_total_reads[sgID]
        else:
            axs[0, 0].bar(df_total_reads.index, df_total_reads[sgID], bottom=bottom1, label=sgID)
            bottom1 += df_total_reads[sgID]

    # Customize plot for total reads per sgID
    axs[0, 0].set_xlabel('Samples')
    axs[0, 0].set_ylabel('Counts')
    axs[0, 0].set_title('Total reads per sgID by Sample')
    axs[0, 0].tick_params(axis='x', rotation=75)

    # Plot stacked bars for unique barcodes per sgID
    for sgID in df_unique_barcodes.columns:
        if bottom2 is None:
            axs[0, 1].bar(df_unique_barcodes.index, df_unique_barcodes[sgID], label=sgID)
            bottom2 = df_unique_barcodes[sgID]
        else:
            axs[0, 1].bar(df_unique_barcodes.index, df_unique_barcodes[sgID], bottom=bottom2, label=sgID)
            bottom2 += df_unique_barcodes[sgID]

    # Customize plot for unique barcodes per sgID
    axs[0, 1].set_xlabel('Samples')
    axs[0, 1].set_ylabel('Counts')
    axs[0, 1].set_title('Unique barcodes per sgID by Sample')
    axs[0, 1].tick_params(axis='x', rotation=45)

    # Plot stacked bars for relative reads per sgID
    for sgID in df_relative_reads.columns:
        if bottom3 is None:
            axs[1, 0].bar(df_relative_reads.index, df_relative_reads[sgID], label=sgID)
            bottom3 = df_relative_reads[sgID]
        else:
            axs[1, 0].bar(df_relative_reads.index, df_relative_reads[sgID], bottom=bottom3, label=sgID)
            bottom3 += df_relative_reads[sgID]
        # Add label if value is greater than 0.1
        for i, value in enumerate(df_relative_reads[sgID]):
            if value > dynamic_threshold if no_dummy else 0.1:
                axs[1, 0].text(i, bottom3[i] - value / 2, sgID, ha='center', va='center',
                            fontsize=8, color='black', rotation=90)

    # Customize plot for relative reads per sgID
    axs[1, 0].set_xlabel('Samples')
    axs[1, 0].set_ylabel('Relative Counts')
    axs[1, 0].set_title('Relative reads per sgID by Sample')
    axs[1, 0].tick_params(axis='x', rotation=75)

    #sgID_to_gene = df_sgID.set_index('sgID')['sgID'].to_dict()  # Create the mapping dictionary

    # Plot stacked bars for relative barcodes per sgID
    for sgID in df_relative_barcodes.columns:
        if bottom4 is None:
            axs[1, 1].bar(df_relative_barcodes.index, df_relative_barcodes[sgID], label=sgID)
            bottom4 = df_relative_barcodes[sgID]
        else:
            axs[1, 1].bar(df_relative_barcodes.index, df_relative_barcodes[sgID], bottom=bottom4, label=sgID)
            bottom4 += df_relative_barcodes[sgID]
        # Add label if value is greater than 0.1
        for i, value in enumerate(df_relative_barcodes[sgID]):
            if value > dynamic_threshold:
                axs[1, 1].text(i, bottom4[i] - value / 2, sgID, ha='center', va='center',
                            fontsize=8, color='black', rotation=90)

    # Customize plot for relative barcodes per sgID
    axs[1, 1].set_xlabel('Samples')
    axs[1, 1].set_ylabel('Relative Counts')
    axs[1, 1].set_title('Relative barcodes per sgID by Sample')
    axs[1, 1].tick_params(axis='x', rotation=75)

    # Create a single legend across the top of the plot
    handles, labels = axs[0, 0].get_legend_handles_labels()
    # Determine the number of columns in the legend (cap at 18)
    num_cols = 18 if len(labels) > 18 else len(labels)
    # Calculate the number of rows in the legend
    legend_rows = (len(labels) + num_cols - 1) // num_cols  # Ceiling division to get the number of rows
    # Create the legend
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.), ncol=num_cols, title='sgIDs',
           fontsize=10, title_fontsize='13', bbox_transform=fig.transFigure)
    
    # Adjust layout and save figure
    # Use fig.subplots_adjust to move the plot area down based on the number of legend rows
    legend_height_adjustment = 0.05 * legend_rows  # Adjust factor (0.05) based on visual spacing

    # Set minimum top margin to avoid it overlapping with the bottom
    min_top_margin = 0.55  # You can adjust this as needed
    top_margin = max(0.85 - legend_height_adjustment, min_top_margin)

    # Set a lower bottom margin to make more space for the plot
    bottom_margin = 0.1  # Fixed lower bound for bottom margin

    # Ensure top > bottom to avoid the error
    if top_margin - bottom_margin < 0.2:  # Ensure a minimum space between the two
        bottom_margin = max(0.05, top_margin - 0.2)  # Adjust bottom margin accordingly

    # Apply the margin adjustments
    fig.subplots_adjust(top=top_margin, bottom=bottom_margin)  # Dynamically adjust both top and bottom margins

    # Adjust layout and save figure
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to make space for the legend

    # Define output CSV file path
    if no_dummy:
        plt.savefig(os.path.join(out_dir, 'sgID_distribution_with_relative_reads_and_barcodes_filtered_noDummy.pdf'), 
                format='pdf', bbox_inches='tight')
        # Show the combined plot
        
    else:
        plt.savefig(os.path.join(out_dir, 'sgID_distribution_with_relative_reads_and_barcodes_filtered.pdf'), 
                format='pdf', bbox_inches='tight')
        # Show the combined plot
        
