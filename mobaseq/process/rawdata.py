import os, sys
import regex as re
import numpy as np
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
from sklearn.linear_model import LinearRegression
from clint.textui import indent
import traceback
import mobaseq.utils.logger as log
import mobaseq.process.tools as tools
import mobaseq.qc.plot as plotting

import signal
from functools import wraps
from contextlib import contextmanager

@contextmanager
def timeout_handler(seconds, description=""):
    def timeout_signal(signum, frame):
        raise TimeoutError(f"Process timed out after {seconds}s: {description}")

    # Set signal handler
    signal.signal(signal.SIGALRM, timeout_signal)
    signal.alarm(seconds)
    
    try:
        yield
    finally:
        signal.alarm(0)

def vectorized_hamming_distance(s1_array, s2_array):
    """Calculate the Hamming distance between arrays of sequences."""
    max_len = max(max(len(s1) for s1 in s1_array), max(len(s2) for s2 in s2_array)) # Get the maximum length of the sequences
    s1_matrix = np.array([list(s1.ljust(max_len, 'N')) for s1 in s1_array]) # Pad the sequences with 'N' to the maximum length
    s2_matrix = np.array([list(s2.ljust(max_len, 'N')) for s2 in s2_array])
    return np.sum(s1_matrix[:, None, :] != s2_matrix[None, :, :], axis=2) # Calculate the Hamming distance

def get_closest_gene_vectorized(sgID_array, query_array):
    """Find the closest gene in the array with the smallest Hamming distance to the query sequences."""
    distances = vectorized_hamming_distance(query_array['sgID'], sgID_array['sgRNA_seq'])
    min_indices = np.argmin(distances, axis=1) # Get the index of the minimum distance
    min_distances = distances[np.arange(len(query_array)), min_indices] # Get the minimum distance

    # Only keep the genes where Hamming Distance is <= 1
    closest_keys = np.where(min_distances <= 1, sgID_array['sgID'][min_indices], query_array['sgID']) # Similar to R ifelse
    min_distances = np.where(min_distances <= 1, min_distances, 'NA')
    
    return closest_keys, min_distances

# For batch processing multiple strings against one reference:
def myDistGB_python_batch(str_array: list[str], barcode: str) -> np.ndarray:
    """Vectorized batch implementation for multiple strings"""
    # Convert strings to 2D array
    arr1 = np.array([list(s) for s in str_array])
    arr2 = np.array(list(barcode))
    # Use broadcasting for comparison
    return np.sum(arr1 != arr2, axis=1) > 2

def split_key(key):
    sgID, barcode = key.split(',')
    return sgID, barcode

def finditer_with_line_numbers(pattern, string, flags=0):
    """
    A version of ``re.finditer`` that returns ``(match, line_number)`` pairs.
    """
    matches = list(re.finditer(pattern, string, flags))
    if matches:
        end = matches[-1].start()
        # -1 so a failed `rfind` maps to the first line.
        newline_table = {-1: 0}
        for i, m in enumerate(re.finditer("\\n", string), 1):
            # Don't find newlines past our last match.
            offset = m.start()
            if offset > end:
                break
            newline_table[offset] = i

        # Failing to find the newline is OK, -1 maps to 0.
        for m in matches:
            newline_offset = string.rfind("\n", 0, m.start())
            line_number = newline_table[newline_offset]
            yield (m, line_number)

# Add helper function
def extract_groups(sequence, pattern):
    match = pattern.search(sequence)
    if match:
        return match.groups()
    return [None] * 3  # Return None for each expected group

def get_regex_patterns(min_length, max_length, all_start_with_G):
    # regex = re.compile('GTT' + (17N random BC) + 'ATGG' + (sgRNA/sgID) + 'GTTTAA') Moba 500
    if min_length > 8:
        # Use regex for guide RNAs longer than 8 if all sequences start with G
        if all_start_with_G:
            regexR1 = re.compile(f'.*GTT(.{{17}})ATG(.{{{min_length},{max_length}}})GTT(TAA){{e<2}}')  # 15 nt for alignment purpose
            #regexR2 = re.compile(f'.*GTT(.{{17}})ATG(.{{{min_length},{max_length}}})GTT(TAA){{e<3}}')  # Looser criteria for R2
            regexR2 = re.compile(f'(TTA){{e<3}}AAC(.{{{min_length},{max_length}}})CAT(.{{17}})AAC.*')  # Reverse complement of R2
        # Regex for when gRNA sequences do not all start with G
        else:
            regexR1 = re.compile(f'.*GTT(.{{17}})ATGG(.{{{min_length},{max_length}}})GTT(TAA){{e<2}}')  # 15 nt for alignment purpose
            #regexR1 = re.compile(f'.*GTT(.{{17}})ATGG(.{{{min_length},{max_length}}})GTT(TAA)')  # 15 nt for alignment purpose without fuzzy matching
            #regexR2 = re.compile(f'.*GTT(.{{17}})ATGG(.{{{min_length},{max_length}}})GTT(TAA){{e<3}}')  # Looser criteria for R2
            regexR2 = re.compile(f'(TTA){{e<3}}AAC(.{{{min_length},{max_length}}})CCAT(.{{17}})AAC.*', flags=re.REVERSE)  # Reverse complement of R2
            #regexR2 = re.compile(f'(TTA)AAC(.{{{min_length},{max_length}}})CCAT(.{{17}})AAC.*')  # Reverse complement of R2 without fuzzy matching
    else:
        # Original 3D Moba regex = re.compile(‘GA’ + (sgID) + (random_barcode) + ‘ATGCCCAAGAAG’)\n”,
        regexR1 = re.compile(r'GA(.{8})(GC.{5}TA.{5}GC.{5}TA.{5}GC)(ATGCCCA){e<2}')  # Allowing 2 errors in the ATGCCCA
        #regexR2 = re.compile(r'A(.{8})(GC.{5}TA.{5}GC.{5}TA.{5}GC)A(TGCCCA){e<3}')  # Allowing 3 errors in the TGCCCA
        regexR2 = re.compile(r'(TGGGCA){e<3}T(GC.{5}TA.{5}GC.{5}TA.{5}GC)(.{8})T')  # Reverse complement of R2
        # The GC.{5}TA.{5}GC.{5}TA.{5}GC is the random barcode that is 17 nt long has specific GC and TA sequences separated by 5 Ns
    return regexR1, regexR2

def countreads(read_one, read_two, sample_name, sgIDs, min_length, max_length, all_start_with_G, out_dir, debug):
    try:
        with timeout_handler(300, f"count-reads({sample_name})"):
            log.logit(f"Starting to get the read counts from: {sample_name}.", color="green")
            log.logit(f"Opening {os.path.basename(read_one)} and {os.path.basename(read_two)}...")
            
            # This returns a list with the sequences from the FASTQ file
            read_one_sequences = tools.read_sequences_from_fastq(read_one, debug)
            #f2 = [reverse_complement(seq) for seq in read_sequences_from_fastq(read_two)]
            read_two_sequences = tools.read_sequences_from_fastq(read_two, debug)    
            log.logit(f"Finished reading sequences from {os.path.basename(read_one)} and {os.path.basename(read_two)}.")
            #read_two_sequences = vectorized_reverse_complement(read_two_sequences)
            #if debug: log.logit(f"Finished performing reverse complement on {os.path.basename(read_two)}.", color="yellow")
            if debug: 
                log.logit(f"There was a total of {len(read_one_sequences)} sequences in {os.path.basename(read_one)}.", color="yellow")
                log.logit(f"There was a total of {len(read_two_sequences)} sequences in {os.path.basename(read_two)}.", color="yellow")
            if len(read_one_sequences) != len(read_two_sequences):
                err_msg = f"ERROR: The number of sequences in {os.path.basename(read_one)} and {os.path.basename(read_two)} do not match. Exiting."
                log.logit(err_msg, color="red")
                sys.exit(f"[err] {err_msg}")
            
            regexR1, regexR2 = get_regex_patterns(min_length, max_length, all_start_with_G)

            sgID_barcodes = {}
            log.logit(f"Starting to match sequences from {os.path.basename(read_one)} and {os.path.basename(read_two)}...")
            if debug: log.logit(f"Using the following regex patterns: {regexR1} and {regexR2}.")

            for seq1, seq2 in tqdm(zip(read_one_sequences, read_two_sequences), total=len(read_one_sequences)):  # Add progress bar
            #for seq1, seq2 in zip(read_one_sequences, read_two_sequences):
                # Check if the sequences match the regex pattern
                r1_match = regexR1.search(seq1)
                #r2_match = regexR2.search(tools.reverse_complement(seq2))
                r2_match = regexR2.search(seq2)
                if r1_match and r2_match:
                    # If we can see patterns in both reads
                    if min_length > 8:
                        r1_sgID = r1_match.group(2)  # For longer guide RNAs, R1sgID is from group 2
                        r1_barcode = r1_match.group(1)    # R1BC is from group 1
                        r2_sgID = tools.reverse_complement(r2_match.group(2))  # For longer guide RNAs, R2sgID is from group 2
                        r2_barcode = tools.reverse_complement(r2_match.group(3))    # R2BC is from group 1
                    else: #TODO: If we do change to not do the reverse complement, need to account for this
                        r1_sgID = r1_match.group(1)  # For shorter guide RNAs, R1sgID is from group 1
                        r1_barcode = r1_match.group(2)    # R1BC is from group 2
                        r2_sgID = r2_match.group(1)  # For shorter guide RNAs, R2sgID is from group 1
                        r2_barcode = r2_match.group(2)    # R2BC is from group 2

                    if r1_barcode == r2_barcode and r1_sgID == r2_sgID and (('N' not in r1_barcode) or ('N' not in r2_barcode)):
                    #if r1_barcode == r2_barcode and ('N' not in r1_barcode):   # This is wrong
                        dict_key = (f"{r1_sgID},{r1_barcode}")
                        # If sgID is not in the dictionary, add it
                        if dict_key not in sgID_barcodes:
                            sgID_barcodes[dict_key] = 1
                        else:
                            sgID_barcodes[dict_key] += 1
            if debug: log.logit(f"Found {len(sgID_barcodes)} unique sgIDs.", color="yellow")
            log.logit(f"Finished matching sequences from {os.path.basename(read_one)} and {os.path.basename(read_two)}.")
            # Convert the sgID_barcodes dictionary to a numpy array of tuples
            sgID_barcodes_list = list(sgID_barcodes.items())
            sgID_barcodes_split = [(split_key(key)[0], split_key(key)[1], value) for key, value in sgID_barcodes_list]
            sgID_barcodes_array = np.array(sgID_barcodes_split, dtype=[('sgID', 'U50'), ('barcode', 'U50'), ('total', 'i4')])
            # Convert the sgIDs dictionary to a numpy array of tuples
            sgID_array = np.array(list(sgIDs.items()), dtype=[('sgID', 'U50'), ('sgRNA_seq', 'U50')])
            # Get the closest gene for each sgID in the sgID_barcodes_array
            log.logit(f"Matching sgIDs to genes...")
            gene, distance = get_closest_gene_vectorized(sgID_array, sgID_barcodes_array)
            log.logit(f"Finished matching sgIDs to genes.")

            log.logit(f"Creating output dictionary...")
            combined_dict = {}
            for gene, distance, barcode, total in zip(gene, distance, sgID_barcodes_array['barcode'], sgID_barcodes_array['total']):
                key = f"{gene},{distance},{barcode}"
                if key not in combined_dict:
                    combined_dict[key] = total
                else:
                    combined_dict[key] += total
            log.logit(f"Finished matching sequences from {os.path.basename(read_one)} and {os.path.basename(read_two)}.")
            log.logit(f"Writing output to {sample_name}_MergeReadOut.csv.")
            with open(out_dir + "/" + sample_name + "_MergeReadOut.csv", "w") as f:
                for key in sorted(combined_dict, key=combined_dict.get, reverse=True):
                    f.write(",".join(map(str, [key, combined_dict[key]])) + "\n")
            log.logit(f"Finished writing output to {sample_name}_MergeReadOut.csv.")
            return True
    except TimeoutError as te:
        log.logit(f"{sample_name} - TimeoutError: {te}", color="red")
        return False
    except Exception as e:
        log.logit(f"Error occurred while processing {sample_name}: {str(e)}", color="red")
        if debug: log.logit(traceback.format_exc())
        return False

def remove_spurious_barcodes(sgid_df, sample_name, out_dir, debug):
    with indent(4, quote=' >'): 
        if debug: log.logit(f"sgID: {sgid_df['sgID'].iloc[0]}")
    sgid_df = sgid_df.sort_values(by='count', ascending=False)
    if len(sgid_df) <= 1:
        return sgid_df
    
    sgid_df = sgid_df.reset_index(drop=True)
    # Convert barcodes to numpy character array
    barcodes = np.array([list(bc) for bc in sgid_df['barcode']])
    # Calculate distances using broadcasting
    distances = barcodes[:, None, :] != barcodes[None, :, :]            # Shape will be (n_barcodes, n_barcodes, barcode_length)
    hamming_distances = np.sum(distances, axis=2)                       # Sum differences along last axis to get Hamming distances
    mask = np.triu(np.ones_like(hamming_distances, dtype=bool), k=1)    # Create mask for upper triangle (we only need to compare each pair once)
    similar_pairs = np.where((hamming_distances <= 2) & mask)           # Get indices where distance <= 2 (similar barcodes)
    # Initialize result array
    bool_all = np.ones(len(sgid_df), dtype=bool)
    # Create hamming_dstances DataFrame
    hamming_distances_df = pd.DataFrame(hamming_distances, index=sgid_df.index, columns=sgid_df.index)
    data = []
    # Mark false using original indices
    for i, j in zip(*similar_pairs):
        # We skip barcodes with less than 5 counts because we don't know which one is the main colony
        if sgid_df.iloc[i]['count'] < 5:
            continue
        already_removed = not bool_all[i] or not bool_all[j]
        if not already_removed:
            data.append({
                'sgID': sgid_df['sgID'].iloc[0],
                'barcode1': sgid_df.iloc[i]['barcode'],
                'barcode2': sgid_df.iloc[j]['barcode'],
                'count1': sgid_df.iloc[i]['count'],
                'count2': sgid_df.iloc[j]['count'],
                'hamming_distance': hamming_distances_df.iloc[i][j],
                'reads_rel_to_parent' : sgid_df.iloc[j]['count'] / sgid_df.iloc[i]['count']
            })
        # However, for all other barcodes, we remove the one with the lower count
        # and since this is organized by count_size, [i] SHOULD always be greater than [j]
        # exception is when i == j, which is the same barcode, then it's a 50/50 toss up...
        bool_all[j] = False
    barcode_df = pd.DataFrame(data) if data else pd.DataFrame()
    if debug and len(barcode_df) > 0:
        output_file = out_dir + "/" + sample_name + "_BarcodeHandshakeResults.txt"
        write_header = not os.path.exists(output_file)
        barcode_df.to_csv(
            output_file,
            mode='a',
            header=write_header,
            index=False,
            sep='\t',
            na_rep='NA'
        )
    # Return a new DataFrame with the filtered barcodes
    return sgid_df[bool_all]
    
def clean_barcodes(merge_reads_csv, sample_name, threads, out_dir, debug):
    try:
        with timeout_handler(300, f"clean-barcodes({sample_name})"):
            log.logit(f"Starting to clean barcodes from: {sample_name}.", color="green")
            log.logit(f"Opening {os.path.basename(merge_reads_csv)}...")
            
            df = pd.read_csv(merge_reads_csv, header=None, names=["sgID", "distance", "barcode", "count"])

            log.logit(f"Finished reading {os.path.basename(merge_reads_csv)}.")
            log.logit(f"Starting to clean barcodes from {os.path.basename(merge_reads_csv)}...")
            # Remove potentially old Handshake results
            if os.path.exists(out_dir + "/" + sample_name + "_BarcodeHandshakeResults.txt") and debug:
                os.remove(out_dir + "/" + sample_name + "_BarcodeHandshakeResults.txt")
            # Remove barcodes where the count is less than 2
            df = df[df['count'] >= 2]
            if threads == 1:
                clean_df = df.groupby('sgID').apply(lambda x: remove_spurious_barcodes(x, sample_name, out_dir, debug))
            else:
                with mp.Pool(threads) as p:
                    results = p.starmap(remove_spurious_barcodes, [(group, sample_name, out_dir, debug) for _, group in df.groupby('sgID')])
                clean_df = pd.concat(results, ignore_index=True)
            log.logit(f"Finished cleaning barcodes from {os.path.basename(merge_reads_csv)}.")
            log.logit(f"Writing output to {out_dir}/{sample_name}_BarcodeClean.txt.")
            clean_df.to_csv(out_dir + "/" + sample_name + "_BarcodeClean.txt", index=False, sep='\t', na_rep = 'NA')
            return True
    except TimeoutError as te:
        log.logit(f"{sample_name} - TimeoutError: {te}", color="red")
        return False
    except Exception as e:
        log.logit(f"Error occurred while processing {sample_name}: {str(e)}", color="red")
        if debug: log.logit(traceback.format_exc())
        return False
    
def cell_num_metrics(df, spike_ins):
    df.loc[(df['sgID'] == 'sgDummy') & (df['barcode'].isin(spike_ins['sequence'])), 'sgID'] += '_spike'
    df['BC_end'] = df['barcode'].str[-4:]
    df['fullBC'] = df['sgID'] + '_' + df['BC_end']
    df['log2_cell_num'] = np.log2(df['cell_num'])
    df['rel_cell_num'] = df.groupby('Sample_ID')['cell_num'].transform(lambda x: x / x.sum())
    return df
    
def get_slope_and_rsq(spike_count):
    log.logit(f"Calculating the slope and R² value for the spike-ins...")
    # Prepare data
    X = spike_count[['count']].values
    y = spike_count['expected_cellnum'].values

    # Fit regression without intercept
    reg = LinearRegression(fit_intercept=False).fit(X, y)

    # Get coefficient (slope)
    slope = reg.coef_[0]

    # Calculate R²
    y_pred = reg.predict(X)
    # https://github.com/scikit-learn/scikit-learn/discussions/21050
    # https://stats.stackexchange.com/questions/26176/removal-of-statistically-significant-intercept-term-increases-r2-in-linear-mo/26205#26205
    # 1 - (((y - y_pred) ** 2).sum() / ((y - np.mean(y)) ** 2).sum())
    rsq = 1 - (np.sum((y - y_pred) ** 2) / np.sum((y - np.mean(y)) ** 2))
    # 1 - (((y - y_pred) ** 2).sum() / (y ** 2).sum())
    rsq = 1 - (np.sum((y - y_pred) ** 2) / np.sum(y ** 2))

    # NOTE: My opinion is to NOT use the correlation of coefficient of determination (R²) as the metric for goodness of fit
    # This is because we assume no intercept, which means that the R² value will be deflated as it is not the best fit, but the best fit with no intercept

    # Calculate residuals
    residuals = y - y_pred
    
    # Calculate z-scores
    z_scores = (residuals - np.mean(residuals)) / np.std(residuals)
    
    # Add to DataFrame
    spike_count = spike_count.copy()
    # Use .loc to set values
    spike_count.loc[:, 'residuals'] = residuals
    spike_count.loc[:, 'z_scores'] = z_scores

    return spike_count, slope, rsq

def cell_number(barcode_clean_txt, sample_name, spike_ins, library_info, out_dir, plot, debug):
    try:
        with timeout_handler(300, f"cell-number({sample_name})"):
            log.logit(f"Starting to get cell numbers from: {sample_name}.", color="green")
            log.logit(f"Opening {os.path.basename(barcode_clean_txt)}...")
            df = pd.read_csv(barcode_clean_txt, sep='\t')
            spike_data = df[df['sgID'] == 'sgDummy']
            log.logit(f"Finished reading {os.path.basename(barcode_clean_txt)}.")
            log.logit(f"Starting to calculate cell numbers from {os.path.basename(barcode_clean_txt)}...")
            # Getting the spike-in information by merging the spike_ins and barcode_clean_txt DataFrames 
            log.logit(f"Determining the raw counts of spike-ins...")
            spike_count = pd.merge(
                spike_data,
                spike_ins[['sequence', 'expected', 'name']],
                left_on='barcode',
                right_on='sequence'
            )
            if len(spike_count) == 0:
                log.logit(f"WARNING: No spike-ins found for sample {sample_name}", color = "yellow")
            spike_num = library_info[library_info['Sample'] == sample_name]['Spike_Num'].iloc[0]
            log.logit(f"Calculating the expected cell numbers based on the spike-ins having a starting value of {spike_num}...")
            spike_count = spike_count.assign(
                Sample_ID=sample_name,
                expected_cellnum=lambda x: x['expected'] * spike_num
            )
            if len(spike_count) > 1: # If there is more than one spike-in, calculate the slope and R² value
                spike_count, slope, rsq = get_slope_and_rsq(spike_count)
                if debug: 
                    log.logit(f"SciKit-Learn Linear Regression Model: {slope}, {rsq}")
                    print(spike_count[['Sample_ID', 'count', 'name', 'expected_cellnum', 'z_scores']])
                if plot: plotting.spike_ins(spike_count, slope, rsq, "red" if rsq < 0.9 else "black", sample_name, out_dir)
                # If the R² value is less than 0.9, remove the Spike-in with the highest z-score
                if rsq < 0.9:
                    log.logit(f"The R² value is {rsq}, trying to remove the Spike-in with the highest z-score...")
                    max_abs_zscore = np.abs(spike_count['z_scores']).max()
                    spike_count_redo = spike_count[np.abs(spike_count['z_scores']) != max_abs_zscore]
                    spike_count_redo, slope_redo, rsq_redo = get_slope_and_rsq(spike_count_redo)
                    if debug: 
                        log.logit(f"SciKit-Learn Linear Regression Model: {slope_redo}, {rsq_redo}")
                        print(spike_count_redo[['Sample_ID', 'count', 'name', 'expected_cellnum', 'z_scores']])
                    if plot: plotting.spike_ins(spike_count_redo, slope_redo, rsq_redo, "red" if rsq < 0.9 else "black", sample_name, out_dir, True)
                    if rsq_redo > rsq:
                        spike_count = spike_count_redo
                        slope = slope_redo
                        rsq = rsq_redo
                spike_count.to_csv(out_dir + "/" + sample_name + "_SpikeInCounts.txt", index=False, sep='\t', na_rep = 'NA')
                log.logit(f"The highest R² value is {rsq}.")
            elif len(spike_count) == 1: # Calculate slope directly from the ratio when only one point is available
                slope = float(spike_count['expected_cellnum'].values / spike_count['count'].values)
            if len(spike_count) != 0:
                log.logit(f"Using the slope of {slope} to calculate the cell numbers.")
                df = df.assign(
                    cell_num=lambda x: x['count'] * slope,
                    reading_depth=1/slope
                )
            else:
                log.logit(f"WARNING: No spike-ins found for sample {sample_name}. Cannot calculate cell numers.", color = "yellow")
                df = df.assign(
                    cell_num=np.nan,
                    reading_depth=np.nan
                )
            log.logit(f"Finished calculating cell numbers from {os.path.basename(barcode_clean_txt)}.")
            # Add additional information to the DataFrame
            df = df.assign(
                Sample_ID = sample_name,
                Tissue = library_info[library_info['Sample'] == sample_name]['Tissue'].iloc[0],
                Time_Point = library_info[library_info['Sample'] == sample_name]['Time_Point'].iloc[0],
                Mouse_ID = library_info[library_info['Sample'] == sample_name]['Mouse_ID'].iloc[0],
                Mouse_Genotype = library_info[library_info['Sample'] == sample_name]['Mouse_Genotype'].iloc[0],
                Sex = library_info[library_info['Sample'] == sample_name]['Sex'].iloc[0],
                Cell_Lib = library_info[library_info['Sample'] == sample_name]['Cell_Library'].iloc[0],
                Body_Weight_Trans = library_info[library_info['Sample'] == sample_name]['BW_Trans'].iloc[0],
                Body_Weight_Coll = library_info[library_info['Sample'] == sample_name]['BW_Coll'].iloc[0],
                Tissue_Weight = library_info[library_info['Sample'] == sample_name]['Tissue_Weight'].iloc[0],
                Spike_Num = spike_num,
                Spike_Rsq = rsq if len(spike_count) > 1 else 0
            )
            df['Sample_Rsq'] = df['Sample_ID'] + "_" + df['Spike_Rsq'].astype(str)
            df = cell_num_metrics(df, spike_ins)
            df.to_csv(out_dir + "/" + sample_name + "_UnfilteredSampleInfo.txt", index=False, sep='\t', na_rep = 'NA')
            filtered_df = df[~df['sgID'].str.contains('_spike') & (df['cell_num'] > 1)]
            filtered_df = filtered_df[~filtered_df['distance'].isna()]
            filtered_df.to_csv(out_dir + "/" + sample_name + "_FilteredSampleInfo.txt", index=False, sep='\t', na_rep = 'NA')
            spike_count.to_csv(out_dir + "/" + sample_name + "_SpikeInInfo.txt", index=False, sep='\t', na_rep = 'NA')
            return True, df, filtered_df, spike_count
    except TimeoutError as te:
        log.logit(f"{sample_name} - TimeoutError: {te}", color="red")
        return False, None, None, None
    except Exception as e:
        log.logit(f"Error occurred while processing {sample_name}: {str(e)}", color="red")
        if debug: log.logit(traceback.format_exc())
        return False, None, None, None