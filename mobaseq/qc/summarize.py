import os, sys
import pandas as pd
#import numpy as np
#from tqdm import tqdm
import traceback
import mobaseq.utils.logger as log
import mobaseq.process.tools as tools
import multiprocessing as mp
import numpy as np

def get_total_reads(csv_file, sgIDs):
    df = pd.read_csv(csv_file, sep=",", header=None, names=['sgID', 'hamming_dist', 'barcode', 'num_reads'])
    total_reads = df['num_reads'].sum() # Total number of reads for this sample
    unmapped_reads = df.loc[~df['sgID'].isin(sgIDs.keys()), 'num_reads'].sum() # Total number of reads that did not map to any sgID
    return total_reads, unmapped_reads

'''Function to output a summary of the sgRNA counts for each sample'''
def sgID_info(merge_read_csv, sample_name, sgIDs, out_dir, debug):
    # We want 4 pieces of information:
    # 1. Total reads per sgID (sum of reads for all barcodes)
    # 2. Number of barcodes per sgID (number of rows for each sgID)
    # 3. Relative reads per sgID (total reads per sgID divided by total reads)
    # 4. Relative barcodes per sgID (number of barcodes per sgID divided by total number of barcodes)
    try:
        log.logit(f"Summarizing sgID information from '{merge_read_csv}'...")
        df = pd.read_csv(merge_read_csv, sep=",", header=None, names=['sgID', 'hamming_dist', 'barcode', 'num_reads'])
        df = df[df['hamming_dist'] == 0] # Filter out rows where hamming_dist != 0

        # Initialize a dictionary to store the results
        res = {
            'sgID': [],
            'total_reads': [],
            'num_barcodes': [],
            'rel_reads': [],
            'rel_barcodes': []
        }

        total_reads, unmapped_reads = get_total_reads(merge_read_csv, sgIDs)
        total_barcodes = len(df) # Total number of barcodes for this sample

        # For each sgID (gene), there will be a group of barcodes
        for sgID, group in df.groupby('sgID'):
            total_reads_per_sgID = group['num_reads'].sum()
            num_barcodes_per_sgID = len(group)
            rel_reads_per_sgID = total_reads_per_sgID / total_reads
            rel_barcodes_per_sgID = num_barcodes_per_sgID / total_barcodes

            res['sgID'].append(sgID)
            res['total_reads'].append(total_reads_per_sgID)
            res['num_barcodes'].append(num_barcodes_per_sgID)
            res['rel_reads'].append(rel_reads_per_sgID)
            res['rel_barcodes'].append(rel_barcodes_per_sgID)

        res_df = pd.DataFrame(res)
        res_df.to_csv(os.path.join(out_dir, f'{sample_name}_sgID_info.csv'), index=False)
        return True, res_df, sample_name
    except Exception as e:
        log.logit(f"Error reading '{merge_read_csv}': {e}")
        if debug:
            traceback.print_exc()
        return False, None, None
