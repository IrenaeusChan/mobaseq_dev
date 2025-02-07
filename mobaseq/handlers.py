import os, sys
import multiprocessing as mp
import pandas as pd
import numpy as np
import mobaseq.process.tools as tools
import mobaseq.process.rawdata as rawdata
import mobaseq.qc.summarize as summarize
import mobaseq.utils.logger as log

def countreads_single(read_one, read_two, sample_name, sgid_file, out_dir, debug):
    sample_name = tools.process_fastq_files(read_one, read_two, debug)
    sgID_dict = tools.process_sgid_file(sgid_file, debug)
    min_length, max_length, all_start_with_G = tools.get_info_about_sgID_file(sgID_dict, debug)
    rawdata.countreads(read_one, read_two, sample_name, sgID_dict, min_length, max_length, all_start_with_G, str(out_dir), debug)

def countreads_batch(input_dir, sgid_file, out_dir, threads, debug):
    list_of_fastqs = tools.get_list_of_files(input_dir, "fastq", debug)
    sgID_dict = tools.process_sgid_file(sgid_file, debug)
    min_length, max_length, all_start_with_G = tools.get_info_about_sgID_file(sgID_dict, debug)

    with mp.Pool(threads) as p:
        success = p.starmap(rawdata.countreads, [
            (fastq1, fastq2, sample_name, sgID_dict, min_length, max_length, all_start_with_G, str(out_dir), debug)
            for fastq1, fastq2, sample_name in list_of_fastqs
        ])
    if not all(success):
        err_msg = f"ERROR: Some processes failed. Check the logs for more details. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")

def clean_barcodes_single(merge_reads_csv, threads, out_dir, debug):
    sample_name = tools.get_sample_name(merge_reads_csv, "merge_read_csv", debug)
    rawdata.clean_barcodes(merge_reads_csv, sample_name, threads, out_dir, debug)

def clean_barcodes_batch(input_dir, out_dir, threads, debug):
    merge_reads_csv_files = tools.get_list_of_files(input_dir, "merge_reads_csv", debug)
    # Each clean_barcodes is also multi-threaded, so if we have 4 threads and 4 samples, we will have 16 threads running
    # therefore, to avoid overloading the system, we will use half the threads for each of the samples
    with mp.Pool(threads) as p:
        success = p.starmap(rawdata.clean_barcodes, [
            (merge_read_csv[0], merge_read_csv[1], 1, out_dir, debug) for merge_read_csv in merge_reads_csv_files
        ])
    if not all(success):
        err_msg = f"ERROR: Some processes failed. Check the logs for more details. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")

def sgID_qc_single(merge_reads_csv, sgid_file, out_dir, debug):
    sgID_dict = tools.process_sgid_file(sgid_file, debug)
    sample_name = tools.get_sample_name(merge_reads_csv, "merge_read_csv", debug)
    summarize.sgID_info(merge_reads_csv, sample_name, sgID_dict, str(out_dir), debug)

def sgID_qc_batch(input_dir, sgid_file, out_dir, threads, debug):
    sgID_dict = tools.process_sgid_file(sgid_file, debug)
    merge_reads_csv_files = tools.get_list_of_files(input_dir, "merge_reads_csv", debug)

    with mp.Pool(threads) as p:
        res = p.starmap(summarize.sgID_info, [
            (merge_reads_csv, sample_name, sgID_dict, str(out_dir), debug)
            for merge_reads_csv, sample_name in merge_reads_csv_files
        ])

    success = [res[0] for res in res]
    results = [[res[1], res[2]] for res in res]

    if not all(success):
        err_msg = f"ERROR: Some processes failed. Check the logs for more details. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    else:
        # Initialize empty DataFrames for combined results using the first result
        combined_total_reads = results[0][0]['sgID']
        combined_num_barcodes = results[0][0]['sgID']
        combined_rel_reads = results[0][0]['sgID']
        combined_rel_barcodes = results[0][0]['sgID']

        # Merge each result DataFrame into the combined DataFrames
        for i, res in enumerate(results):
            df, sample_name = res[0], res[1]
            combined_total_reads = pd.merge(combined_total_reads, df[['sgID', 'total_reads']].rename(columns={'total_reads': f'{sample_name}'}), on='sgID', how='outer').fillna(0)
            combined_num_barcodes = pd.merge(combined_num_barcodes, df[['sgID', 'num_barcodes']].rename(columns={'num_barcodes': f'{sample_name}'}), on='sgID', how='outer').fillna(0)
            combined_rel_reads = pd.merge(combined_rel_reads, df[['sgID', 'rel_reads']].rename(columns={'rel_reads': f'{sample_name}'}), on='sgID', how='outer').fillna(0)
            combined_rel_barcodes = pd.merge(combined_rel_barcodes, df[['sgID', 'rel_barcodes']].rename(columns={'rel_barcodes': f'{sample_name}'}), on='sgID', how='outer').fillna(0)

        # Save each combined DataFrame to a separate CSV file
        combined_total_reads.to_csv(os.path.join(out_dir, "total_reads_per_sgid.csv"), index=False)
        combined_num_barcodes.to_csv(os.path.join(out_dir, "barcodes_per_sgid.csv"), index=False)
        combined_rel_reads.to_csv(os.path.join(out_dir, "relative_reads_per_sgid.csv"), index=False)
        combined_rel_barcodes.to_csv(os.path.join(out_dir, "relative_barcodes_per_sgid.csv"), index=False)

def mapped_reads_single(merge_reads_csv, sgid_file, out_dir, debug):
    log.logit(f"Summarizing mapped reads for {merge_reads_csv}...")
    sample_name = tools.get_sample_name(merge_reads_csv, "merge_read_csv", debug)
    sgID_dict = tools.process_sgid_file(sgid_file, debug)
    total_reads, unmapped_reads = summarize.get_total_reads(merge_reads_csv, sgID_dict)
    mapped_percentage = [((total_reads - unmapped_reads) / total_reads) * 100]
    unmapped_percentage = [(unmapped_reads / total_reads) * 100]
    log.logit(f"Writing mapped percentages to csv...")
    df = pd.DataFrame(
        {
            "Sample": sample_name,
            "Mapped %": mapped_percentage,
            "Unmapped %": unmapped_percentage,
        }
    )
    df.to_csv(f"{out_dir}/{sample_name}_mapped_percentages.csv", index=False)
    log.logit(f"Results written to {out_dir}/{sample_name}_mapped_percentages.csv")

def mapped_reads_batch(input_dir, sgid_file, out_dir, threads, debug):
    log.logit(f"Summarizing mapped reads...")
    merge_reads_csv_files = tools.get_list_of_files(input_dir, "merge_reads_csv", debug)
    sgID_dict = tools.process_sgid_file(sgid_file, debug)

    with mp.Pool(threads) as p:
        results = p.starmap(summarize.get_total_reads, [
            (merge_read_csv[0], sgID_dict) for merge_read_csv in merge_reads_csv_files
        ])

    sample_names = [merge_read_csv[1] for merge_read_csv in merge_reads_csv_files]
    total_reads = np.array([res[0] for res in results])
    unmapped_reads = np.array([res[1] for res in results])
    total_reads_safe = np.where(total_reads == 0, 1, total_reads)     # Avoid division by zero by setting total_reads to 1 where it is 0
    mapped_percentages = ((total_reads - unmapped_reads) / total_reads_safe) * 100
    unmapped_percentages = (unmapped_reads / total_reads_safe) * 100
    # Set percentages to 0 where total_reads was originally 0
    mapped_percentages[total_reads == 0] = 0
    unmapped_percentages[total_reads == 0] = 0

    mapped_percentages = mapped_percentages.tolist()
    unmapped_percentages = unmapped_percentages.tolist()

    log.logit(f"Writing mapped percentages to csv...")
    df = pd.DataFrame(
        {
            "Sample": sample_names,
            "Mapped %": mapped_percentages,
            "Unmapped %": unmapped_percentages,
        }
    )
    df.to_csv(f"{out_dir}/mapped_percentages.csv", index=False)
    log.logit(f"Results written to {out_dir}/mapped_percentages.csv")
    return df, mapped_percentages, unmapped_percentages