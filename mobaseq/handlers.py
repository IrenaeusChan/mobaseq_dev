import os, sys
import multiprocessing as mp
import pandas as pd
import numpy as np
import traceback
import mobaseq.process.tools as tools
import mobaseq.process.rawdata as rawdata
import mobaseq.qc.summarize as summarize
import mobaseq.qc.plot as plotting
import mobaseq.utils.logger as log

def process_with_timeout(pool, func, args_list, timeout=3600):
    results = []
    for args in args_list:
        try:
            # Start async process
            async_result = pool.apply_async(func, args)
            # Wait for result with timeout
            result = async_result.get(timeout=timeout)
            results.append(result)
        except mp.TimeoutError:
            log.logit(f"Process timed out after {timeout}s: {args[1]}", color="red")
            results.append(False)
        except Exception as e:
            log.logit(f"Process failed: {args[1]} - {str(e)}", color="red")
            results.append(False)
    return results

def parallel_process_with_timeout(pool, func, args_list, timeout=1800):
    # Start all processes asynchronously
    async_results = [
        pool.apply_async(func, args) 
        for args in args_list
    ]
    
    # Collect results with timeout
    results = []
    for i, async_result in enumerate(async_results):
        try:
            result = async_result.get(timeout=timeout)
            results.append(result)
        except mp.TimeoutError:
            log.logit(f"Process {i} timed out after {timeout}s", color="red")
            results.append(False)
        except Exception as e:
            log.logit(f"Process {i} failed: {str(e)}", color="red")
            results.append(False)
    
    return results

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
        # success = p.starmap(rawdata.countreads, [
        #     (fastq1, fastq2, sample_name, sgID_dict, min_length, max_length, all_start_with_G, str(out_dir), debug)
        #     for fastq1, fastq2, sample_name in list_of_fastqs
        # ])
        success = parallel_process_with_timeout(
            pool=p,
            func=rawdata.countreads,
            args_list=[
                (fastq1, fastq2, sample_name, sgID_dict, min_length, max_length, all_start_with_G, str(out_dir), debug)
                for fastq1, fastq2, sample_name in list_of_fastqs
            ],
            timeout=1800  # 30 min timeout
        )
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
        # success = p.starmap(rawdata.clean_barcodes, [
        #     (merge_read_csv, sample_name, 1, out_dir, debug) for merge_read_csv, sample_name in merge_reads_csv_files
        # ])
        success = parallel_process_with_timeout(
            pool=p,
            func=rawdata.clean_barcodes,
            args_list=[
                (merge_read_csv, sample_name, 1, out_dir, debug)
                for merge_read_csv, sample_name in merge_reads_csv_files
            ],
            timeout=1800
        )

    if not all(success):
        err_msg = f"ERROR: Some processes failed. Check the logs for more details. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")

def cell_number_single(barcode_clean_txt, spike_ins, library_info, out_dir, plot, debug):
    sample_name = tools.get_sample_name(barcode_clean_txt, "barcode_clean_txt", debug)
    spike_ins = tools.process_spike_ins(spike_ins, debug)
    library_info = pd.read_excel(library_info)
    rawdata.cell_number(barcode_clean_txt, sample_name, spike_ins, library_info, out_dir, plot, debug)

def cell_number_batch(input_dir, spike_ins, library_info, out_dir, threads, plot, debug):
    barcode_clean_txt_files = tools.get_list_of_files(input_dir, "barcode_clean_txt", debug)
    spike_ins = tools.process_spike_ins(spike_ins, debug)
    library_info = pd.read_excel(library_info)
    
    with mp.Pool(threads) as p:
        # results = p.starmap(rawdata.cell_number, [
        #     (barcode_clean_txt, sample_name, spike_ins, library_info, out_dir, plot, debug)
        #     for barcode_clean_txt, sample_name in barcode_clean_txt_files
        # ])
        results = parallel_process_with_timeout(
            pool=p,
            func=rawdata.cell_number,
            args_list=[
                (barcode_clean_txt, sample_name, spike_ins, library_info, out_dir, plot, debug)
                for barcode_clean_txt, sample_name in barcode_clean_txt_files
            ],
            timeout=1800
        )
    success = [res[0] for res in results]
    dfs = [res[1] for res in results]
    filtered_dfs = [res[2] for res in results]
    spike_counts_dfs = [res[3] for res in results]
    if not all(success):
        err_msg = f"ERROR: Some processes failed. Check the logs for more details. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    else:
        log.logit(f"Cell Number Calculations Complete. Combining results...", color="green")
        # Initialize empty DataFrames for combined results using the first result
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_filtered_df = pd.concat(filtered_dfs, ignore_index=True)
        combined_spike_counts = pd.concat(spike_counts_dfs, ignore_index=True)
        
        # Get the unique Sample_ID and Spike_Rsq values from combined_df
        unique_samples = combined_df['Sample_Rsq'].unique()
        split_samples = [x.split('_') for x in unique_samples]
        Rsq = pd.DataFrame(split_samples, columns=['Sample_ID', 'Spike_Rsq'])
        Rsq['Spike_Rsq'] = pd.to_numeric(Rsq['Spike_Rsq'])
        if debug: print(Rsq)
        if plot: plotting.rsq_per_sample(Rsq, out_dir)
        
        # Select columns and get unique combinations
        unique_df = combined_df[['reading_depth', 'Sample_ID']].drop_duplicates().reset_index(drop=True)
        if debug: print(unique_df)
        if plot: plotting.reading_depth_per_sample(unique_df, out_dir)
        
        aggr_df = (combined_df.groupby('Sample_ID')
           .agg({'reading_depth': 'mean', 'Spike_Rsq': 'mean'})
           .rename(columns={'reading_depth': 'mean_ReadingDepth', 
                          'Spike_Rsq': 'mean_Rsq'})
           .reset_index())
        aggr_df['cellsPerRead'] = 1/aggr_df['mean_ReadingDepth']
        if plot: plotting.mean_rsq_distribution(aggr_df, out_dir)
        if plot: plotting.mean_reading_depth_distribution(aggr_df, out_dir)
        
        log.logit(f"Results written to {out_dir}/combined_results/")
        out_dir = tools.ensure_abs_path(out_dir + "/combined_results/")
        Rsq.to_csv(f"{out_dir}/RsqPerSample.csv", index=False)
        unique_df.to_csv(f"{out_dir}/ReadingDepthPerSample.csv", index=False)
        aggr_df.to_csv(f"{out_dir}/CountInfoPerSample.csv", index=False)
        combined_df.to_csv(f"{out_dir}/UnfilteredSampleInfo.csv", index=False)
        combined_filtered_df.to_csv(f"{out_dir}/FilteredSampleInfo.csv", index=False)
        combined_spike_counts.to_csv(f"{out_dir}/SpikeInInfo.csv", index=False)

def sgID_qc_single(merge_reads_csv, sgid_file, out_dir, debug):
    sgID_dict = tools.process_sgid_file(sgid_file, debug)
    sample_name = tools.get_sample_name(merge_reads_csv, "merge_read_csv", debug)
    summarize.sgID_info(merge_reads_csv, sample_name, sgID_dict, str(out_dir), debug)

def sgID_qc_batch(input_dir, sgid_file, out_dir, threads, debug):
    sgID_dict = tools.process_sgid_file(sgid_file, debug)
    merge_reads_csv_files = tools.get_list_of_files(input_dir, "merge_reads_csv", debug)

    with mp.Pool(threads) as p:
        # res = p.starmap(summarize.sgID_info, [
        #     (merge_reads_csv, sample_name, sgID_dict, str(out_dir), debug)
        #     for merge_reads_csv, sample_name in merge_reads_csv_files
        # ])
        results = parallel_process_with_timeout(
            pool=p,
            func=summarize.sgID_info,
            args_list=[
                (merge_reads_csv, sample_name, sgID_dict, str(out_dir), debug)
                for merge_reads_csv, sample_name in merge_reads_csv_files
            ],
            timeout=1800
        )

    success = [res[0] for res in results]
    result = [[res[1], res[2]] for res in results]

    if not all(success):
        err_msg = f"ERROR: Some processes failed. Check the logs for more details. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    else:
        # Initialize empty DataFrames for combined results using the first result
        combined_total_reads = result[0][0]['sgID']
        combined_num_barcodes = result[0][0]['sgID']
        combined_rel_reads = result[0][0]['sgID']
        combined_rel_barcodes = result[0][0]['sgID']

        # Merge each result DataFrame into the combined DataFrames
        for res in result:
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
        # results = p.starmap(summarize.get_total_reads, [
        #     (merge_read_csv[0], sgID_dict) for merge_read_csv in merge_reads_csv_files
        # ])
        results = parallel_process_with_timeout(
            pool=p,
            func=summarize.get_total_reads,
            args_list=[
                (merge_read_csv[0], sgID_dict) for merge_read_csv in merge_reads_csv_files
            ],
            timeout=1800
        )

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

def plot_spike_ins(spike_count, out_dir, debug):
    log.logit(f"Plotting spike-in counts vs. expected cell numbers", color="green")
    out_dir = tools.ensure_abs_path(out_dir)
    sample_name = tools.get_sample_name(spike_count, "spike_count", debug)
    spike_count = pd.read_csv(spike_count, sep="\t")
    spike_count, slope, rsq = rawdata.get_slope_and_rsq(spike_count)
    plotting.spike_ins(spike_count, slope, rsq, "black", sample_name, str(out_dir))

def plot_rsq_per_sample(rsq, out_dir, debug):
    out_dir = tools.ensure_abs_path(out_dir)
    rsq = pd.read_csv(rsq)
    plotting.rsq_per_sample(rsq, str(out_dir))

def plot_read_depth_per_sample(reading_depth, out_dir, debug):
    out_dir = tools.ensure_abs_path(out_dir)
    reading_depth = pd.read_csv(reading_depth)
    plotting.reading_depth_per_sample(reading_depth, str(out_dir))

def plot_mean_read_depth_distribution(aggr_df, out_dir, debug):
    out_dir = tools.ensure_abs_path(out_dir)
    aggr_df = pd.read_csv(aggr_df)
    plotting.mean_reading_depth_distribution(aggr_df, str(out_dir))

def plot_mean_rsq_distribution(aggr_df, out_dir, debug):
    out_dir = tools.ensure_abs_path(out_dir)
    aggr_df = pd.read_csv(aggr_df)
    plotting.mean_rsq_distribution(aggr_df, str(out_dir))

def check_sgid_files(sgID_csv_dir):
    required_files = [
        'total_reads_per_sgid.csv',
        'barcodes_per_sgid.csv',
        'relative_reads_per_sgid.csv',
        'relative_barcodes_per_sgid.csv'
    ]
    
    missing_files = []
    for file in required_files:
        file_path = os.path.join(sgID_csv_dir, file)
        if not os.path.exists(file_path):
            missing_files.append(file)
    
    if missing_files:
        raise FileNotFoundError(
            f"Missing required files in {sgID_csv_dir}: {', '.join(missing_files)}"
        )
    
    log.logit(f"All required sgID files found in {sgID_csv_dir}")
    return True

def process_dataframe(df, rows_to_exclude=None):
    """Helper function to process each DataFrame"""
    df = df.rename(columns={'preTran1': 'pretransplantation'})
    if rows_to_exclude:
        df = df[~df['sgID'].isin(rows_to_exclude)]
    return df

def calculate_relative_values(df, total_per_sample):
    """Calculate relative values for reads/barcodes"""
    sgid_col = df['sgID']
    relative_df = df.iloc[:, 1:].div(total_per_sample[1:], axis=1)
    relative_df.insert(0, 'sgID', sgid_col)
    return relative_df

def plot_per_sgid(sgID_csv_dir, out_dir, no_dummy, debug):
    log.logit(f"Starting to process sgID data...", color = "green")
    out_dir = tools.ensure_abs_path(out_dir)
    try:
        log.logit(f"Looking for the four CSV files in {sgID_csv_dir}...")
        check_sgid_files(sgID_csv_dir)
        log.logit(f"Reading in CSV Files...")
        dataframes = {}
        for file_type in ['total_reads', 'barcodes', 'relative_reads', 'relative_barcodes']:
            filename = f'{file_type}_per_sgid.csv'
            dataframes[file_type] = pd.read_csv(os.path.join(sgID_csv_dir, filename))

        # Process all DataFrames
        rows_to_exclude = ['sgDummy', 'SpikeIn'] if no_dummy else None
        for key in dataframes:
            dataframes[key] = process_dataframe(dataframes[key], rows_to_exclude)
        
        if no_dummy:
            log.logit(f"Recalculating relative reads and barcodes after filtering out sgDummy...")
            # Pre-filter out rows with sgID == 'sgDummy'
            total_reads = dataframes['total_reads'].sum(axis=0)
            total_barcodes = dataframes['barcodes'].sum(axis=0)

            dataframes['relative_reads'] = calculate_relative_values(dataframes['total_reads'], total_reads)
            dataframes['relative_barcodes'] = calculate_relative_values(dataframes['barcodes'], total_barcodes)


        # At the moment, this cut-off is an arbitrary cut-off... need to think of a better method
        log.logit(f"Determining the most significant sgIDs per sample...")
        threshold = 0.1 if no_dummy else 0.02
        sig_entries = []
        for df_type in ['relative_reads', 'relative_barcodes']:
            sig_entries.extend(filter_significant_entries(dataframes[df_type], threshold))

        # Save results
        sig_df = pd.DataFrame(sig_entries)
        suffix = '_noDummy' if no_dummy else ''
        csv_filename = os.path.join(out_dir, f'significant_relative_values_filtered{suffix}.csv')
        log.logit(f"Writing the most significant sgIDs to {csv_filename}")
        sig_df.to_csv(csv_filename, index=False)
    
        # Plot results
        plotting.per_sgid(
            dataframes['total_reads'],
            dataframes['barcodes'],
            dataframes['relative_reads'],
            dataframes['relative_barcodes'],
            out_dir, no_dummy, debug
        )
    except Exception as e:
        log.logit(f"Error processing sgID data: {str(e)}", color="red")
        if debug: log.logit(traceback.format_exc())
        raise

def filter_significant_entries(df, value_threshold):
    significant_entries = []
    # Skip first column by using iloc[:, 1:]
    data_cols = df.iloc[:, 1:]
    
    for sgID in df['sgID']:
        for sample in data_cols.columns:
            value = data_cols.loc[df['sgID'] == sgID, sample].iloc[0]
            if value > value_threshold:
                significant_entries.append({
                    'Sample': sample, 
                    'sgID': sgID, 
                    'Value': value
                })
    return significant_entries

def run_pipeline(input_dir, sgid_file, spike_ins, library_info, plot, out_dir, threads, debug):
    log.logit(f"Reading in files from {input_dir}...")
    # This returns (fastq1, fastq2, sample_name)
    list_of_fastqs = tools.get_list_of_files(input_dir, "fastq", debug)

    # Make a list of sample_names to check for which files might fail
    sample_names = [sample_name for _, _, sample_name in list_of_fastqs]
    
    log.logit(f"Preparing all necessary files...")
    sgID_dict = tools.process_sgid_file(sgid_file, debug)
    min_length, max_length, all_start_with_G = tools.get_info_about_sgID_file(sgID_dict, debug)
    spike_ins = tools.process_spike_ins(spike_ins, debug)
    library_info = pd.read_excel(library_info)

    # Count Reads
    log.logit(f"Counting reads for all samples...", color="green")
    count_reads_out_dir = out_dir + "/MergeReadOutFiles"
    count_reads_out_dir = tools.ensure_abs_path(count_reads_out_dir)
    with mp.Pool(threads) as p:
        # success = p.starmap(rawdata.countreads, [
        #     (fastq1, fastq2, sample_name, sgID_dict, min_length, max_length, all_start_with_G, str(count_reads_out_dir), debug)
        #     for fastq1, fastq2, sample_name in list_of_fastqs
        # ])
        success = parallel_process_with_timeout(
            pool=p,
            func=rawdata.countreads,
            args_list=[
                (fastq1, fastq2, sample_name, sgID_dict, min_length, max_length, all_start_with_G, str(count_reads_out_dir), debug)
                for fastq1, fastq2, sample_name in list_of_fastqs
            ],
            timeout=1800
        )
    if not all(success):
        err_msg = f"ERROR: Some processes during Count-Reads Task failed. Check the logs for more details. Exiting."
        merge_reads_csv = tools.get_list_of_files(count_reads_out_dir, "merge_reads_csv", debug)
        finished_samples = [sample_name for _, sample_name in merge_reads_csv]
        failed_samples = [sample_name for sample_name in sample_names if sample_name not in finished_samples]
        log.logit(f"Failed samples: {failed_samples}", color="red")
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    log.logit(f"Counting reads complete. Moving on to the next step...", color="green")

    # Clean Barcodes
    log.logit(f"Cleaning barcodes for all samples...", color="green")
    merge_reads_csv_files = tools.get_list_of_files(count_reads_out_dir, "merge_reads_csv", debug)
    clean_barcodes_out_dir = out_dir + "/BarcodeCleanOutFiles"
    clean_barcodes_out_dir = tools.ensure_abs_path(clean_barcodes_out_dir)
    with mp.Pool(threads) as p:
        # success = p.starmap(rawdata.clean_barcodes, [
        #     (merge_read_csv, sample_name, 1, str(clean_barcodes_out_dir), debug) for merge_read_csv, sample_name in merge_reads_csv_files
        # ])
        success = parallel_process_with_timeout(
            pool=p,
            func=rawdata.clean_barcodes,
            args_list=[
                (merge_read_csv, sample_name, 1, str(clean_barcodes_out_dir), debug) for merge_read_csv, sample_name in merge_reads_csv_files
            ],
            timeout=1800
        )
    if not all(success):
        err_msg = f"ERROR: Some processes during Clean-Barcodes Task failed. Check the logs for more details. Exiting."
        barcode_clean_txt_files = tools.get_list_of_files(clean_barcodes_out_dir, "barcode_clean_txt", debug)
        finished_samples = [sample_name for _, sample_name in barcode_clean_txt_files]
        failed_samples = [sample_name for sample_name in sample_names if sample_name not in finished_samples]
        log.logit(f"Failed samples: {failed_samples}", color="red")
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    log.logit(f"Cleaning barcodes complete. Moving on to the next step...", color="green")

    # sgID QC
    log.logit(f"Summarizing sgID information for all samples...", color="green")
    sgid_qc_out_dir = out_dir + "/sgIDQCOutFiles"
    sgid_qc_out_dir = tools.ensure_abs_path(sgid_qc_out_dir)
    with mp.Pool(threads) as p:
        # results = p.starmap(summarize.sgID_info, [
        #     (merge_reads_csv, sample_name, sgID_dict, str(sgid_qc_out_dir), debug) for merge_reads_csv, sample_name in merge_reads_csv_files
        # ])
        results = parallel_process_with_timeout(
            pool=p,
            func=summarize.sgID_info,
            args_list=[
                (merge_reads_csv, sample_name, sgID_dict, str(sgid_qc_out_dir), debug) for merge_reads_csv, sample_name in merge_reads_csv_files
            ],
            timeout=1800
        )
    success = [res[0] for res in results]
    sgID_qc_df_sample_name = [[res[1], res[2]] for res in results]
    if not all(success):
        err_msg = f"ERROR: Some processes failed. Check the logs for more details. Exiting."
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    else:
        # Initialize empty DataFrames for combined results using the first result
        combined_total_reads = sgID_qc_df_sample_name[0][0]['sgID']
        combined_num_barcodes = sgID_qc_df_sample_name[0][0]['sgID']
        combined_rel_reads = sgID_qc_df_sample_name[0][0]['sgID']
        combined_rel_barcodes = sgID_qc_df_sample_name[0][0]['sgID']

        # Merge each result DataFrame into the combined DataFrames
        for sgID_qc_df, sample_name in sgID_qc_df_sample_name:
            combined_total_reads = pd.merge(combined_total_reads, sgID_qc_df[['sgID', 'total_reads']].rename(columns={'total_reads': f'{sample_name}'}), on='sgID', how='outer').fillna(0)
            combined_num_barcodes = pd.merge(combined_num_barcodes, sgID_qc_df[['sgID', 'num_barcodes']].rename(columns={'num_barcodes': f'{sample_name}'}), on='sgID', how='outer').fillna(0)
            combined_rel_reads = pd.merge(combined_rel_reads, sgID_qc_df[['sgID', 'rel_reads']].rename(columns={'rel_reads': f'{sample_name}'}), on='sgID', how='outer').fillna(0)
            combined_rel_barcodes = pd.merge(combined_rel_barcodes, sgID_qc_df[['sgID', 'rel_barcodes']].rename(columns={'rel_barcodes': f'{sample_name}'}), on='sgID', how='outer').fillna(0)

        log.logit(f"Results written to {out_dir}")
        # Save each combined DataFrame to a separate CSV file
        combined_total_reads.to_csv(os.path.join(out_dir, "total_reads_per_sgid.csv"), index=False)
        combined_num_barcodes.to_csv(os.path.join(out_dir, "barcodes_per_sgid.csv"), index=False)
        combined_rel_reads.to_csv(os.path.join(out_dir, "relative_reads_per_sgid.csv"), index=False)
        combined_rel_barcodes.to_csv(os.path.join(out_dir, "relative_barcodes_per_sgid.csv"), index=False)
    log.logit(f"Summarizing sgID information complete. Moving on to the next step...", color="green")

    # Mapped Reads
    log.logit(f"Summarizing mapped reads for all samples...", color="green")
    with mp.Pool(threads) as p:
        # results = p.starmap(summarize.get_total_reads, [
        #     (merge_read_csv, sgID_dict) for merge_read_csv, _ in merge_reads_csv_files
        # ])
        results = parallel_process_with_timeout(
            pool=p,
            func=summarize.get_total_reads,
            args_list=[
                (merge_read_csv, sgID_dict) for merge_read_csv, _ in merge_reads_csv_files
            ],
            timeout=1800
        )
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

    log.logit(f"Results written to {out_dir}/mapped_percentages.csv")
    mapped_df = pd.DataFrame(
        {
            "Sample": sample_names,
            "Mapped %": mapped_percentages,
            "Unmapped %": unmapped_percentages,
        }
    )
    mapped_df.to_csv(f"{out_dir}/mapped_percentages.csv", index=False)
    log.logit(f"Summarizing mapped reads complete. Moving on to the next step...", color="green")

    # Cell Counts
    log.logit(f"Calculating cell counts for all samples...", color="green")
    barcode_clean_txt_files = tools.get_list_of_files(clean_barcodes_out_dir, "barcode_clean_txt", debug)
    cell_number_out_dir = out_dir + "/CellCountOutFiles"
    cell_number_out_dir = tools.ensure_abs_path(cell_number_out_dir)
    with mp.Pool(threads) as p:
        # results = p.starmap(rawdata.cell_number, [
        #     (barcode_clean_txt, sample_name, spike_ins, library_info, str(cell_number_out_dir), plot, debug)
        #     for barcode_clean_txt, sample_name in barcode_clean_txt_files
        # ])
        results = parallel_process_with_timeout(
            pool=p,
            func=rawdata.cell_number,
            args_list=[
                (barcode_clean_txt, sample_name, spike_ins, library_info, str(cell_number_out_dir), plot, debug)
                for barcode_clean_txt, sample_name in barcode_clean_txt_files
            ],
            timeout=1800
        )
    success = [res[0] for res in results]
    dfs = [res[1] for res in results]
    filtered_dfs = [res[2] for res in results]
    spike_counts_dfs = [res[3] for res in results]
    if not all(success):
        err_msg = f"ERROR: Some processes during Cell-Number Task failed. Check the logs for more details. Exiting."
        cell_count_files = tools.get_list_of_files(cell_number_out_dir, "cell_count", debug)
        finished_samples = [sample_name for _, sample_name in cell_count_files]
        failed_samples = [sample_name for sample_name in sample_names if sample_name not in finished_samples]
        log.logit(f"Failed samples: {failed_samples}", color="red")
        log.logit(err_msg, color="red")
        sys.exit(f"[err] {err_msg}")
    else:
        log.logit(f"Cell Number Calculations Complete. Combining results...", color="green")
        # Initialize empty DataFrames for combined results using the first result
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_filtered_df = pd.concat(filtered_dfs, ignore_index=True)
        combined_spike_counts = pd.concat(spike_counts_dfs, ignore_index=True)
        
        # Get the unique Sample_ID and Spike_Rsq values from combined_df
        unique_samples = combined_df['Sample_Rsq'].unique()
        split_samples = [x.split('_') for x in unique_samples]
        Rsq = pd.DataFrame(split_samples, columns=['Sample_ID', 'Spike_Rsq'])
        Rsq['Spike_Rsq'] = pd.to_numeric(Rsq['Spike_Rsq'])
        if debug: print(Rsq)
        if plot: plotting.rsq_per_sample(Rsq, out_dir)
        
        # Select columns and get unique combinations
        unique_df = combined_df[['reading_depth', 'Sample_ID']].drop_duplicates().reset_index(drop=True)
        if debug: print(unique_df)
        if plot: plotting.reading_depth_per_sample(unique_df, out_dir)
        
        aggr_df = (combined_df.groupby('Sample_ID')
           .agg({'reading_depth': 'mean', 'Spike_Rsq': 'mean'})
           .rename(columns={'reading_depth': 'mean_ReadingDepth', 
                          'Spike_Rsq': 'mean_Rsq'})
           .reset_index())
        aggr_df['cellsPerRead'] = 1/aggr_df['mean_ReadingDepth']
        if plot: plotting.mean_rsq_distribution(aggr_df, out_dir)
        if plot: plotting.mean_reading_depth_distribution(aggr_df, out_dir)
        
        log.logit(f"Results written to {out_dir}")
        Rsq.to_csv(f"{out_dir}/RsqPerSample.csv", index=False)
        unique_df.to_csv(f"{out_dir}/ReadingDepthPerSample.csv", index=False)
        aggr_df.to_csv(f"{out_dir}/CountInfoPerSample.csv", index=False)
        combined_df.to_csv(f"{out_dir}/UnfilteredSampleInfo.csv", index=False)
        combined_filtered_df.to_csv(f"{out_dir}/FilteredSampleInfo.csv", index=False)
        combined_spike_counts.to_csv(f"{out_dir}/SpikeInInfo.csv", index=False)
    log.logit(f"Calculating cell counts complete.", color="green")
        