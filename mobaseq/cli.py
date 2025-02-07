import os, signal
import click
from clint.textui import puts, colored, indent
import mobaseq.handlers as handlers
import mobaseq.qc.plot as plotting
import mobaseq.utils.logger as log
import mobaseq.process.tools as tools

from mobaseq.version import __version__
#import ch.utils.logger as log

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def cli():
    '''Description.'''
    # to make this script/module behave nicely with unix pipes
    # http://newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

@cli.group('process', short_help="Set of tools to perform basic read processing for MOBAseq.")
def process():
    """
    Suite of processing tools available for use to process MOBAseq data.
    """
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

@process.command('count-reads', short_help="Perform countreads on a single file or batch.")
@click.option('--read_one', '-r1', type=click.Path(exists=True), required=False, help="Individual R1 FASTQ File")
@click.option('--read_two', '-r2', type=click.Path(exists=True), required=False, help="Individual R2 FASTQ File")
@click.option('--input_dir', '-i', type=click.Path(exists=True), required=False, help="Directory with all R1/R2 FASTQ files. Use this option for BATCH Processing")
@click.option('--sample_name', '-s', type=click.STRING, required=False, help="Sample name")
@click.option('--sgid_file', '-f', type=click.Path(exists=True), required=True, help="sgID Excel File")
@click.option('--out_dir', '-o', type=click.Path(), required=True, help="Output directory to place CountRead results")
@click.option('--threads', '-t', type=click.INT, default=1, required=False, show_default=True, help="Number of threads to use for parallel processing")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=False, help="Print extra debugging output")
def countreads(read_one, read_two, input_dir, sample_name, sgid_file, out_dir, threads, debug):
    out_dir = tools.ensure_abs_path(out_dir) # Ensure out_dir is an absolute path
    if (read_one and read_two) and not input_dir: # Single file processing
        handlers.countreads_single(read_one, read_two, sample_name, sgid_file, str(out_dir), debug)
        log.logit(f"---> Successfully performed read counts for sample: {os.path.basename(read_one)} and {os.path.basename(read_two)}", color="green")
    elif input_dir and not (read_one or read_two): # Batch processing
        handlers.countreads_batch(input_dir, sgid_file, str(out_dir), threads, debug)
        log.logit(f"---> Successfully performed read counts for batch inside {input_dir}", color="green")
    else:
        raise click.UsageError("You must provide either --read_one and --read_two or --input_dir, but not both.")
    
@process.command('clean-barcodes', short_help="Clean barcodes from count-read files.")
@click.option('--merge_reads_csv', '-m', type=click.Path(exists=True), required=False, help="Individual _MergeReadOut.csv file")
@click.option('--input_dir', '-i', type=click.Path(exists=True), required=False, help="Input directory containing all processed CountRead files (_MergeReadOut.csv)")
@click.option('--out_dir', '-o', type=click.Path(), required=True, help="Output directory to place cleaned results")
@click.option('--threads', '-t', type=click.INT, default=1, required=False, show_default=True, help="Number of threads to use for parallel processing")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=False, help="Print extra debugging output")
def clean_barcodes(merge_reads_csv, input_dir, out_dir, threads, debug):
    out_dir = tools.ensure_abs_path(out_dir)
    if merge_reads_csv and not input_dir:
        handlers.clean_barcodes_single(merge_reads_csv, threads, str(out_dir), debug)
        log.logit(f"---> Successfully cleaned barcodes for file: {merge_reads_csv}", color="green")
    elif input_dir and not merge_reads_csv:
        handlers.clean_barcodes_batch(input_dir, str(out_dir), threads, debug)
        log.logit(f"---> Successfully cleaned barcodes inside {input_dir}", color="green")
    else:
        raise click.UsageError("You must provide either --merge_reads_csv or --input_dir, but not both.")
    
@process.command('sgID-qc', short_help="Summarize basic information from raw count files.")
@click.option('--merge_reads_csv', '-m', type=click.Path(exists=True), required=False, help="Individual _MergeReadOut.csv file")
@click.option('--input_dir', '-i', type=click.Path(exists=True), required=False, help="Input directory containing all processed CountRead files (_MergeReadOut.csv)")
@click.option('--sgid_file', '-f', type=click.Path(exists=True), required=True, help="sgID Excel File")
@click.option('--out_dir', '-o', type=click.Path(), required=True, help="Output directory to place summarized results")
@click.option('--threads', '-t', type=click.INT, default=1, required=False, show_default=True, help="Number of threads to use for parallel processing")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=False, help="Print extra debugging output")
def sgID_qc(merge_reads_csv, input_dir, sgid_file, out_dir, threads, debug):
    out_dir = tools.ensure_abs_path(out_dir) # Ensure out_dir is an absolute path
    if merge_reads_csv and not input_dir: # Single file processing
        handlers.sgID_qc_single(merge_reads_csv, sgid_file, str(out_dir), debug)
        log.logit(f"---> Successfully summarized raw counts for file: {merge_reads_csv}", color="green")
    elif input_dir and not merge_reads_csv: # Batch processing
        handlers.sgID_qc_batch(input_dir, sgid_file, str(out_dir), threads, debug)
        log.logit(f"---> Successfully summarized raw counts inside {input_dir}", color="green")
    else:
        raise click.UsageError("You must provide either --merge_reads_csv or --input_dir, but not both.")
    
@process.command('mapped-reads', short_help="Get mapped and unmapped reads.")
@click.option('--merge_reads_csv', '-m', type=click.Path(exists=True), required=False, help="Individual _MergeReadOut.csv file")
@click.option('--input_dir', '-i', type=click.Path(exists=True), required=False, help="Input directory containing all processed CountRead files (_MergeReadOut.csv)")
@click.option('--sgid_file', '-f', type=click.Path(exists=True), required=True, help="sgID Excel File")
@click.option('--out_dir', '-o', type=click.Path(), required=True, help="Output directory to place summarized results")
@click.option('--threads', '-t', type=click.INT, default=1, required=False, show_default=True, help="Number of threads to use for parallel processing")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=False, help="Print extra debugging output")
def mapped_reads(merge_reads_csv, input_dir, sgid_file, out_dir, threads, debug):
    out_dir = tools.ensure_abs_path(out_dir)
    if merge_reads_csv and not input_dir: # Single file processing
        handlers.mapped_reads_single(merge_reads_csv, sgid_file, str(out_dir), debug)
        log.logit(f"---> Successfully summarized mapped and unmapped reads for file: {merge_reads_csv}", color="green")
    elif input_dir and not merge_reads_csv: # Batch processing
        handlers.mapped_reads_batch(input_dir, sgid_file, str(out_dir), threads, debug)
        log.logit(f"---> Successfully summarized mapped and unmapped reads inside {input_dir}", color="green")
    else:
        raise click.UsageError("You must provide either --merge_reads_csv or --input_dir, but not both.")

# Add the first-level subcommand group to the main command group
cli.add_command(process)

@cli.group('plot', short_help="Set of tools to plot and visualize MOBAseq data.")
def plot():
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

@plot.command('mapped-reads', short_help="Plot mapped and unmapped reads.")
@click.option('--mapped_percentages_csv', '-i', type=click.Path(exists=True), required=True, help="The _mapped_percentages.csv file to plot")
@click.option('--out_dir', '-o', type=click.Path(), required=True, help="Output directory to place summarized results")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=False, help="Print extra debugging output")
def plot_mapped_reads(mapped_percentages_csv, out_dir, debug):
    out_dir = tools.ensure_abs_path(out_dir) # Ensure out_dir is an absolute path
    plotting.mapped_reads(mapped_percentages_csv, str(out_dir), debug)
    log.logit(f"---> Successfully plotted mapped and unmapped reads. Plot can be found: {out_dir}/mapped_percentages.png", color="green")

@plot.command('per-sgID', short_help="Plot per sgID plots.")
@click.option('--sgID_info_csv', '-i', type=click.Path(exists=True), required=True, help="The _sgID_info.csv file to plot")
@click.option('--out_dir', '-o', type=click.Path(), required=True, help="Output directory to place summarized results")
@click.option('--no_dummy', '-nd', is_flag=True, show_default=True, default=False, required=False, help="Do not include dummy sgIDs in the plots")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=False, help="Print extra debugging output")
def plot_per_sgid(sgID_info_csv, out_dir, no_dummy, debug):
    out_dir = tools.ensure_abs_path(out_dir)
    plotting.per_sgid(sgID_info_csv, str(out_dir), no_dummy, debug)
    log.logit(f"---> Successfully plotted per sgID plots. Plots can be found inside: {out_dir}", color="green")


cli.add_command(plot)

@cli.group('legacy', short_help="Original code written by Andy Xu.")
def legacy():
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

@legacy.command('count-reads', short_help="Perform countreads on a single file or batch.")
@click.option('--read_one', '-r1', type=click.Path(exists=True), required=False, help="Individual R1 FASTQ File")
@click.option('--read_two', '-r2', type=click.Path(exists=True), required=False, help="Individual R2 FASTQ File")
@click.option('--input_dir', '-i', type=click.Path(exists=True), required=False, help="Directory with all R1/R2 FASTQ files. Use this option for BATCH Processing")
@click.option('--sample_name', '-s', type=click.STRING, required=False, help="Sample name")
@click.option('--sgid_file', '-f', type=click.Path(exists=True), required=True, help="sgID Excel File")
@click.option('--out_dir', '-o', type=click.Path(), required=True, help="Output directory to place CountRead results")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=False, help="Print extra debugging output")
def countreads(read_one, read_two, input_dir, sample_name, sgid_file, out_dir, debug):
    import mobaseq.process.legacy as legacy
    import mobaseq.process.tools as tools
    
    out_dir = tools.ensure_abs_path(out_dir) # Ensure out_dir is an absolute path
    if (read_one and read_two) and not input_dir:
        # Single file processing
        sample_name = tools.process_fastq_files(read_one, read_two, debug)
        sgID_dict = tools.process_sgid_file(sgid_file, debug)
        min_length, max_length, all_start_with_G = tools.get_info_about_sgID_file(sgID_dict, debug)
        legacy.countreads(read_one, read_two, sample_name, sgID_dict, min_length, max_length, all_start_with_G, str(out_dir), debug)
        log.logit(f"---> Successfully performed read counts for sample: {sample_name}", color="green")
    elif input_dir and not (read_one or read_two):
        # Batch processing
        list_of_fastqs = tools.get_list_of_files(input_dir, "fastq", debug)
        sgID_dict = tools.process_sgid_file(sgid_file, debug)
        min_length, max_length, all_start_with_G = tools.get_info_about_sgID_file(sgID_dict, debug)
        for fastq1, fastq2, sample_name in list_of_fastqs:
            legacy.countreads(fastq1, fastq2, sample_name, sgID_dict, min_length, max_length, all_start_with_G, str(out_dir), debug)
        log.logit(f"---> Successfully performed read counts for batch inside {input_dir}", color="green")
    else:
        raise click.UsageError("You must provide either --read_one and --read_two or --input_dir, but not both.")
    
@legacy.command('clean-barcodes', short_help="Clean barcodes from count-read files.")
@click.option('--input_dir', '-i', type=click.Path(exists=True), required=True, help="Input directory containing all processed CountRead files (_MergeReadOut.csv)")
@click.option('--out_dir', '-o', type=click.Path(), required=True, help="Output directory to place cleaned results")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=False, help="Print extra debugging output")
def clean_barcodes(input_dir, out_dir, debug):
    import mobaseq.process.legacy as legacy
    import mobaseq.process.tools as tools
    if not os.path.isabs(out_dir): out_dir = os.path.abspath(out_dir)
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    if not os.path.isabs(input_dir): input_dir = os.path.abspath(input_dir)
    legacy.clean_barcodes(input_dir, out_dir, debug)
    log.logit(f"---> Successfully cleaned barcodes inside {input_dir}", color="green")

@legacy.command('sgID-qc', short_help="Summarize basic information from raw count files.")
@click.option('--input_dir', '-i', type=click.Path(exists=True), required=True, help="Input directory containing all processed CountRead files (_MergeReadOut.csv)")
@click.option('--sgid_file', '-f', type=click.Path(exists=True), required=True, help="sgID Excel File")
@click.option('--out_dir', '-o', type=click.Path(), required=True, help="Output directory to place summarized results")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=False, help="Print extra debugging output")
def summarize(input_dir, sgid_file, out_dir, debug):
    import mobaseq.process.legacy as legacy
    import mobaseq.process.tools as tools
    sgID_dict = tools.process_sgid_file(sgid_file, debug)
    if not os.path.isabs(out_dir): out_dir = os.path.abspath(out_dir) # Make sure the output directory is a full path
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    if not os.path.isabs(input_dir): input_dir = os.path.abspath(input_dir)
    legacy.summarize(input_dir, sgID_dict, out_dir, debug)
    log.logit(f"---> Successfully summarized raw counts inside {input_dir}", color="green")

@legacy.command('plot-mapped-and-unmapped', short_help="Plot mapped and unmapped reads.")
@click.option('--input_dir', '-i', type=click.Path(exists=True), required=True, help="Input directory containing all processed CountRead files (_MergeReadOut.csv)")
@click.option('--sgid_file', '-f', type=click.Path(exists=True), required=True, help="sgID Excel File")
@click.option('--out_dir', '-o', type=click.Path(), required=True, help="Output directory to place summarized results")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=False, help="Print extra debugging output")
def plot_mapped_and_unmapped(input_dir, sgid_file, out_dir, debug):
    import mobaseq.process.legacy as legacy
    import mobaseq.process.tools as tools
    sgID_dict = tools.process_sgid_file(sgid_file, debug)
    if not os.path.isabs(out_dir): out_dir = os.path.abspath(out_dir)
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    if not os.path.isabs(input_dir): input_dir = os.path.abspath(input_dir)
    legacy.plot_mapped_and_unmapped(input_dir, sgID_dict, out_dir, debug)
    log.logit(f"---> Successfully plotted mapped and unmapped reads inside {input_dir}", color="green")

@legacy.command('plot-per-sgID', short_help="Plot per sgID plots.")
@click.option('--reads_csv', '-r', type=click.Path(exists=True), required=True, help="Total reads CSV file")
@click.option('--barcodes_csv', '-b', type=click.Path(exists=True), required=True, help="Barcodes CSV file")
@click.option('--rel_reads_csv', '-rr', type=click.Path(exists=True), required=True, help="Relative reads CSV file")
@click.option('--rel_barcodes_csv', '-rb', type=click.Path(exists=True), required=True, help="Relative barcodes CSV file")
@click.option('--sgid_file', '-f', type=click.Path(exists=True), required=True, help="sgID Excel File")
@click.option('--out_dir', '-o', type=click.Path(), required=True, help="Output directory to place summarized results")
@click.option('--no_dummy', '-nd', is_flag=True, show_default=True, default=False, required=False, help="Do not include dummy sgIDs in the plots")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=False, help="Print extra debugging output")
def plot_per_sgid(reads_csv, barcodes_csv, rel_reads_csv, rel_barcodes_csv, sgid_file, out_dir, no_dummy, debug):
    import mobaseq.process.legacy as legacy
    import mobaseq.process.tools as tools
    if not os.path.isabs(out_dir): out_dir = os.path.abspath(out_dir)
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    sgIDs = tools.process_sgid_file(sgid_file, debug)
    legacy.plot_per_sgid(reads_csv, barcodes_csv, rel_reads_csv, rel_barcodes_csv, sgIDs, out_dir, no_dummy, debug)
    log.logit(f"---> Successfully plotted per sgID plots inside {out_dir}", color="green")

cli.add_command(legacy)

@cli.command('do-something', short_help="Should check a few things and output: Hello World! :D")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=False, help="Print extra debugging output")
def do_something(debug):
    """Test Command."""
    log.logit(f"---> Successfully ran your first CLI script!", color="green")
    puts(colored.green(f"Hello World! :D"))
    #import ch.vdbtools.importer as importer
    #importer.import_samples(samples, sample_duckdb, batch_number, debug, clobber)
    #log.logit(f"---> Successfully ran your first CLI script!", color="green")
    #puts(colored.green(f"Hello World! :D"))
