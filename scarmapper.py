"""
@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2023
"""
import csv
import datetime
import glob
import itertools
import os
import collections
import subprocess
import argparse
import sys
import time
import pathos
from natsort import natsort
from scipy import stats
from distutils.util import strtobool
from scipy.stats import gmean
from Valkyries import Tool_Box, Version_Dependencies as VersionDependencies, FASTQ_Tools
from scarmapper import ScarMapperPlot
import re
import pathlib
import platform

# This is a seriously ugly hack to check the existence and age of the compiled file.
folder_content = os.listdir("{0}{1}scarmapper{1}".format(pathlib.Path(__file__).parent.absolute(), os.sep))
python_ver = platform.python_version()

subver = 7
if "3.8.0" <= python_ver < "3.9.0":
    subver = 8
elif "3.9.0" <= python_ver < "3.10.0":
    subver = 9
elif "3.10.0" <= python_ver < "3.11.0":
    subver = 10

regex = re.compile("SlidingWindow.cpython-3{}.*.so".format(subver))
cfile = ""
old_file = False
for f in folder_content:
    if regex.search(f):
        cfile = f
        break
if cfile:
    cpath = "{0}{1}scarmapper{1}{2}".format(pathlib.Path(__file__).parent.absolute(), os.sep, cfile)
    pyx_file = "{0}{1}scarmapper{1}SlidingWindow.pyx".format(pathlib.Path(__file__).parent.absolute(), os.sep)
    compiled_time = pathlib.Path(cpath).stat().st_ctime
    pyx_module_time = pathlib.Path(pyx_file).stat().st_ctime

    if pyx_module_time >= compiled_time:
        old_file = True

if not cfile or old_file:
    print("Compiled Module Doesn't Exist or is Old; Compiling New SlidingWindow Module")
    
    # Save current directory
    original_dir = os.getcwd()
    
    # Change to the parent directory (where scarmapper_AL.py is)
    parent_dir = pathlib.Path(__file__).parent.absolute()
    os.chdir(parent_dir)
    
    # Run setup.py from parent directory, pointing to the setup.py inside scarmapper
    setup_command = f"python3.{subver} scarmapper/setup.py build_ext --inplace"
    subprocess.run(setup_command, shell=True)
    
    # Return to original directory
    os.chdir(original_dir)
    
    # The sleep is to allow for network or disk latency.
    time.sleep(5.0)
    
from scarmapper import INDEL_Processing as Indel_Processing, TargetMapper as Target_Mapper

__author__ = 'Dennis A. Simpson'
__version__ = '2.0.0 BETA'
__package__ = 'ScarMapper'


def pear_consensus(args, log, fq1=None, fq2=None, sample_prefix=None):
    """
    This will take the input FASTQ files and use PEAR to generate a consensus file.
    @param fq2:
    @param fq1:
    @param args:
    @param log:
    @param sample_prefix: Optional prefix for output files (for batch mode)
    @return:
    """
    log.info("Beginning PEAR Consensus")
    
    # Check if FASTQ files were provided as parameters (batch mode)
    if fq1:
        fastq1 = fq1
        fastq2 = fq2
    else:
        # Fall back to args attributes (single sample mode)
        fastq1 = args.FASTQ1
        fastq2 = args.FASTQ2
    
    # Use custom prefix if provided (batch mode), otherwise use Job_Name
    if sample_prefix:
        fastq_consensus_prefix = "{}{}".format(args.WorkingFolder, sample_prefix)
    else:
        fastq_consensus_prefix = "{}{}".format(args.WorkingFolder, args.Job_Name)
    fastq_consensus_file = "{}.assembled.fastq".format(fastq_consensus_prefix)
    discarded_fastq = "{}.discarded.fastq".format(fastq_consensus_prefix)
    r1_unassembled = "{}.unassembled.forward.fastq".format(fastq_consensus_prefix)
    r2_unassembled = "{}.unassembled.reverse.fastq".format(fastq_consensus_prefix)

    y = "-y {} ".format(args.Memory)
    j = "-j {} ".format(args.Spawn)

    p_value = ''
    if args.PValue:
        p_value = "-p {} ".format(args.PValue)
    min_overlap = ''
    if args.MinOverlap:
        min_overlap = "-v {} ".format(args.MinOverlap)
    quality_threshold = ""
    if args.QualityThreshold:
        quality_threshold = "-q {} ".format(args.QualityThreshold)
    phred_value = ""
    if args.PhredValue:
        phred_value = "-b {} ".format(args.PhredValue)
    test_method = ""
    if args.TestMethod:
        test_method = "-g {}".format(args.TestMethod)
    n = ""
    if args.MinConsensusLength:
        n = "-n {} ".format(args.MinConsensusLength)

    proc = subprocess.run(
        "{0}{1}Pear{1}bin{1}./pear -f {2} -r {3} -o {4} {5}{6}{7}{8}{9}{10}{11}"
        .format(pathlib.Path(__file__).parent.absolute(), os.sep, fastq1, fastq2, fastq_consensus_prefix, y, j, n,
                p_value, min_overlap, quality_threshold, phred_value, test_method),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    if proc.stderr:
        log.error("{}\n{}\n".format(proc.stderr.decode(), proc.stdout.decode()))
        return
    else:
        log.info(
        "Begin PEAR Output\n"
        "----------------------------------------------------------------------------------------------------------\n{}"
        "\n----------------------------------------------------------------------------------------------------------\n"
         .format(proc.stdout.decode()))

    file_list = [fastq_consensus_file, r1_unassembled, r2_unassembled]

    if pathlib.Path(discarded_fastq).exists():
        file_list.append(discarded_fastq)
    else:
        Tool_Box.delete([discarded_fastq])

    return file_list


def preprocess_bad_fastq(args, log):

    fq2 = FASTQ_Tools.FASTQ_Reader(args.FASTQ2, log)
    fq1 = FASTQ_Tools.FASTQ_Reader(args.FASTQ1, log)

    # MSK specific code
    if pathlib.Path("{}{}_corrected_R1.fastq".format(args.WorkingFolder, args.Job_Name)).exists():
        Tool_Box.delete(["{}{}_corrected_R1.fastq".format(args.WorkingFolder, args.Job_Name),
                         "{}{}_corrected_R2.fastq".format(args.WorkingFolder, args.Job_Name)])
    outfile1 = open("{}{}_corrected_R1.fastq".format(args.WorkingFolder, args.Job_Name), "a")
    outfile2 = open("{}{}_corrected_R2.fastq".format(args.WorkingFolder, args.Job_Name), "a")
    outstring1 = ""
    outstring2 = ""
    read_counter = 0
    good_reads = 0
    filtered_count = 0
    eof = False

    # This block was created to clean up some FASTQ files that were not filtered well off the sequencer.
    while not eof:
        try:
            read_counter += 1
            fastq1_read = next(fq1.seq_read())
            fastq2_read = next(fq2.seq_read())
            len_r1 = len(fastq1_read.seq)
            len_r2 = len(fastq2_read.seq)

            # Unbalanced short reads are crashing Pear.
            if len_r1 != len_r2 or fastq1_read.seq.count("N") > 5 or fastq2_read.seq.count("N") > 5:
                filtered_count += 1
                continue

            # if "GGCTCCATCG" in fastq1_read.seq:
            # elif "ATGTGGCTCT" in fastq1_read.seq:

            if "CTGGCTCCA" in fastq1_read.seq:
                # These are our forward reads in R1
                outstring1 += "@{}\n{}\n+\n{}\n".format(fastq1_read.name, fastq1_read.seq, fastq1_read.qual)
                outstring2 += "@{}\n{}\n+\n{}\n".format(fastq2_read.name, fastq2_read.seq, fastq2_read.qual)
                good_reads += 1

            elif "TGTGGCTCTG" in fastq1_read.seq:
                # These are reverse reads in R1
                good_reads += 1
                outstring2 += "@{}\n{}\n+\n{}\n".format(fastq1_read.name, fastq1_read.seq, fastq1_read.qual)
                outstring1 += "@{}\n{}\n+\n{}\n".format(fastq2_read.name, fastq2_read.seq, fastq2_read.qual)

        except StopIteration:
            eof = True

        # Write to FASTQ files every 1 million reads.
        if read_counter % 1000000 == 0:
            log.info("Processed {} reads.  {} were filtered and {} were accepted."
                     .format(read_counter, filtered_count, good_reads))
            outfile1.write(outstring1)
            outfile2.write(outstring2)
            outstring1 = ""
            outstring2 = ""

    # Write any remaining reads to FASTQ files and close them.

    outfile1.write(outstring1)
    outfile2.write(outstring2)
    outfile1.close()
    outfile2.close()

    fq1 = "{}{}_corrected_R1.fastq".format(args.WorkingFolder, args.Job_Name)
    fq2 = "{}{}_corrected_R2.fastq".format(args.WorkingFolder, args.Job_Name)
    log.info("{} Filtered Reads.  {} Accepted reads. {} Total reads."
             .format(filtered_count, good_reads, read_counter))
    return fq1, fq2

def process_single_sample(args, log, sample_info, target_mapper, version, run_start):
    """
    Process a single pre-demultiplexed sample.
    
    @param args: Command line arguments
    @param log: Logger object
    @param sample_info: Dictionary with sample information
    @param target_mapper: TargetMapper object
    @param version: Version string
    @param run_start: Run start timestamp
    @return: Summary data for the sample
    """
    import pathlib
    
    index_name = sample_info['Index_Name']
    fastq1_path = sample_info['FASTQ1_Path']
    fastq2_path = sample_info.get('FASTQ2_Path', '')
    sample_hr_donor = sample_info.get('HR_Donor', '')
    
    log.info(f"Processing sample: {index_name}")
    if sample_hr_donor:
        log.info(f"Sample-specific HR_Donor: {sample_hr_donor}")
    
    # Check if paired-end or single-end
    is_paired_end = fastq2_path and pathlib.Path(fastq2_path).exists()
    
    file_list = []
    
    if is_paired_end:
        log.info(f"Paired-end reads detected for {index_name}. Running PEAR consensus.")
        # Run PEAR with custom file names and unique prefix for this sample
        file_list = pear_consensus(args, log, fq1=fastq1_path, fq2=fastq2_path, sample_prefix=index_name)
        
        if not file_list:
            log.error(f"PEAR failed for {index_name}. Skipping.")
            return None
            
        fastq_consensus = file_list[0]
        fq1 = FASTQ_Tools.FASTQ_Reader(fastq_consensus, log)
        fq2 = None
    else:
        log.info(f"Single-end reads detected for {index_name}. Skipping PEAR.")
        fq1 = FASTQ_Tools.FASTQ_Reader(fastq1_path, log)
        fq2 = None
    
    # Create a mock index_dict for this single sample
    # Include sample-specific HR_Donor in position 8
    index_dict = {
        index_name: [
            '',  # right_index_sequence (not used in pre-demux mode)
            0,   # right_index_count
            '',  # left_index_sequence (not used in pre-demux mode)
            0,   # left_index_count
            index_name,
            sample_info['Sample_Name'],
            sample_info['Replicate'],
            sample_info['Locus'],
            sample_hr_donor  # Sample-specific HR_Donor
        ]
    }
    
    # Create sequence list from FASTQ file
    sequence_list = []
    log.info(f"Reading sequences from {index_name}")
    
    eof = False
    read_count = 0
    while not eof:
        try:
            fastq_read = next(fq1.seq_read())
            sequence_list.append(fastq_read.seq)
            read_count += 1
            
            if read_count % 100000 == 0:
                log.info(f"Read {read_count} sequences for {index_name}")
                
        except StopIteration:
            eof = True
    
    log.info(f"Total sequences read for {index_name}: {read_count}")
    
    # Process the sample using ScarSearch
    scar_search = Indel_Processing.ScarSearch(
        log, args, version, run_start,
        target_mapper.targets,
        index_dict,
        index_name,
        sequence_list,
        read_count,  # indexed_read_count
        args.PatternThreshold * 0.00001  # lower_limit_count
    )
    
    # Clean up PEAR files if needed
    if is_paired_end and file_list:
        if args.DeleteConsensusFASTQ:
            log.info(f"Deleting PEAR FASTQ Files for {index_name}")
            Tool_Box.delete(file_list)
        else:
            log.info(f"Compressing PEAR FASTQ Files for {index_name}")
            import pathos
            import itertools
            p = pathos.multiprocessing.Pool(args.Spawn)
            p.starmap(Tool_Box.compress_files, zip(file_list, itertools.repeat(log)))
    
    return scar_search.summary_data


def batch_process_samples(args, log, version, run_start):
    """
    Process multiple pre-demultiplexed samples from Sample Manifest.
    
    @param args: Command line arguments
    @param log: Logger object
    @param version: Version string
    @param run_start: Run start timestamp
    """
    log.info("Starting batch processing of pre-demultiplexed samples")
    
    # Parse sample manifest
    sample_manifest = Tool_Box.FileParser.indices(log, args.SampleManifest)
    
    # Create target mapper
    target_mapper = Target_Mapper.TargetMapper(log, args, sample_manifest)
    
    # Build list of samples to process
    samples_to_process = []
    for sample_row in sample_manifest:
        # Determine column positions based on manifest format
        # Format: Index_Name, Sample_Name, Replicate, Forward_Phase, Reverse_Phase, Locus, FASTQ1, FASTQ2, HR_Donor
        sample_info = {
            'Index_Name': sample_row[0],
            'Sample_Name': sample_row[1],
            'Replicate': sample_row[2],
            'Forward_Phase': sample_row[3] if len(sample_row) > 3 else '',
            'Reverse_Phase': sample_row[4] if len(sample_row) > 4 else '',
            'Locus': sample_row[5] if len(sample_row) > 5 else '',
            'FASTQ1_Path': sample_row[6] if len(sample_row) > 6 else '',
            'FASTQ2_Path': sample_row[7] if len(sample_row) > 7 else '',
            'HR_Donor': sample_row[8].upper().strip() if len(sample_row) > 8 and sample_row[8] else ''
        }
        
        # Validate FASTQ1 exists
        if not pathlib.Path(sample_info['FASTQ1_Path']).exists():
            log.error(f"FASTQ1 file not found for {sample_info['Index_Name']}: {sample_info['FASTQ1_Path']}")
            continue
            
        samples_to_process.append(sample_info)
    
    log.info(f"Found {len(samples_to_process)} samples to process")
    
    # Process samples (can be done in parallel or serial)
    all_summary_data = []
    
    if args.Spawn > 1:
        log.info(f"Processing samples in parallel with {args.Spawn} processes")
        import pathos
        p = pathos.multiprocessing.Pool(args.Spawn)
        
        # Prepare arguments for parallel processing
        process_args = [
            (args, log, sample_info, target_mapper, version, run_start)
            for sample_info in samples_to_process
        ]
        
        all_summary_data = p.starmap(process_single_sample, process_args)
    else:
        log.info("Processing samples serially")
        for sample_info in samples_to_process:
            summary_data = process_single_sample(args, log, sample_info, target_mapper, version, run_start)
            if summary_data:
                all_summary_data.append(summary_data)
    
    # Filter out None results (failed samples)
    all_summary_data = [sd for sd in all_summary_data if sd is not None]
    
    # Write combined summary file
    write_batch_summary(args, log, all_summary_data, version, run_start, samples_to_process)
    
    log.info("Batch processing complete")


def write_batch_summary(args, log, all_summary_data, version, run_start, samples_info):
    """
    Write a summary file for batch-processed samples.
    
    @param args: Command line arguments
    @param log: Logger object
    @param all_summary_data: List of summary data from all samples
    @param version: Version string
    @param run_start: Run start timestamp
    @param samples_info: List of sample information dictionaries
    """
    import datetime
    
    log.info("Writing batch summary file")
    
    summary_file = open(f"{args.WorkingFolder}{args.Job_Name}_Batch_Summary.txt", "w")
    
    run_stop = datetime.datetime.today().strftime("%a %b %d %H:%M:%S %Y")
    
    summary_outstring = f"ScarMapper Batch Processing v{version}\n"
    summary_outstring += f"Start: {run_start}\nEnd: {run_stop}\n\n"
    
    # Check if any sample has HR_Donor
    has_hr = any(info.get('HR_Donor', '') for info in samples_info)
    
    summary_outstring += "Index_Name\tSample_Name\tReplicate\tTarget\tTotal_Reads\t"
    summary_outstring += "Passing_Filters\tScar_Count\tScar_Fraction\t"
    if has_hr:
        summary_outstring += "HR_Donor\tHR_Count\tHR_Fraction\t"
    summary_outstring += "TMEJ\tNHEJ\tNon-MH_Deletion\tInsertion\tSNV\n"
    
    # Create lookup dictionaries for sample info by index_name
    hr_donor_lookup = {info['Index_Name']: info.get('HR_Donor', '') for info in samples_info}
    sample_name_lookup = {info['Index_Name']: info.get('Sample_Name', info['Index_Name']) for info in samples_info}
    replicate_lookup = {info['Index_Name']: info.get('Replicate', '1') for info in samples_info}
    
    for summary_data in all_summary_data:
        if not summary_data:
            continue
            
        index_name = summary_data[0]
        passing_filters = summary_data[1]
        target = summary_data[9]
        junction_data = summary_data[8]
        
        scar_count = passing_filters - summary_data[6][0] - summary_data[6][1]
        scar_fraction = scar_count / passing_filters if passing_filters > 0 else 0
        
        # Extract repair pathway counts
        tmej = junction_data[0] if junction_data else 0
        nhej = junction_data[1] if junction_data else 0
        non_mh_del = junction_data[4] if junction_data else 0
        insertion = junction_data[2] if junction_data else 0
        snv = junction_data[5] if junction_data else 0
        
        # Get sample name and replicate from manifest
        sample_name = sample_name_lookup.get(index_name, index_name)
        replicate = replicate_lookup.get(index_name, '1')
        
        summary_outstring += f"{index_name}\t{sample_name}\t{replicate}\t{target}\t{passing_filters}\t"
        summary_outstring += f"{passing_filters}\t{scar_count}\t{scar_fraction:.4f}\t"
        
        if has_hr:
            hr_donor = hr_donor_lookup.get(index_name, '')
            hr_left = summary_data[10][0] if len(summary_data) > 10 else 0
            hr_right = summary_data[10][1] if len(summary_data) > 10 else 0
            hr_total = hr_left + hr_right
            hr_freq = hr_total / passing_filters if passing_filters > 0 else 0
            summary_outstring += f"{hr_donor}\t{hr_total}\t{hr_freq:.4f}\t"
        
        summary_outstring += f"{tmej}\t{nhej}\t{non_mh_del}\t{insertion}\t{snv}\n"
    
    summary_file.write(summary_outstring)
    summary_file.close()
    
    log.info(f"Batch summary written to {args.WorkingFolder}{args.Job_Name}_Batch_Summary.txt")


def main(command_line_args=None):
    """
    Let's get this party started.
    @param command_line_args
    """
    start_time = time.time()
    VersionDependencies.python_check()

    if not command_line_args:
        command_line_args = sys.argv

    run_start = datetime.datetime.today().strftime("%H:%M:%S %Y  %a %b %d")
    parser = argparse.ArgumentParser(description="A package to map genomic repair scars at defined loci.\n {} v{}"
                                     .format(__package__, __version__), formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--options_file', action='store', dest='options_file', required=True,
                        help='File containing program parameters.')

    # Check options file for errors and return object.
    args = error_checking(string_to_boolean(parser))

    log = Tool_Box.Logger(args)
    Tool_Box.log_environment_info(log, args, command_line_args)

    module_name = ""
    log.info("{} v{}".format(__package__, __version__))

    # Check if we're in batch mode (pre-demultiplexed samples)
    if hasattr(args, 'BatchMode') and args.BatchMode:
        log.info("Running in Batch Mode for pre-demultiplexed samples")
        batch_process_samples(args, log, __version__, run_start)
        
    elif args.IndelProcessing:
        # Original processing logic for multiplexed FASTQ files
        file_list = []
        if args.Platform == "Illumina" or args.Platform == "Ramsden" or args.Platform == "TruSeq":
            log.info("Sending FASTQ files to FASTQ preprocessor.")

            if args.PEAR:
                file_list = pear_consensus(args, log)
                if not file_list:
                    log.error("PEAR failed.  Check logs.")
                    raise SystemExit(1)
                fastq_consensus = file_list[0]
                fq1 = FASTQ_Tools.FASTQ_Reader(fastq_consensus, log)
                fq2 = None

            else:
                fq2 = FASTQ_Tools.FASTQ_Reader(args.FASTQ2, log)
                fq1 = FASTQ_Tools.FASTQ_Reader(args.FASTQ1, log)

            sample_manifest = Tool_Box.FileParser.indices(log, args.SampleManifest)
            indel_processing = \
                Indel_Processing.DataProcessing(log, args, run_start, __version__,
                                                Target_Mapper.TargetMapper(log, args, sample_manifest), fq1, fq2)

            indel_processing.main_loop()

            # Compress or delete PEAR files.
            if args.PEAR and file_list:
                if args.DeleteConsensusFASTQ:
                    log.info("Deleting PEAR FASTQ Files.")
                    Tool_Box.delete(file_list)
                else:
                    log.info("Compressing {} FASTQ Files Generated by PEAR.".format(len(file_list)))
                    p = pathos.multiprocessing.Pool(args.Spawn)
                    p.starmap(Tool_Box.compress_files, zip(file_list, itertools.repeat(log)))
        else:
            log.error("Only 'Illumina', 'TruSeq' or 'Ramsden' --Platform methods currently allowed.")
            raise SystemExit(1)

    elif not args.IndelProcessing:
        # Run frequency file Combine module
        run_start = datetime.datetime.today().strftime("%a %b %d %H:%M:%S %Y")
        log.info("Process Replicates.")
        data_dict = collections.defaultdict(list)
        file_list = [f for f in glob.glob("{}*ScarMapper_Frequency.txt".format(args.DataFiles, ))]
        file_count = len(file_list)
        page_header = "# ScarMapper File Merge v{}\n# Run: {}\n# Sample Name: {}\n" \
            .format(__version__, run_start, args.SampleName)

        line_num = 0
        index_file = list(csv.reader(open(file_list[0]), delimiter='\t'))
        for line in index_file:
            if not line:
                break
            elif line_num > 3:
                page_header += "{}\n".format(line[0])

            line_num += 1
        page_header += "\n\n"

        for file_name in file_list:
            freq_file_data = Tool_Box.FileParser.indices(log, file_name)

            for row in freq_file_data:
                key = "{}|{}|{}|{}".format(row[3], row[4], row[6], row[8])
                row_data = row[2:]

                if key in data_dict:
                    data_dict[key][0].append(float(row[1]))
                else:
                    data_dict[key] = [[float(row[1])], row_data]

        # Process Data and Write Combined Frequency results file

        plot_data_dict = collections.defaultdict(list)
        label_dict = collections.defaultdict(float)
        output_data_dict = collections.defaultdict(list)
        marker_list = []

        for key, row_list in data_dict.items():
            # Require pattern to be in at least half of the files.
            if len(row_list[0]) / file_count >= 0.5:
                row_string = "\t".join(row_list[1])
                freq = gmean(row_list[0])
                sem = stats.sem(row_list[0])
                freq_results_outstring = "{}\t{}\t{}\n".format(freq, sem, row_string)
                output_key = freq

                # Freq is a 17 digit float making it very unlikely to be duplicated but if it is this increments it by
                # a small number and checks the uniqueness again.
                if output_key in output_data_dict:
                    output_key = output_key+1e-16
                    if output_key in output_data_dict:
                        output_key = output_key + 1e-16

                scar_type = row_list[1][0]
                label_dict[scar_type] += freq

                # Gather our data for plotting
                lft_del = int(row_list[1][1])
                rt_del = int(row_list[1][2])
                mh_size = int(row_list[1][5])
                ins_size = int(row_list[1][7])

                output_data_dict[output_key] = \
                    [(freq, lft_del, rt_del, mh_size, ins_size, scar_type), freq_results_outstring]

        freq_results_outstring = \
            "{}# Frequency\tSEM\tScar Type\tLeft Deletions\tRight Deletions\tDeletion Size\tMicrohomology\t" \
            "Microhomology Size\tInsertion\tInsertion Size\tLeft Template\tRight Template\tConsensus Left Junction\t" \
            "Consensus Right Junction\tTarget Left Junction\tTarget Right Junction\tConsensus\tTarget Region\n" \
            .format(page_header)

        # Now draw a pretty graph of the data if we are not dealing with a negative control.
        for k in natsort.natsorted(output_data_dict, reverse=True):
            data_list = output_data_dict[k]
            freq_results_outstring += data_list[1]

            freq = data_list[0][0]
            lft_del = data_list[0][1]
            rt_del = data_list[0][2]
            mh_size = data_list[0][3]
            ins_size = data_list[0][4]
            scar_type = data_list[0][5]

            # Plotting all scar patterns is messy.  This provides a cutoff.
            if freq < 0.0025:
                continue

            y_value = freq * 0.5
            lft_ins_width = freq
            rt_ins_width = freq

            # This is to find the largest value.  Used to set the x-axis limits.
            marker_list.extend([lft_del + (mh_size * 0.5), rt_del + (mh_size * 0.5), ins_size])

            # Deletion size included half the size of any microhomology present.
            lft_del_plot_value = (lft_del + (mh_size * 0.5)) * -1
            rt_del_plot_value = rt_del + (mh_size * 0.5)

            # Insertions are centered on 0, so we need to take half the value for each side.
            lft_ins_plot_value = (ins_size * 0.5) * -1
            rt_ins_plot_value = ins_size * 0.5

            # Scale the width of bars for insertions inside of deletions
            if lft_del + (mh_size * 0.5) != 0:
                lft_ins_width = freq * 0.5
            if rt_del + (mh_size * 0.5) != 0:
                rt_ins_width = freq * 0.5

            if scar_type not in plot_data_dict:
                plot_data_dict[scar_type] = \
                    [[freq], [lft_del_plot_value], [rt_del_plot_value], [lft_ins_plot_value],
                     [rt_ins_plot_value], [lft_ins_width], [rt_ins_width], [y_value]]
            else:
                # Get some previous plot data
                count = len(plot_data_dict[scar_type][0])
                previous_freq = plot_data_dict[scar_type][0][count - 1]
                previous_y = plot_data_dict[scar_type][7][count - 1]

                plot_data_dict[scar_type][0].append(freq)
                plot_data_dict[scar_type][1].append(lft_del_plot_value)
                plot_data_dict[scar_type][2].append(rt_del_plot_value)
                plot_data_dict[scar_type][3].append(lft_ins_plot_value)
                plot_data_dict[scar_type][4].append(rt_ins_plot_value)
                plot_data_dict[scar_type][5].append(lft_ins_width)
                plot_data_dict[scar_type][6].append(rt_ins_width)

                # Use the previous plot data to find the y-value of the current bar.
                plot_data_dict[scar_type][7] \
                    .append(previous_y + 0.002 + (0.5 * previous_freq) + y_value)

        plot_data_dict['Marker'] = [(max(marker_list)) * -1, max(marker_list)]
        # sample_name = "{}.{}".format(args.Job_Name, args.SampleName)

        ScarMapperPlot.scarmapperplot(args, datafile=None, sample_name=args.SampleName, plot_data_dict=plot_data_dict,
                                      label_dict=label_dict)

        freq_results_file = \
            open("{}{}_ScarMapper_Combined_Frequency.txt".format(args.WorkingFolder, args.SampleName), "w")

        freq_results_file.write(freq_results_outstring)
        freq_results_file.close()

    warning = "\033[1;31m **See warnings above**\033[m" if log.warning_occurred else ''
    elapsed_time = int(time.time() - start_time)

    log.info("****ScarMapper {0} complete ({1} seconds, {2} Mb peak memory).****"
             .format(module_name, elapsed_time, Tool_Box.peak_memory(), warning))

    exit(0)


def error_checking(options_parser):
    """
    Check parameter file for errors.
    :return:
    :param args:
    """

    args = options_parser.parse_args()

    options_parser.set_defaults(Search_KMER=int(args.Search_KMER))
    options_parser.set_defaults(Spawn=int(args.Spawn))
    options_parser.set_defaults(N_Limit=float(args.N_Limit))
    options_parser.set_defaults(PatternThreshold=float(args.PatternThreshold))
    options_parser.set_defaults(Minimum_Length=int(args.Minimum_Length))
    if getattr(args, "HR_Donor", False):
        options_parser.set_defaults(HR_Donor=args.HR_Donor).upper()

    args = options_parser.parse_args()

    if not pathlib.Path(args.WorkingFolder).exists():
        print("\033[1;31mERROR:\n\tWorking Folder Path: {} Not Found.  Check Options File."
              .format(args.WorkingFolder))
        raise SystemExit(1)

    if getattr(args, "FASTQ1", False) and getattr(args, "ConsensusSequence", False):
        print("\033[1;31mERROR:\n\t--FASTQ1 and --ConsensusSequence both set.  Pick one or the other and try again.")
        raise SystemExit(1)

    if getattr(args, "FASTQ2", False) and getattr(args, "ConsensusSequence", False):
        print("\033[1;31mERROR:\n\t--FASTQ2 and --ConsensusSequence both set.  Pick one or the other and try again.")
        raise SystemExit(1)

    if getattr(args, "FASTQ1", False) and not pathlib.Path(args.FASTQ1).exists():
        print("\033[1;31mERROR:\n\t--FASTQ1: {} Not Found.  Check Options File."
              .format(args.FASTQ1))
        raise SystemExit(1)

    if getattr(args, "FASTQ2", False) and not pathlib.Path(args.FASTQ2).exists():
        print("\033[1;31mERROR:\n\t--FASTQ2: {} Not Found.  Check Options File."
              .format(args.FASTQ2))
        raise SystemExit(1)

    if getattr(args, "ConsensusSequence", False) and not pathlib.Path(args.ConsensusSequence).exists():
        print("\033[1;31mERROR:\n\t--ConsensusSequence: {} Not Found.  Check Options File."
              .format(args.FASTQ2))
        raise SystemExit(1)

    if not getattr(args, "RefSeq", False):
        print("\033[1;31mERROR:\n\t--RefSeq Not Set.  Check Options File.")
        raise SystemExit(1)
    elif "{}".format(os.sep) in args.RefSeq and not os.path.exists(args.RefSeq):
        print("\033[1;31mERROR:\n\t--RefSeq: {} Not Found.  Check Options File."
              .format(args.RefSeq))
        raise SystemExit(1)

    if getattr(args, "Master_Index_File", False) and not pathlib.Path(args.Master_Index_File).exists():
        print("\033[1;31mERROR:\n\t--Master_Index_File: {} Not Found.  Check Options File."
              .format(args.Master_Index_File))
        raise SystemExit(1)

    if getattr(args, "SampleManifest", False) and not pathlib.Path(args.SampleManifest).exists():
        print("\033[1;31mERROR:\n\t--SampleManifest: {} Not Found.  Check Options File."
              .format(args.SampleManifest))
        raise SystemExit(1)

    if getattr(args, "TargetFile", False) and not pathlib.Path(args.TargetFile).exists():
        print("\033[1;31mERROR:\n\t--TargetFile: {} Not Found.  Check Options File."
              .format(args.TargetFile))
        raise SystemExit(1)

    return args


def string_to_boolean(parser):
    """
    Converts strings to boolean.  Done to keep the eval() function out of the code.
    param parser:
    return:
    """
    options_parser = Tool_Box.options_file(parser)
    args = options_parser.parse_args()

    if args.IndelProcessing == "True":
        options_parser.set_defaults(PEAR=True)
        options_parser.set_defaults(Demultiplex=bool(strtobool(args.Demultiplex)))
        options_parser.set_defaults(OutputRawData=bool(strtobool(args.OutputRawData)))
        options_parser.set_defaults(DeleteConsensusFASTQ=bool(strtobool(args.DeleteConsensusFASTQ)))

    options_parser.set_defaults(IndelProcessing=bool(strtobool(args.IndelProcessing)))
    options_parser.set_defaults(Verbose=args.Verbose.upper())
    
    # Add batch mode setting
    if hasattr(args, 'BatchMode'):
        options_parser.set_defaults(BatchMode=bool(strtobool(args.BatchMode)))
    else:
        options_parser.set_defaults(BatchMode=False)

    return options_parser


if __name__ == '__main__':
    main()
