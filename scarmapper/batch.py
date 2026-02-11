"""
Batch processing for pre-demultiplexed samples.

Each sample has its own FASTQ file(s) and is processed independently,
with results collected into a single summary file.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2023
"""

import datetime
import itertools
import pathlib

import pathos

from Valkyries import Tool_Box, FASTQ_Tools
from scarmapper import INDEL_Processing as Indel_Processing, TargetMapper as Target_Mapper
from scarmapper.pear import pear_consensus

__author__ = 'Dennis A. Simpson'
__version__ = '2.0.0'


# ---------------------------------------------------------------------------
# Manifest helpers
# ---------------------------------------------------------------------------

# Column indices in the sample manifest
_COL = dict(
    index_name=0, sample_name=1, replicate=2, fwd_phase=3, rev_phase=4,
    locus=5, fastq1=6, fastq2=7, hr_donor=8,
)


def _parse_manifest_row(row):
    """
    Turn a single manifest row into a normalised dictionary.
    Missing optional columns default to empty strings.
    """

    def _col(key):
        idx = _COL[key]
        return row[idx].strip() if len(row) > idx and row[idx] else ''

    return {
        'Index_Name': _col('index_name'),
        'Sample_Name': _col('sample_name'),
        'Replicate': _col('replicate'),
        'Forward_Phase': _col('fwd_phase'),
        'Reverse_Phase': _col('rev_phase'),
        'Locus': _col('locus'),
        'FASTQ1_Path': _col('fastq1'),
        'FASTQ2_Path': _col('fastq2'),
        'HR_Donor': _col('hr_donor').upper(),
    }


# ---------------------------------------------------------------------------
# Per-sample processing
# ---------------------------------------------------------------------------

def _process_single_sample(args, log, sample_info, target_mapper, version, run_start):
    """
    Process one pre-demultiplexed sample through PEAR → ScarSearch.

    Returns the ``summary_data`` list produced by ScarSearch, or None on failure.
    """
    index_name = sample_info['Index_Name']
    fastq1 = sample_info['FASTQ1_Path']
    fastq2 = sample_info.get('FASTQ2_Path', '')
    sample_hr = sample_info.get('HR_Donor', '')

    log.info(f"Processing sample: {index_name}")
    if sample_hr:
        log.info(f"  Sample-specific HR_Donor: {sample_hr}")

    is_paired = fastq2 and pathlib.Path(fastq2).exists()
    file_list = []

    # ------------------------------------------------------------------
    # Read input (paired-end via PEAR, or single-end directly)
    # ------------------------------------------------------------------
    if is_paired:
        log.info(f"  Paired-end reads detected. Running PEAR consensus.")
        file_list = pear_consensus(args, log, fq1=fastq1, fq2=fastq2, sample_prefix=index_name)
        if not file_list:
            log.error(f"  PEAR failed for {index_name}. Skipping.")
            return None
        fq1 = FASTQ_Tools.FASTQ_Reader(file_list[0], log)
    else:
        log.info(f"  Single-end reads detected. Skipping PEAR.")
        fq1 = FASTQ_Tools.FASTQ_Reader(fastq1, log)

    # ------------------------------------------------------------------
    # Read sequences into memory
    # ------------------------------------------------------------------
    sequence_list = []
    eof = False
    while not eof:
        try:
            read = next(fq1.seq_read())
            sequence_list.append(read.seq)
        except StopIteration:
            eof = True

    read_count = len(sequence_list)
    log.info(f"  Total sequences for {index_name}: {read_count}")

    # ------------------------------------------------------------------
    # Build a minimal index dict for ScarSearch
    # ------------------------------------------------------------------
    index_dict = {
        index_name: [
            '',                            # right_index_sequence (unused)
            0,                             # right_index_count
            '',                            # left_index_sequence  (unused)
            0,                             # left_index_count
            index_name,
            sample_info['Sample_Name'],
            sample_info['Replicate'],
            sample_info['Locus'],
            sample_hr,                     # sample-specific HR_Donor
        ]
    }

    # ------------------------------------------------------------------
    # Run ScarSearch
    # ------------------------------------------------------------------
    scar_search = Indel_Processing.ScarSearch(
        log, args, version, run_start,
        target_mapper.targets,
        index_dict,
        index_name,
        sequence_list,
        read_count,
        args.PatternThreshold * 0.00001,
    )

    # ------------------------------------------------------------------
    # Cleanup PEAR artefacts
    # ------------------------------------------------------------------
    if is_paired and file_list:
        if args.DeleteConsensusFASTQ:
            log.info(f"  Deleting PEAR FASTQ files for {index_name}.")
            Tool_Box.delete(file_list)
        else:
            log.info(f"  Compressing PEAR FASTQ files for {index_name}.")
            pool = pathos.multiprocessing.Pool(args.Spawn)
            pool.starmap(Tool_Box.compress_files, zip(file_list, itertools.repeat(log)))

    return scar_search.summary_data


# ---------------------------------------------------------------------------
# Top-level batch entry point
# ---------------------------------------------------------------------------

def batch_process_samples(args, log, version, run_start):
    """
    Iterate over all samples listed in the manifest, process each, and
    write a combined summary file.
    """
    log.info("Starting batch processing of pre-demultiplexed samples")

    sample_manifest = Tool_Box.FileParser.indices(log, args.SampleManifest)
    target_mapper = Target_Mapper.TargetMapper(log, args, sample_manifest)

    # Validate and collect samples
    samples = []
    for row in sample_manifest:
        info = _parse_manifest_row(row)
        fq1 = info['FASTQ1_Path']
        if not fq1 or not pathlib.Path(fq1).exists():
            log.error(f"FASTQ1 not found for {info['Index_Name']}: {fq1}")
            continue
        samples.append(info)

    log.info(f"Found {len(samples)} valid samples")

    # Process (parallel or serial)
    if args.Spawn > 1:
        log.info(f"Processing in parallel ({args.Spawn} workers)")
        pool = pathos.multiprocessing.Pool(args.Spawn)
        all_summaries = pool.starmap(
            _process_single_sample,
            [(args, log, s, target_mapper, version, run_start) for s in samples],
        )
    else:
        log.info("Processing serially")
        all_summaries = [
            _process_single_sample(args, log, s, target_mapper, version, run_start)
            for s in samples
        ]

    all_summaries = [s for s in all_summaries if s is not None]
    _write_batch_summary(args, log, all_summaries, version, run_start, samples)
    log.info("Batch processing complete")


# ---------------------------------------------------------------------------
# Summary output
# ---------------------------------------------------------------------------

def _write_batch_summary(args, log, all_summaries, version, run_start, samples):
    """Write the combined batch-summary TSV."""
    run_stop = datetime.datetime.today().strftime("%a %b %d %H:%M:%S %Y")
    has_hr = any(s.get('HR_Donor') for s in samples)

    # Build lookup tables
    hr_lookup = {s['Index_Name']: s.get('HR_Donor', '') for s in samples}
    name_lookup = {s['Index_Name']: s.get('Sample_Name', s['Index_Name']) for s in samples}
    rep_lookup = {s['Index_Name']: s.get('Replicate', '1') for s in samples}

    lines = [
        f"ScarMapper Batch Processing v{version}",
        f"Start: {run_start}",
        f"End: {run_stop}",
        "",
    ]

    # Header row
    header_cols = [
        "Index_Name", "Sample_Name", "Replicate", "Target",
        "Total_Reads", "Passing_Filters", "Scar_Count", "Scar_Fraction",
    ]
    if has_hr:
        header_cols += ["HR_Donor", "HR_Count", "HR_Fraction"]
    header_cols += ["TMEJ", "NHEJ", "Non-MH_Deletion", "Insertion", "SNV"]
    lines.append("\t".join(header_cols))

    # Data rows
    for sd in all_summaries:
        if not sd:
            continue

        idx = sd[0]
        passing = sd[1]
        target = sd[9]
        jdata = sd[8]

        scar_count = passing - sd[6][0] - sd[6][1]
        scar_frac = scar_count / passing if passing > 0 else 0

        row = [
            idx,
            name_lookup.get(idx, idx),
            rep_lookup.get(idx, '1'),
            target,
            str(passing),
            str(passing),
            str(scar_count),
            f"{scar_frac:.4f}",
        ]

        if has_hr:
            donor = hr_lookup.get(idx, '')
            hr_total = (sd[10][0] + sd[10][1]) if len(sd) > 10 else 0
            hr_frac = hr_total / passing if passing > 0 else 0
            row += [donor, str(hr_total), f"{hr_frac:.4f}"]

        row += [
            str(jdata[0] if jdata else 0),   # TMEJ
            str(jdata[1] if jdata else 0),   # NHEJ
            str(jdata[4] if jdata else 0),   # Non-MH Deletion
            str(jdata[2] if jdata else 0),   # Insertion
            str(jdata[5] if jdata else 0),   # SNV
        ]
        lines.append("\t".join(row))

    out_path = f"{args.WorkingFolder}{args.Job_Name}_Batch_Summary.txt"
    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    log.info(f"Batch summary written to {out_path}")
