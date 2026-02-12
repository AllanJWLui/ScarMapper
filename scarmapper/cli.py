"""CLI entry point, argument parsing, validation, and top-level dispatch."""

import argparse
import datetime
import itertools
import os
import pathlib
import sys
import time

import pathos

from scarmapper import INDEL_Processing as Indel_Processing, TargetMapper as Target_Mapper, tools
from scarmapper.pear import pear_consensus

__version__ = '3.0.0'
__package__ = 'ScarMapper'


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _strtobool(val: str) -> bool:
    """
    Replace deprecated distutils.util.strtobool.
    Accepts 'true', '1', 'yes' (case-insensitive) as True; everything else as False.
    """
    return val.strip().lower() in ('true', '1', 'yes')


def _parse_args(command_line_args=None):
    """
    Build the base ArgumentParser and layer options-file values on top.
    Returns the fully resolved *options_parser* (not yet finalized with defaults).
    """
    parser = argparse.ArgumentParser(
        description=f"A package to map genomic repair scars at defined loci.\n {__package__} v{__version__}",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        '--options_file', action='store', dest='options_file', required=True,
        help='File containing program parameters.',
    )

    options_parser = tools.options_file(parser)
    return options_parser


def _apply_boolean_defaults(options_parser):
    """
    Convert string options to proper Python types and set sensible defaults.
    """
    args = options_parser.parse_args()

    if args.IndelProcessing == "True":
        options_parser.set_defaults(PEAR=True)
        options_parser.set_defaults(Demultiplex=_strtobool(args.Demultiplex))
        options_parser.set_defaults(OutputRawData=_strtobool(args.OutputRawData))
        options_parser.set_defaults(DeleteConsensusFASTQ=_strtobool(args.DeleteConsensusFASTQ))

    options_parser.set_defaults(IndelProcessing=_strtobool(args.IndelProcessing))
    options_parser.set_defaults(Verbose=args.Verbose.upper())
    options_parser.set_defaults(BatchMode=_strtobool(args.BatchMode) if hasattr(args, 'BatchMode') else False)

    return options_parser


def _apply_numeric_defaults(options_parser):
    """
    Cast numeric options from strings and apply an upper-case HR_Donor if present.
    """
    args = options_parser.parse_args()

    options_parser.set_defaults(
        Search_KMER=int(args.Search_KMER),
        Spawn=int(args.Spawn),
        N_Limit=float(args.N_Limit),
        PatternThreshold=float(args.PatternThreshold),
        Minimum_Length=int(args.Minimum_Length),
    )

    if getattr(args, "HR_Donor", False):
        options_parser.set_defaults(HR_Donor=args.HR_Donor.upper())

    return options_parser


def _validate_paths(args):
    """
    Verify that every user-supplied path actually exists.
    Raises SystemExit with a descriptive message on failure.
    """
    _require_dir(args.WorkingFolder, "Working Folder")

    # Mutually exclusive: raw FASTQ vs. pre-built consensus
    if getattr(args, "FASTQ1", False) and getattr(args, "ConsensusSequence", False):
        _die("--FASTQ1 and --ConsensusSequence both set. Pick one or the other.")
    if getattr(args, "FASTQ2", False) and getattr(args, "ConsensusSequence", False):
        _die("--FASTQ2 and --ConsensusSequence both set. Pick one or the other.")

    _require_file_if_set(args, "FASTQ1")
    _require_file_if_set(args, "FASTQ2")
    _require_file_if_set(args, "ConsensusSequence")
    _require_file_if_set(args, "Master_Index_File")
    _require_file_if_set(args, "SampleManifest")
    _require_file_if_set(args, "TargetFile")

    if not getattr(args, "RefSeq", False):
        _die("--RefSeq not set. Check Options File.")
    elif os.sep in args.RefSeq and not os.path.exists(args.RefSeq):
        _die(f"--RefSeq: {args.RefSeq} not found. Check Options File.")


def _require_dir(path, label):
    if not pathlib.Path(path).exists():
        _die(f"{label} path: {path} not found. Check Options File.")


def _require_file_if_set(args, attr):
    value = getattr(args, attr, False)
    if value and not pathlib.Path(value).exists():
        _die(f"--{attr}: {value} not found. Check Options File.")


def _die(message):
    print(f"\033[1;31mERROR:\n\t{message}\033[m")
    raise SystemExit(1)


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

def main(command_line_args=None):
    """
    Top-level entry point.
    """
    start_time = time.time()

    if not command_line_args:
        command_line_args = sys.argv

    run_start = datetime.datetime.today().strftime("%H:%M:%S %Y  %a %b %d")

    # Parse and validate
    options_parser = _parse_args(command_line_args)
    options_parser = _apply_boolean_defaults(options_parser)
    options_parser = _apply_numeric_defaults(options_parser)
    args = options_parser.parse_args()
    _validate_paths(args)

    log = tools.Logger(args)
    tools.log_environment_info(log, args, command_line_args)
    log.info(f"{__package__} v{__version__}")

    # Dispatch to the appropriate processing mode
    if args.BatchMode:
        _run_batch(args, log, run_start)
    elif args.IndelProcessing:
        _run_indel_processing(args, log, run_start)
    else:
        _run_combine(args, log, run_start)

    warning = "\033[1;31m **See warnings above**\033[m" if log.warning_occurred else ''
    elapsed = int(time.time() - start_time)
    log.info(f"****ScarMapper complete ({elapsed} seconds, {tools.peak_memory()} Mb peak memory).****{warning}")


def _run_batch(args, log, run_start):
    """Dispatch: batch mode for pre-demultiplexed samples."""
    from scarmapper.batch import batch_process_samples

    log.info("Running in Batch Mode for pre-demultiplexed samples")
    batch_process_samples(args, log, __version__, run_start)


def _run_indel_processing(args, log, run_start):
    """Dispatch: standard indel-processing pipeline (multiplexed FASTQ)."""
    valid_platforms = ("Illumina", "TruSeq", "Ramsden")
    if args.Platform not in valid_platforms:
        log.error(f"--Platform must be one of {valid_platforms}.")
        raise SystemExit(1)

    log.info("Sending FASTQ files to FASTQ preprocessor.")

    file_list = []
    if args.PEAR:
        file_list = pear_consensus(args, log)
        if not file_list:
            log.error("PEAR failed. Check logs.")
            raise SystemExit(1)
        fq1 = tools.FASTQ_Reader(file_list[0], log)
        fq2 = None
    else:
        fq1 = tools.FASTQ_Reader(args.FASTQ1, log)
        fq2 = tools.FASTQ_Reader(args.FASTQ2, log)

    sample_manifest = tools.FileParser.indices(log, args.SampleManifest)
    targeting = Target_Mapper.TargetMapper(log, args, sample_manifest)

    processor = Indel_Processing.DataProcessing(
        log, args, run_start, __version__, targeting, fq1, fq2,
    )
    processor.main_loop()

    # Compress or delete PEAR files
    if args.PEAR and file_list:
        if args.DeleteConsensusFASTQ:
            log.info("Deleting PEAR FASTQ files.")
            tools.delete(file_list)
        else:
            log.info(f"Compressing {len(file_list)} FASTQ files generated by PEAR.")
            pool = pathos.multiprocessing.Pool(args.Spawn)
            pool.starmap(tools.compress_files, zip(file_list, itertools.repeat(log)))


def _run_combine(args, log, run_start):
    """Dispatch: combine replicate frequency files and plot."""
    from scarmapper.combine import combine_frequency_files

    log.info("Process Replicates.")
    combine_frequency_files(args, log, __version__, run_start)
