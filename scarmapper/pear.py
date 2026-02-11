"""PEAR consensus generation and FASTQ preprocessing utilities."""

import os
import pathlib
import subprocess

from scarmapper import tools


def pear_consensus(args, log, fq1=None, fq2=None, sample_prefix=None):
    """
    Run PEAR to merge paired-end reads into consensus sequences.

    Parameters
    ----------
    args : argparse.Namespace
        Program arguments.
    log : tools.Logger
        Logger instance.
    fq1 : str, optional
        Path to R1 FASTQ (overrides args.FASTQ1 for batch mode).
    fq2 : str, optional
        Path to R2 FASTQ (overrides args.FASTQ2 for batch mode).
    sample_prefix : str, optional
        Output-file prefix (overrides args.Job_Name for batch mode).

    Returns
    -------
    list[str] or None
        Paths to generated FASTQ files, or None on failure.
    """
    log.info("Beginning PEAR Consensus")

    fastq1 = fq1 or args.FASTQ1
    fastq2 = fq2 or args.FASTQ2
    prefix = f"{args.WorkingFolder}{sample_prefix or args.Job_Name}"

    output_files = {
        'assembled': f"{prefix}.assembled.fastq",
        'discarded': f"{prefix}.discarded.fastq",
        'unassembled_fwd': f"{prefix}.unassembled.forward.fastq",
        'unassembled_rev': f"{prefix}.unassembled.reverse.fastq",
    }

    # Build the PEAR command from optional parameters
    pear_binary = f"{pathlib.Path(__file__).parent.parent.absolute()}{os.sep}Pear{os.sep}bin{os.sep}pear"
    cmd_parts = [
        pear_binary,
        f"-f {fastq1}",
        f"-r {fastq2}",
        f"-o {prefix}",
        f"-y {args.Memory}",
        f"-j {args.Spawn}",
    ]

    _optional = [
        ('PValue', '-p'),
        ('MinOverlap', '-v'),
        ('QualityThreshold', '-q'),
        ('PhredValue', '-b'),
        ('TestMethod', '-g'),
        ('MinConsensusLength', '-n'),
    ]
    for attr, flag in _optional:
        value = getattr(args, attr, None)
        if value:
            cmd_parts.append(f"{flag} {value}")

    proc = subprocess.run(
        " ".join(cmd_parts),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
    )

    if proc.stderr:
        log.error(f"{proc.stderr.decode()}\n{proc.stdout.decode()}\n")
        return None

    log.info(
        "Begin PEAR Output\n"
        f"{'-' * 106}\n{proc.stdout.decode()}\n{'-' * 106}\n"
    )

    file_list = [
        output_files['assembled'],
        output_files['unassembled_fwd'],
        output_files['unassembled_rev'],
    ]

    discarded = output_files['discarded']
    if pathlib.Path(discarded).exists():
        file_list.append(discarded)
    else:
        tools.delete([discarded])

    return file_list


def preprocess_misoriented_reads(args, log):
    """
    Re-orient reads where R1 contains a mix of forward and reverse orientations.

    This was created for a specific sequencing run where reads were not properly
    oriented off the sequencer.  Forward reads in R1 are identified by the presence
    of ``CTGGCTCCA``; reverse reads by ``TGTGGCTCTG``.

    Parameters
    ----------
    args : argparse.Namespace
    log : tools.Logger

    Returns
    -------
    tuple[str, str]
        Paths to the corrected R1 and R2 FASTQ files.
    """
    # TODO: Make the orientation-detection sequences configurable via the options file.
    FORWARD_MARKER = "CTGGCTCCA"
    REVERSE_MARKER = "TGTGGCTCTG"

    fq1 = tools.FASTQ_Reader(args.FASTQ1, log)
    fq2 = tools.FASTQ_Reader(args.FASTQ2, log)

    corrected_r1_path = f"{args.WorkingFolder}{args.Job_Name}_corrected_R1.fastq"
    corrected_r2_path = f"{args.WorkingFolder}{args.Job_Name}_corrected_R2.fastq"

    # Remove stale output if it exists
    tools.delete([corrected_r1_path, corrected_r2_path])

    read_counter = 0
    good_reads = 0
    filtered_count = 0
    buf1 = ""
    buf2 = ""

    with open(corrected_r1_path, "a") as out1, open(corrected_r2_path, "a") as out2:
        eof = False
        while not eof:
            try:
                read_counter += 1
                r1 = next(fq1.seq_read())
                r2 = next(fq2.seq_read())
            except StopIteration:
                eof = True
                continue

            # Filter unbalanced or N-heavy read pairs
            if len(r1.seq) != len(r2.seq) or r1.seq.count("N") > 5 or r2.seq.count("N") > 5:
                filtered_count += 1
                continue

            r1_line = f"@{r1.name}\n{r1.seq}\n+\n{r1.qual}\n"
            r2_line = f"@{r2.name}\n{r2.seq}\n+\n{r2.qual}\n"

            if FORWARD_MARKER in r1.seq:
                buf1 += r1_line
                buf2 += r2_line
                good_reads += 1
            elif REVERSE_MARKER in r1.seq:
                # Swap: R1 is actually the reverse read
                buf1 += r2_line
                buf2 += r1_line
                good_reads += 1

            if read_counter % 1_000_000 == 0:
                log.info(f"Processed {read_counter} reads. {filtered_count} filtered, {good_reads} accepted.")
                out1.write(buf1)
                out2.write(buf2)
                buf1 = ""
                buf2 = ""

        # Flush remaining
        out1.write(buf1)
        out2.write(buf2)

    log.info(f"{filtered_count} filtered. {good_reads} accepted. {read_counter} total.")
    return corrected_r1_path, corrected_r2_path
