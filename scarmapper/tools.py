"""
Consolidated utilities for ScarMapper.

Combines the subset of Valkyries functions previously spread across Tool_Box,
FASTQ_Tools, and Sequence_Magic that ScarMapper actually uses.
"""

import csv
import getpass
import gzip
import inspect
import logging
import ntpath
import os
import pathlib
import platform
import re
import resource
import socket
import sys
from contextlib import suppress
from datetime import datetime

import magic
from Levenshtein import distance


# ---------------------------------------------------------------------------
# General utilities
# ---------------------------------------------------------------------------

def peak_memory():
    """Return peak memory usage in MB."""
    peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    if platform.system() == 'darwin':
        peak_mb /= 1024
    return int(peak_mb)


def delete(file_list):
    """Delete one or more files, ignoring missing files."""
    for file in file_list:
        if os.path.isfile(file):
            with suppress(FileNotFoundError):
                os.remove(str(file))


def compress_files(file, log):
    """Gzip-compress a file in place (max compression)."""
    if os.path.isfile(file):
        delete([file + ".gz"])
    cmd = "gzip -9 " + file
    os.system(cmd)
    log.debug(f"{file} Compressed")


def debug_messenger(reason=None):
    """Print filename and line number for debugging."""
    if reason is None:
        reason = "Programmer Neglected to Enlighten Us About the Need for Debugging This Section."
    frameinfo = inspect.getframeinfo(inspect.currentframe().f_back)
    print(f"\033[1;31m***WARNING: Debugging Module {frameinfo.filename} at Line {frameinfo.lineno}.\n\t-->REASON: {reason}\033[m")


# ---------------------------------------------------------------------------
# Options-file parser
# ---------------------------------------------------------------------------

def options_file(options_parser):
    """
    Parse a tab-delimited options file and add each key-value pair to an
    argparse parser as a default.
    """
    count = 0
    config_file = options_parser.parse_args().options_file

    if not os.path.isfile(config_file):
        print(f"\033[1;31mWARNING:\n\tOptions_File {config_file} Not Found.  Check File Name and Path.")
        raise SystemExit(1)

    options = csv.reader(open(config_file), delimiter='\t')

    for line in options:
        count += 1
        if len(line) < 2 or "#" in str(line[0]):
            continue
        try:
            value = re.sub(r"[\s]", "", line[1].split("#")[0])
        except IndexError:
            raise SystemExit(f"There is a syntax error in the options file on line {count}")

        key = line[0].strip('--')
        options_parser.add_argument(line[0], dest=key, default=value)

    return options_parser


# ---------------------------------------------------------------------------
# File parsing
# ---------------------------------------------------------------------------

class FileParser:
    @staticmethod
    def indices(log, input_file):
        """
        Parse a tab-delimited index/manifest/target file and return a list of
        rows.  Blank lines and comment lines (starting with ``#``) are skipped.
        """
        log.info(f"Parsing {input_file}")
        if not os.path.isfile(input_file):
            log.error(f"{input_file} Not Found.  Check File Name and Path.")
            raise SystemExit(1)

        index_list = []
        line_num = 0
        index_file = list(csv.reader(open(input_file), delimiter='\t'))

        for line in index_file:
            line_num += 1
            col_count = len(line)

            if col_count > 0 and len(line[0].split("#")[0]) > 0:
                tmp_line = []
                for i in range(col_count):
                    try:
                        line[i] = line[i].split("#")[0]
                    except IndexError:
                        raise SystemExit(
                            f"There is a syntax error in file {input_file} on line {line_num}, column {i} ")
                    line[i] = re.sub(",", '', line[i])
                    tmp_line.append(line[i])
                index_list.append(tmp_line)

        log.debug(f"Parsing Complete for  {input_file}")
        return index_list


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

class Logger:
    _DATE_FORMAT = '%Y-%m-%d %H:%M:%S'
    _FILE_LOG_FORMAT = '%(asctime)s|%(levelname)s|%(start_time)s|%(host)s|%(user)s|%(message)s'
    _CONSOLE_LOG_FORMAT = '%(asctime)s|%(levelname)s|%(message)s'

    def __init__(self, args, console_stream=None, parellel_id=None):
        self._verbose = args.Verbose
        if parellel_id:
            log_file = f"{args.Job_Name}_{parellel_id}"
        else:
            log_file = args.Job_Name

        self._log_filename = f"{args.WorkingFolder}{log_file}.log"

        try:
            log = open(self._log_filename, "w")
            log.close()
        except IOError:
            raise SystemExit(
                f'Cannot create log file [{self._log_filename}]. Review inputs and try again.')

        if console_stream:
            self._console_stream = console_stream
        else:
            self._console_stream = sys.stderr

        user = getpass.getuser()
        host = socket.gethostname()
        start_time = datetime.now().strftime(Logger._DATE_FORMAT)

        self._logging_dict = {'user': user, 'host': host, 'start_time': start_time}

        logging.basicConfig(
            format=Logger._FILE_LOG_FORMAT,
            level=args.Verbose,
            datefmt=Logger._DATE_FORMAT,
            filename=self._log_filename,
        )
        self._file_logger = logging
        self.warning_occurred = False

    def _print(self, level, message, args):
        now = datetime.now().strftime(Logger._DATE_FORMAT)
        print(
            Logger._CONSOLE_LOG_FORMAT % {
                'asctime': now, 'levelname': level,
                'message': self._format(message, args),
            },
            file=self._console_stream,
        )
        self._console_stream.flush()

    @staticmethod
    def _format(message, args):
        try:
            log_message = message.format(*[i for i in args])
        except IndexError as err:
            log_message = (
                "Malformed log message ({}: {})""|{}|{}"
                .format(type(err).__name__, err, message, [str(i) for i in args])
            )
        return log_message

    def debug(self, message, *args):
        if self._verbose == "DEBUG":
            self._print("\033[96mDEBUG\033[m", message, args)
        self._file_logger.debug(self._format(message, args), extra=self._logging_dict)

    def _log(self, msg_type, method, message, *args):
        self._print(msg_type, message, args)
        method(self._format(message, args), extra=self._logging_dict)

    def error(self, message, *args):
        self._log("\033[38;5;202mERROR\033[m", self._file_logger.error, message, *args)

    def info(self, message, *args):
        self._log("\033[38;5;220mINFO\033[m", self._file_logger.info, message, *args)

    def warning(self, message, *args):
        self._log("\033[1;31mWARNING\033[m", self._file_logger.warning, message, *args)
        self.warning_occurred = True


def log_environment_info(log, args, command_line_args):
    log.info(f"original_command_line|{' '.join(command_line_args)}")
    log.info(f'command_options|{args}')
    log.info(f'WorkingFolder|{args.WorkingFolder}')
    log.info(f'platform_uname|{platform.uname()}')
    log.info(f'platform_python_version|{platform.python_version()}')


# ---------------------------------------------------------------------------
# FASTQ I/O
# ---------------------------------------------------------------------------

class FASTQ_Reader:
    """Read a FASTQ file (plain text or gzipped) one record at a time via a generator."""

    __slots__ = ['input_file', 'log', 'name', 'seq', 'index', 'qual', 'read_block', 'file_name', 'fq_file']

    def __init__(self, input_file, log=None):
        self.name = None
        self.seq = None
        self.index = None
        self.qual = None
        self.input_file = input_file
        self.log = log
        self.read_block = []
        self.file_name = ntpath.basename(input_file)
        self.fq_file = self._open_fastq()

    def _open_fastq(self):
        if len(self.input_file) < 3:
            self.log.warning("FASTQ file parameter missing from options file. Correct error and try again.")
            raise SystemExit(1)

        if not pathlib.Path(self.input_file).is_file():
            self.log.warning(f"FASTQ file {self.input_file} not found.  Correct error and run again.")
            raise SystemExit(1)

        try:
            mime_type = magic.from_file(self.input_file, mime=True).decode()
        except AttributeError:
            mime_type = magic.from_file(self.input_file, mime=True)

        if "text" in mime_type:
            return open(self.input_file, 'rU')
        elif "gzip" in mime_type:
            return gzip.open(self.input_file, 'rt', encoding='utf-8')
        else:
            self.log.warning(f"Unsupported file-type for {self.input_file}.  Only TEXT or GZIP Allowed.")
            raise SystemExit(1)

    def line_reader(self):
        for line in self.fq_file:
            while True:
                yield line

    def seq_read(self):
        """Generator that yields one FASTQ record at a time."""
        read_block = []
        count = 0
        eof = False
        try:
            while count < 4:
                read_block.append(next(FASTQ_Reader.line_reader(self)))
                count += 1
        except StopIteration:
            eof = True

        if len(read_block) == 4 and not eof:
            self.name = read_block[0].strip("\n").strip("@")
            self.seq = read_block[1].strip("\n").strip()
            self.index = read_block[2].strip("\n").strip()
            self.qual = read_block[3].strip("\n").strip()

            if len(self.seq) != len(self.qual):
                self.log.error(
                    f"Sequence and quality scores of different lengths! \n{self.name}\n{self.seq}\n{self.index}\n{self.qual}")
                raise ValueError(
                    f"Sequence and quality scores of different lengths! \n{self.name}\n{self.seq}\n{self.index}\n{self.qual}")
            yield self

        self.name = None


class Writer:
    """Write FASTQ records to a file."""

    __slots__ = ['log', 'file']

    def __init__(self, log, out_file_string):
        self.file = open(out_file_string, "w")
        self.log = log

    def write(self, read_list):
        outstring = ""
        for read in read_list:
            try:
                assert len(read[1]) == len(read[2])
            except AssertionError:
                self.log.error(
                    f"Sequence and quality scores of different lengths! Read Name {read[0]}; Seq Length {len(read[1])}; "
                    f"Qual Length {len(read[2])}")
                raise SystemExit(1)
            outstring += f"@{read[0]}\n{read[1]}\n+\n{read[2]}\n"

        self.file.write(outstring)
        read_list.clear()
        return True

    def close(self):
        self.file.close()
        return True


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

_DNA_COMPLEMENT = str.maketrans(
    "ACGTMRWSYKVHDBXNacgtmrwsykvhdbxn",
    "TGCAKYWSRMBDHVXNtgcakywsrmbdhvxn",
)


def rcomp(seq):
    """Return the reverse complement of a DNA sequence."""
    return seq[::-1].translate(_DNA_COMPLEMENT)


def match_maker(query, unknown):
    """
    Compute the Levenshtein distance between *query* and *unknown*, adjusted
    for any length difference.
    """
    query_mismatch = distance(query, unknown)
    adjusted = query_mismatch - (len(unknown) - len(query))
    return adjusted
