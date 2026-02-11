"""
Combine replicate frequency files and generate merged output with plots.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2023
"""

import collections
import csv
import glob

from natsort import natsort
from scipy import stats
from scipy.stats import gmean

from Valkyries import Tool_Box
from scarmapper import ScarMapperPlot

__author__ = 'Dennis A. Simpson'
__version__ = '2.0.0'

# Minimum frequency for a scar pattern to appear in the plot
_PLOT_FREQ_FLOOR = 0.0025


def combine_frequency_files(args, log, version, run_start):
    """
    Read all ``*ScarMapper_Frequency.txt`` files from ``args.DataFiles``,
    merge replicates via geometric mean, write a combined frequency file,
    and produce a summary plot.
    """
    file_list = sorted(glob.glob(f"{args.DataFiles}*ScarMapper_Frequency.txt"))
    file_count = len(file_list)
    if file_count == 0:
        log.error(f"No frequency files found in {args.DataFiles}")
        raise SystemExit(1)

    log.info(f"Found {file_count} frequency files to merge")

    page_header = _build_page_header(file_list[0], version, run_start, args.SampleName)
    data_dict = _collect_replicate_data(file_list, log)
    _write_combined_output(args, log, data_dict, file_count, page_header)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_page_header(first_file, version, run_start, sample_name):
    """Extract run metadata from the first frequency file to form a header."""
    header = f"# ScarMapper File Merge v{version}\n# Run: {run_start}\n# Sample Name: {sample_name}\n"

    with open(first_file) as fh:
        for line_num, line in enumerate(fh):
            if not line.strip():
                break
            if line_num > 3:
                header += f"{line.split(chr(9))[0]}\n"

    return header + "\n\n"


def _collect_replicate_data(file_list, log):
    """
    Walk every frequency file and group rows by a composite key
    (scar_type | lft_del | insertion | insertion_size).

    Returns a dict mapping that key to ``[[freq1, freq2, …], row_data]``.
    """
    data_dict = collections.defaultdict(list)

    for file_name in file_list:
        freq_data = Tool_Box.FileParser.indices(log, file_name)
        for row in freq_data:
            key = f"{row[3]}|{row[4]}|{row[6]}|{row[8]}"
            row_data = row[2:]

            if key in data_dict:
                data_dict[key][0].append(float(row[1]))
            else:
                data_dict[key] = [[float(row[1])], row_data]

    return data_dict


def _write_combined_output(args, log, data_dict, file_count, page_header):
    """
    From the grouped replicate data, compute geometric means, build the
    combined frequency file and a plot.
    """
    plot_data = collections.defaultdict(list)
    label_dict = collections.defaultdict(float)
    output_data = {}
    marker_list = []

    for key, (freq_list, row_data) in data_dict.items():
        # Require the pattern to appear in at least half of replicates
        if len(freq_list) / file_count < 0.5:
            continue

        freq = gmean(freq_list)
        sem = stats.sem(freq_list)
        row_string = "\t".join(row_data)

        # Use frequency as output key; nudge on the rare collision
        output_key = freq
        while output_key in output_data:
            output_key += 1e-16

        scar_type = row_data[0]
        label_dict[scar_type] += freq

        lft_del = int(row_data[1])
        rt_del = int(row_data[2])
        mh_size = int(row_data[5])
        ins_size = int(row_data[7])

        output_data[output_key] = (
            (freq, lft_del, rt_del, mh_size, ins_size, scar_type),
            f"{freq}\t{sem}\t{row_string}\n",
        )

    # Assemble output string
    out_string = (
        f"{page_header}"
        "# Frequency\tSEM\tScar Type\tLeft Deletions\tRight Deletions\tDeletion Size\t"
        "Microhomology\tMicrohomology Size\tInsertion\tInsertion Size\t"
        "Left Template\tRight Template\tConsensus Left Junction\t"
        "Consensus Right Junction\tTarget Left Junction\tTarget Right Junction\t"
        "Consensus\tTarget Region\n"
    )

    for k in natsort.natsorted(output_data, reverse=True):
        metrics, line = output_data[k]
        out_string += line

        freq, lft_del, rt_del, mh_size, ins_size, scar_type = metrics

        if freq < _PLOT_FREQ_FLOOR:
            continue

        _accumulate_plot_data(plot_data, marker_list, freq, lft_del, rt_del, mh_size, ins_size, scar_type)

    # Finalise marker range for x-axis limits
    if marker_list:
        plot_data['Marker'] = [max(marker_list) * -1, max(marker_list)]

    ScarMapperPlot.scarmapperplot(
        args, datafile=None, sample_name=args.SampleName,
        plot_data_dict=plot_data, label_dict=label_dict,
    )

    out_path = f"{args.WorkingFolder}{args.SampleName}_ScarMapper_Combined_Frequency.txt"
    with open(out_path, "w") as fh:
        fh.write(out_string)

    log.info(f"Combined frequency file written to {out_path}")


def _accumulate_plot_data(plot_data, marker_list, freq, lft_del, rt_del, mh_size, ins_size, scar_type):
    """
    Append one scar-pattern's metrics into *plot_data* (a dict of lists),
    computing the stacked y-position along the way.
    """
    y_value = freq * 0.5
    lft_ins_width = freq
    rt_ins_width = freq

    marker_list.extend([lft_del + mh_size * 0.5, rt_del + mh_size * 0.5, ins_size])

    lft_del_val = (lft_del + mh_size * 0.5) * -1
    rt_del_val = rt_del + mh_size * 0.5
    lft_ins_val = (ins_size * 0.5) * -1
    rt_ins_val = ins_size * 0.5

    if lft_del + mh_size * 0.5 != 0:
        lft_ins_width = freq * 0.5
    if rt_del + mh_size * 0.5 != 0:
        rt_ins_width = freq * 0.5

    if scar_type not in plot_data:
        plot_data[scar_type] = [
            [freq], [lft_del_val], [rt_del_val], [lft_ins_val],
            [rt_ins_val], [lft_ins_width], [rt_ins_width], [y_value],
        ]
    else:
        prev = plot_data[scar_type]
        n = len(prev[0])
        prev_freq = prev[0][n - 1]
        prev_y = prev[7][n - 1]

        prev[0].append(freq)
        prev[1].append(lft_del_val)
        prev[2].append(rt_del_val)
        prev[3].append(lft_ins_val)
        prev[4].append(rt_ins_val)
        prev[5].append(lft_ins_width)
        prev[6].append(rt_ins_width)
        prev[7].append(prev_y + 0.002 + 0.5 * prev_freq + y_value)
