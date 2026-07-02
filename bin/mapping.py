#!/usr/bin/env python

import argparse
import pandas as pd
import mappy as mp
import math
from collections import defaultdict
import os
import plotly.graph_objects as go
import plotly.express as px
from plotly.io import templates
import sys

templates.default = "plotly_white"


def map_to_refs(query: str, reference: str, preset=None) -> tuple[int, dict, list]:
    """Map reads to specified reference set and use
    mappy to calculate mapping stats
    Parameters
        query (str): path to reads
        reference (str): path to refrence set
    Returns:
        int: number of reads that mapped
        dict: ranges of start and end of each read
        along the genome, with read name as the key
        list: list of dicts with a mapping stats dict per
        read
    """

    counts = defaultdict(int)
    ranges = defaultdict(list)
    sys.stdout.write(query)
    sys.stdout.write(reference)
    a = mp.Aligner(reference, best_n=1, preset=preset)  # load or build index
    # if not a:
    #     raise Exception("ERROR: failed to load/build index")

    mapping_stats = []
    read_count = 0
    for name, seq, qual in mp.fastx_read(query):  # read a fasta/q sequence
        read_count += 1
        stats_dict = {}
        for hit in a.map(seq):  # traverse alignments
            counts[hit.ctg] += 1
            ranges[hit.ctg].append([hit.r_st, hit.r_en])
            stats_dict["name"] = name
            stats_dict["reference"] = hit.ctg
            stats_dict["mapq"] = hit.mapq
            stats_dict["start"] = hit.r_st
            stats_dict["end"] = hit.r_en
            mapping_stats.append(stats_dict)
            break
    return read_count, ranges, mapping_stats


def check_pileup(
    ref: str, ref_ranges: list, reference_file: str, min_coverage=0
) -> list:
    """Iterate over read coordinates and +1 for each position covered to
     build pileup
    Parameters
        ref (str): name of refrence you want to build pileup for
        ref_ranges (list): list of lists specify start and end
        point of each read
        reference_file (str): name of reference file mapped to
    Returns:
        list: coverage at each position of the genome
        position
    """

    if len(ref_ranges) == 0:
        return 0
    for name, seq, qual in mp.fastx_read(reference_file):
        if name == ref:
            coverages = [0] * len(seq)
            for r in ref_ranges:
                for i in range(r[0], r[1]):
                    coverages[i] += 1

    return coverages


def plot_mapping_stats(
    mapping_stats_df: pd.DataFrame, ref_file: str, competitive: bool
) -> px.histogram:
    """Plot histogram of mapping stats
    Parameters
      mapping_stats_df (df): information for each read that mapped
      to any reference provided
      ref_file (str): name of reference file provided
    Returns:
      px.histogram: plotly go scatter plot
    """

    if competitive:
        colours = px.colors.qualitative.Pastel
        total_values = mapping_stats_df["reference"].nunique()
        if total_values > len(colours):
            scale = total_values / len(colours)
            colours = colours * int(math.ceil(scale))

            fig = px.histogram(
                mapping_stats_df,
                x="mapq",
                nbins=6,
                range_x=[0, 60],
                log_y=True,
                color="reference",
                color_discrete_sequence=colours,
                title=f"Read mapping quality to {ref_file}",
            )
    else:
        colours = ["rgb(180, 151, 231)"]
        fig = px.histogram(
            mapping_stats_df,
            x="mapq",
            nbins=6,
            range_x=[0, 60],
            log_y=True,
            color_discrete_sequence=colours,
            title=f"Read mapping quality to {ref_file}",
        )

    return fig


def mapping_stats_text(total: int, ref: str, mapping_stats_df: pd.DataFrame) -> str:
    """Calculates % reads mapped and formats text
    Parameters
        total (int): total number of reads in sample
        mapping_stats_df (df): information for each read that mapped
        to any reference provided
        ref (str): accession of top hit reference
    Returns:
        str: text detailing mapping stats
    """

    count = mapping_stats_df["name"].nunique()
    mapped_pct = count / total * 100
    mapq_count = mapping_stats_df[mapping_stats_df["mapq"] > 50]["name"].nunique()
    text = f"""Read mapping quality to {ref}

{count} of {total} ({round(mapped_pct, 3)}%) reads mapped to a reference in {ref}.
Of the reads that mapped, {mapq_count} mapped with a MAPQ score >50."""

    return text


def plot_pileup(pileup_df: pd.DataFrame, accession: str) -> go.Figure:
    """Plot pileup
    Parameters
        pileup_df (df): coverage at each genome position
        accession (str): accession of top hit reference
    Returns:
        go.Figure: plotly go scatter plot
    """

    fig = go.Figure(
        go.Scatter(
            x=pileup_df["pos"],
            y=pileup_df["coverage"],
            marker={"color": "rgb(102, 197, 204)"},
        )
    )
    fig.update_layout(
        title=f"Coverage of reads mapped to {accession}",
        xaxis_title="Genome Position",
        yaxis_title="Coverage",
        showlegend=False,
    )

    return fig


def pileup_text(pileup_df: pd.DataFrame, ref: str) -> str:
    """Calculates coverage metrics and formats legend
    Parameters
        pileup_df (df): coverage at each genome position
        ref (str): accession of top hit reference
    Returns:
        str: string detailing coverage metrics
        of pileup
    """

    max_x = len(pileup_df.index)

    covered_sites = len(pileup_df[pileup_df["coverage"] > 0])
    pct_ref_cov = (covered_sites / max_x) * 100

    covered_sites_5x = len(pileup_df[pileup_df["coverage"] >= 5])
    pct_ref_cov_5x = (covered_sites_5x / max_x) * 100

    text = f"""Coverage of reads mapped to {ref}

The most reads mapped to {ref}. The sample mapped to {ref} with a coverage of {round(pct_ref_cov, 3)}% across
the genome at 1x and {round(pct_ref_cov_5x, 3)}% at 5x."""

    return text


def commandline():
    parser = argparse.ArgumentParser(
        prog="Mapping",
        epilog="""
            Map reads to reference set using mappy. Plots mapping stats and pileup of top reference hit (optional)
            python mapping.py --sample --fastq --references --pileup
            """,
    )
    parser.add_argument(
        "--sample",
        "-s",
        type=str,
        required=True,
        help="Required: CLIMB-ID of sample",
    )
    parser.add_argument(
        "--fastq",
        "-f",
        type=str,
        required=True,
        help="Required: path to reads",
    )
    parser.add_argument(
        "--references",
        "-r",
        type=str,
        required=True,
        help="Required: path to reference set to map to",
    )
    parser.add_argument(
        "--pileup",
        "-p",
        action="store_true",
        help="Optional: specify if you want to plot pileup of topreference hit",
    )
    parser.add_argument(
        "--competitive",
        "-c",
        action="store_true",
        help="Optional: specify if you are doing competitive mapping",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        required=False,
        default=".",
        help="Optional: directory to save outputs",
    )

    return parser.parse_args()


def main():
    # arguments
    args = commandline()
    sample = args.sample
    refs = args.references
    reads = args.fastq
    pileup = args.pileup
    is_competitive = args.competitive
    save_path = args.output_dir

    sys.stdout.write(mp.__version__)

    # run mapping
    sys.stdout.write(f"mapping {os.path.basename(reads)} to {os.path.basename(refs)}\n")
    total_reads, ranges, mapping_stats = map_to_refs(reads, refs)
    mapping_stats_df = pd.DataFrame.from_dict(mapping_stats, orient="columns")

    # save stats
    ref_file_name = os.path.basename(refs)
    ref_file_name = ref_file_name.split(".")[0]
    stats_prefix = f"{sample}_{ref_file_name}"

    # stats_save_name = os.path.join(save_path, f'{stats_prefix}_mapping_stats.csv')
    stats_save_name = os.path.join(save_path, f"{sample}_mapping_stats.csv")
    sys.stdout.write("saving mapping stats\n")
    mapping_stats_df.to_csv(stats_save_name)

    if not mapping_stats_df.empty:
        # plot mapping stats
        sys.stdout.write("plotting mapping stats\n")
        stats_fig = plot_mapping_stats(mapping_stats_df, ref_file_name, is_competitive)
        stats_text = mapping_stats_text(total_reads, ref_file_name, mapping_stats_df)

        # save outputs
        # stats_fig_save_name = os.path.join(save_path, f'{stats_prefix}_mapping_stats_histogram.html')
        stats_fig_save_name = os.path.join(
            save_path, f"{sample}_mapping_stats_histogram.html"
        )
        stats_fig.write_html(stats_fig_save_name)
    else:
        stats_text = f"No reads mapped to {ref_file_name}"

    # stats_text_save_name = os.path.join(save_path, f'{stats_prefix}_mapping_stats_text.txt')
    stats_text_save_name = os.path.join(save_path, f"{sample}_mapping_stats_text.txt")
    with open(stats_text_save_name, "w") as stats_text_file:
        stats_text_file.write(stats_text)

    # if pileup selected do the plotting
    if pileup:
        sys.stdout.write("plotting pileup\n")
        # find the accession with the most hits
        top_accessions = mapping_stats_df.reference.mode()

        # plot pileup for top accessions
        for accession in top_accessions:
            coverage_list = check_pileup(accession, ranges[accession], refs)
            coverages_df = pd.DataFrame({"coverage": coverage_list})
            coverages_df = coverages_df.reset_index().rename(columns={"index": "pos"})

            # save output
            cov_prefix = f"{sample}_{accession}"
            # cov_save_name = os.path.join(save_path, f'{cov_prefix}_pileup.csv')
            cov_save_name = os.path.join(save_path, f"{sample}_pileup.csv")
            coverages_df.to_csv(cov_save_name)

            # plot pileup
            cov_fig = plot_pileup(coverages_df, accession)
            cov_text = pileup_text(coverages_df, accession)

            # save outputs
            # cov_fig_save_name = os.path.join(save_path, f'{cov_prefix}_pileup.html')
            cov_fig_save_name = os.path.join(save_path, f"{sample}_pileup.html")
            cov_fig.write_html(cov_fig_save_name)

            # cov_text_save_name = os.path.join(save_path, f'{cov_prefix}_pileup_text.txt')
            cov_text_save_name = os.path.join(save_path, f"{sample}_pileup_text.txt")
            with open(cov_text_save_name, "w") as cov_text_file:
                cov_text_file.write(cov_text)


if __name__ == "__main__":
    main()
