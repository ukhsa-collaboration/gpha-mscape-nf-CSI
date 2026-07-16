#!/usr/bin/env python

import argparse
import pandas as pd
import os
import plotly.express as px
from plotly.io import templates
import mapping


templates.default = "plotly_white"


def map_reads_to_contigs(reads: str, contigs: str) -> tuple["int", pd.DataFrame]:
    """Map reads to contigs to calculate how
    many reads assembled
    Parameters
        reads (str): path to reads fastq
        contigs (str): path to contigs
    Returns:
        int: total reads in sample
        pd.DataFrame: information for each read that mapped
        to any reference provided
    """

    total_reads, ranges, mapping_stats = mapping.map_to_refs(reads, contigs)
    mapping_stats_df = pd.DataFrame.from_dict(mapping_stats, orient="columns")

    return total_reads, mapping_stats_df


def parse_flye_stats(path: str) -> pd.DataFrame:
    """Parse assembly_info.txt from Flye
    Parameters
        path (str): path to assembly_info.txt
    Returns:
        pd.DataFrame: length of each contig
    """

    contig_names = []
    contig_length = []
    with open(path, "r") as file:
        flye_results = file.readlines()
        for i in range(1, len(flye_results)):
            contig_names.append(flye_results[i].strip().split("\t")[0])
            contig_length.append(int(flye_results[i].strip().split("\t")[1]))

    contig_info = pd.DataFrame({"contig": contig_names, "length": contig_length})

    return contig_info


def plot_assembly_stats(contig_df: pd.DataFrame) -> px.histogram:
    """Plot contig lengths
    Parameters
        contig_df (df): length of each contig
    Returns:
        px.histogram: histogram of contig lengths
    """

    colours = ["rgb(158, 185, 243)"]
    fig = px.histogram(contig_df, x="length", color_discrete_sequence=colours)
    fig.update_layout(title="Contig Length")

    return fig


def assembly_stats_text(
    contig_df: pd.DataFrame, total_reads: int, mapping_stats_df: pd.DataFrame
) -> str:
    """Calculate numner of reads assembled to contigs and other stats,
    then format text.
    Parameters
        contig_df (df): length of each contig
        total reads (int): total reads in sample
        mapping_stats_df (df): information for each read that mapped
        to any reference provided
    Returns:
        str: text detailing fly assembly results
    """

    count = mapping_stats_df["name"].nunique()
    assembled_pct = count / total_reads * 100
    total_contigs = contig_df["contig"].nunique()
    longest = contig_df["length"].max()
    shortest = contig_df["length"].min()
    mean = contig_df["length"].mean()

    text = f"""Contig Assembly Statistics

{count} reads assembled into {total_contigs} contigs ({round(assembled_pct)}%) with an average contig length of {round(mean)} bases.
The longest contig was {longest} bases in length, whilst the shortest contig was {shortest} bases in length."""

    return text


def commandline():
    parser = argparse.ArgumentParser(
        prog="Flye de novo assembly summary",
        epilog="""
            Calculate how many reads assembled into contigs and plot
            assembly information.
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
        "--reads",
        "-r",
        type=str,
        required=True,
        help="Required: path to reads",
    )
    parser.add_argument(
        "--assembly_fasta",
        "-f",
        type=str,
        required=True,
        help="Required: path to Flye output assembly.fasta",
    )
    parser.add_argument(
        "--assembly_info",
        "-i",
        type=str,
        required=True,
        help="Required: path to Flye output assembly_info.txt.Includes contig lengths.",
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
    args = commandline()
    sample = args.sample
    assembly_info = args.assembly_info
    assembly_fasta = args.assembly_fasta
    reads = args.reads
    save_path = args.output_dir

    # map reads to contigs to see if they all assembled
    total_reads, mapping_stats_df = map_reads_to_contigs(reads, assembly_fasta)

    # save outputs
    table_save_name = os.path.join(save_path, f"{sample}_read2contig_mapping_stats.csv")
    mapping_stats_df.to_csv(table_save_name)

    # look at contig info given by flye
    contig_info = parse_flye_stats(assembly_info)
    stats_fig = plot_assembly_stats(contig_info)

    # combine information together into text
    stats_text = assembly_stats_text(contig_info, total_reads, mapping_stats_df)

    # save outputs
    stats_fig_save_name = os.path.join(
        save_path, f"{sample}_de_novo_assembly_stats.html"
    )
    stats_fig.write_html(stats_fig_save_name)

    stats_text_save_name = os.path.join(
        save_path, f"{sample}_de_novo_assembly_stats_text.txt"
    )
    with open(stats_text_save_name, "w") as stats_text_file:
        stats_text_file.write(stats_text)


if __name__ == "__main__":
    main()
