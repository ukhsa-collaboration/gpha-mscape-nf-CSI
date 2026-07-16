#!/usr/bin/env python

import argparse
import pandas as pd
import os
import plotly.graph_objects as go
from plotly.io import templates

templates.default = "plotly_white"


def parse_all_kraken_results(outputs):
    """Parse kraken summary outputs
    from kraken_summary.py

    Parameters:
    outputs (list): paths to kraken summaries

    Returns:
    pd.DataFrame: df of kraken2 classifications
    """

    all_results_list = []
    for path in outputs:
        df = pd.read_csv(path)
        filename = os.path.basename(path)
        df["database"] = filename.split("_")[3]
        df["input"] = filename.split("_")[1]

        all_results_list.append(df)

    all_results = pd.concat(all_results_list)

    return all_results


def process_read_class(
    tobacco_check: pd.DataFrame, all_read_results_df: pd.DataFrame
) -> pd.DataFrame:
    """Formats read classifications from different databases,
    combines into single df with 1 row per path in Sankey
    and a column per x point. Then counts each path.

    Parameters:
    tobacco_check (pd.DataFrame): mapping to tobacco genome
    results
    all_read_results_df (pd.DataFrame): per read classification
    results for all databases used

    Returns:
    pd.DataFrame: df with all database results formatted, combined
    and counted (each row will be a path in the Sankey
    """
    viraldb = all_read_results_df[all_read_results_df["database"] == "viper"]
    pluspf23 = all_read_results_df[all_read_results_df["database"] == "pluspf2023"]
    pluspf25 = all_read_results_df[all_read_results_df["database"] == "pluspf2025"]

    # each result needs unique field name for each x position
    viraldb_renamed = viraldb.rename({"result_grouped": "viper_result"}, axis=1)

    pluspf23_renamed = pluspf23.rename({"result_grouped": "pluspf23_result"}, axis=1)

    pluspf25_renamed = pluspf25.rename({"result_grouped": "pluspf25_result"}, axis=1)

    # tobacco needs a 'result'
    tobacco_check = tobacco_check.assign(tobacco_mapping="Mapped to tobacco genome")

    # put them all to together in single df
    results_merge = pluspf23_renamed[["read", "pluspf23_result"]].merge(
        viraldb_renamed[["read", "viper_result"]], on="read", how="left"
    )

    results_merge = results_merge.merge(
        pluspf25_renamed[["read", "pluspf25_result"]], on="read", how="left"
    )
    results_merge = results_merge.merge(
        tobacco_check[["name", "tobacco_mapping"]],
        left_on="read",
        right_on="name",
        how="left",
    )

    # bit more formatting
    results_merge.loc[
        results_merge["pluspf23_result"] == "Unclassified", "tobacco_mapping"
    ] = "Mapped to tobacco genome"
    results_merge.loc[
        results_merge["pluspf23_result"] != "Unclassified", "tobacco_mapping"
    ] = "Not mapped to tobacco genome"

    # only want to keep unclassified pluspf 23 results
    results_merge = results_merge[results_merge["pluspf23_result"] == "Unclassified"]

    results_merge_count = (
        results_merge.groupby(
            ["pluspf23_result", "tobacco_mapping", "viper_result", "pluspf25_result"]
        )
        .size()
        .reset_index(name="count")
    )

    return results_merge_count


def process_contig_class(
    pluspf23_read: pd.DataFrame,
    contig_mapping: pd.DataFrame,
    all_contig_results_df: pd.DataFrame,
) -> pd.DataFrame:
    """Formats contig classifications from different databases,
    combines into single df with 1 row per path in Sankey
    and a column per x point. Then counts each path.

    Parameters:
    pluspf23_read (pd.DataFrame): df of pluspf 2023 results per read
    contig_mapping (pd.DataFrame): df of results mapping reads -> contig
    all_contig_results_df (pd.DataFrame): per contig classification
    results for all databases used

    Returns:
    pd.DataFrame: df with all database results formatted, combined
    and counted (each row will be a path in the Sankey)
    """

    viraldb = all_contig_results_df[all_contig_results_df["database"] == "viper"]
    pluspf23 = all_contig_results_df[all_contig_results_df["database"] == "pluspf2023"]
    pluspf25 = all_contig_results_df[all_contig_results_df["database"] == "pluspf2025"]

    viraldb_renamed = viraldb.rename(
        {"result_grouped": "viper_result", "read": "contig"}, axis=1
    )

    pluspf23_read_renamed = pluspf23_read.rename(
        {"result_grouped": "pluspf23_read_result"}, axis=1
    )

    pluspf23_renamed = pluspf23.rename(
        {"result_grouped": "pluspf23_result", "read": "contig"}, axis=1
    )

    pluspf25_renamed = pluspf25.rename(
        {"result_grouped": "pluspf25_result", "read": "contig"}, axis=1
    )

    contig_mapping_renamed = contig_mapping.rename({"reference": "contig"}, axis=1)

    # put them all to together in single df
    results_merge = pluspf23_read_renamed[["read", "pluspf23_read_result"]].merge(
        contig_mapping_renamed[["name", "contig"]],
        left_on="read",
        right_on="name",
        how="left",
    )

    test = results_merge.merge(
        viraldb_renamed[["contig", "viper_result"]], on="contig", how="left"
    )

    results_merge = results_merge.merge(
        viraldb_renamed[["contig", "viper_result"]], on="contig", how="left"
    )

    results_merge = results_merge.merge(
        pluspf23_renamed[["contig", "pluspf23_result"]], on="contig", how="left"
    )

    results_merge = results_merge.merge(
        pluspf25_renamed[["contig", "pluspf25_result"]], on="contig", how="left"
    )

    # only want to keep unclassified pluspf 23 results
    results_merge = results_merge[
        results_merge["pluspf23_read_result"] == "Unclassified"
    ]  # would have to change this for the HCID one

    results_merge_count = (
        results_merge.groupby(
            [
                "pluspf23_read_result",
                "contig",
                "pluspf23_result",
                "viper_result",
                "pluspf25_result",
            ],
            dropna=False,
        )
        .size()
        .reset_index(name="count")
    )
    # makes sure we keep reads that did not assemble into contigs
    results_merge_count = results_merge_count.fillna({"contig": "Not assembled"})
    results_merge_count = results_merge_count.fillna("empty node")

    return results_merge_count


def filter_col(series: pd.Series) -> pd.Series:
    """remove any nodes that are empty at
    each x position"""
    return series[series != "empty node"]


def generate_sankey_chart_data(
    df: pd.DataFrame, columns: list, sankey_link_weight: str
):
    """Find indices for each node and creates links between source and targets.
    Removes any empty nodes so they are not included in the plot.

    Parameters:
    df (pd.DataFrame): from process_x_class functions. 1 row per path in
    sankey with each column corresponding to an x position. Count column
    used as weight for sankey.
    columns (list): list of df column names (each x position)
    sankey_link_weight (str): name of column specifying link weight.

    Returns:
    A list for node labels, source and target indices, and link weight.
    Each list to be used in plotly go sankey function.
    """

    # list of list: each list are the set of nodes in each tier/column
    column_values_all = [df[col] for col in columns]

    column_values = [filter_col(col) for col in column_values_all]

    # this generates the labels for the sankey by taking all the unique values
    labels = sum([list(node_values.unique()) for node_values in column_values], [])

    # initializes a dict of dicts (one dict per tier)
    link_mappings = {col: {} for col in columns}

    # each dict maps a node to unique number value (same node in different tiers
    # will have different nubmer values
    i = 0
    for col, nodes in zip(columns, column_values):
        for node in nodes.unique():
            link_mappings[col][node] = i
            i = i + 1

    # specifying which coluns are serving as sources and which as sources
    # ie: given 3 df columns (col1 is a source to col2, col2 is target to col1 and
    # a source to col 3 and col3 is a target to col2
    source_nodes = column_values[: len(columns) - 1]
    target_nodes = column_values[1:]
    source_cols = columns[: len(columns) - 1]
    target_cols = columns[1:]
    links = []

    # loop to create a list of links in the format [((src,tgt),wt),(),()...]
    for source, target, source_col, target_col in zip(
        source_nodes, target_nodes, source_cols, target_cols
    ):
        for val1, val2, link_weight in zip(source, target, df[sankey_link_weight]):
            links.append(
                (
                    (link_mappings[source_col][val1], link_mappings[target_col][val2]),
                    link_weight,
                )
            )

    # creating a dataframe with 2 columns: for the links (src, tgt) and weights
    df_links = pd.DataFrame(links, columns=["link", "weight"])

    # aggregating the same links into a single link (by weight)
    df_links = df_links.groupby(by=["link"], as_index=False).agg({"weight": "sum"})

    # generating three lists needed for the sankey visual
    sources = [val[0] for val in df_links["link"]]
    targets = [val[1] for val in df_links["link"]]
    weights = df_links["weight"]

    return labels, sources, targets, weights


def make_sankey(
    df: pd.DataFrame, title: str, taxids_dict: dict, rank_dict: dict, colour_dict: dict
):
    """Calls function to find labels, sources, targets, weights and
    uses to plot sankey. Dicts are filtered to the availables labels
    and used to plot colours/hover data.

    Parameters:
    df (pd.DataFrame): from process_x_class functions. 1 row per path in
    sankey with each column corresponding to an x position. Count column
    used as weight for sankey. Fed to generate_sankey_chart_data method.
    title (str): name of plot (read or contig?)
    taxids_dict (dict): dict of node labels and corresponding taxid (None
    where no taxid).
    rank_dict (dict): dict of node labels and corresponding taxonomic rank
    (None where no rank).
    rank_dict (dict): dict of node labels and corresponding colour based on
    a category.

    Returns:
    Sankey diagram figure
    """

    cols = list(df.columns)[:-1]
    labels, sources, targets, weights = generate_sankey_chart_data(df, cols, "count")

    # filter dicts using labels
    taxids = [taxids_dict[name] for name in labels]
    ranks = [rank_dict[name] for name in labels]
    custom = [[taxids[val], ranks[val]] for val in range(len(labels))]
    node_colours = [colour_dict[name] for name in labels]

    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    # line = dict(width = 0.5),
                    label=labels,
                    customdata=custom,
                    hovertemplate="Result: %{label}<br />"
                    + "Total Reads: %{value}<br />"
                    + "Taxid: %{customdata[0]}<br />"
                    + "Rank: %{customdata[1]}<extra></extra>",
                    color=node_colours,
                ),
                link=dict(
                    source=sources,  # indices correspond to labels, eg A1, A2, A1, B1, ...
                    target=targets,
                    value=weights,
                ),
            )
        ]
    )

    for x_coordinate, column_name in enumerate(cols):
        fig.add_annotation(
            x=x_coordinate,
            y=1.05,
            xref="x",
            yref="paper",
            text=column_name.replace("_", " "),
            showarrow=False,
            font=dict(
                size=13,
            ),
            align="center",
        )

    fig.update_layout(
        xaxis={
            "showgrid": False,
            "zeroline": False,
            "visible": False,
        },
        yaxis={
            "showgrid": False,
            "zeroline": False,
            "visible": False,
        },
        plot_bgcolor="rgba(0,0,0,0)",
        font_size=10,
        hovermode="x",
        title={"text": title, "font": dict(size=16)},
    )

    fig.update_layout(autosize=False, height=800, width=1100)

    # display(HTML(fig.to_html()))

    return fig


def colours(df: pd.DataFrame) -> dict:
    """Assigns a category and colour based on
    classification/node labels. Mapping results
    manually added.

    Parameters:
    df (pd.DataFrame): all classification results concatenated

    Returns:
    dict: dictionary with node labels and assigned colours
    """
    # colours for both read and contig plots
    colour_cat_dict = {
        "Unclassified": "rgb(158, 185, 243)",
        "contig": "rgb(102, 197, 204)",
        "Not assembled": "rgb(139, 224, 164)",
        "Mapped to tobacco genome": "rgb(139, 224, 164)",
        "Not mapped to tobacco genome": "rgb(102, 197, 204)",
        "Virus": "rgb(246, 207, 113)",
        "Phage": "rgb(248, 156, 116)",
        "Bacteria": "rgb(180, 151, 231)",
        "Archaea": "rgb(201, 219, 116)",
        "Eukaryote": "rgb(254, 136, 177)",
        "Other": "rgb(201, 219, 116)",
    }

    df["colour"] = df["category"].map(colour_cat_dict)

    # make dictionary for classifications
    colours_mapped = df[["result_grouped", "colour"]].drop_duplicates()
    colours_dict = dict(zip(colours_mapped["result_grouped"], colours_mapped["colour"]))

    # need colours for other categories eg., tobacco mapping, not assembled etc.
    # check if it's contig data (flye calls each contig 'contig_n')
    contigs = df[df["read"].str.contains("contig")][["read"]].drop_duplicates()
    if not contigs.empty:
        contig_colours_dict = dict(
            zip(contigs["read"], [colour_cat_dict["contig"]] * len(contigs))
        )
        colours_dict = colours_dict | contig_colours_dict
        colours_dict["Not assembled"] = colour_cat_dict["Not assembled"]
    else:  # if no reads labelled with 'contig' must be reads only
        colours_dict["Mapped to tobacco genome"] = colour_cat_dict[
            "Mapped to tobacco genome"
        ]
        colours_dict["Not mapped to tobacco genome"] = colour_cat_dict[
            "Not mapped to tobacco genome"
        ]

    return colours_dict


def sankey_text(input_type: str, rank: str) -> str:
    """Write text to accomany sankey diagrams

    Parameters:
    input_type (str): read or contigs
    rank (str): taxonomic rank results have been collapsed
    upto in previous parts of pipeline

    str (text): legend to go with sankey
    """

    text = f"""{input_type} Kraken2 Classification Results

Kraken2 can classify {input_type}s using different databases. Here, {input_type}s are classified using
using PlusPF 2023, PlusPF 2025 and a curated viral database (viper). Taxonomic classifications
assigned by Kraken2 have been collapsed up to {rank} to simplify the diagram. Where a {input_type} has a higher
taxonomic classification than {rank} level, the initial taxonomic classification is used.
{input_type.capitalize()}s classified as any bacteriophage are labelled as 'Phage'. The thickness of
each path represents the number of reads whilst the colour of each node represents a category: unclassified,
contigs, unassembled reads, reads that mapped to the tobacco genome, reads that did not map to the tobacco
genome, virus, phage, bacteria, archaea, eukaryote and other."""

    return text


def commandline():
    parser = argparse.ArgumentParser(
        prog="Kraken2 Classification Sankey",
        epilog="""
            Plot classifications
            """,
    )
    parser.add_argument(
        "--sample",
        "-s",
        type=str,
        required=True,
        help="Required: CLIMB-ID of sample.",
    )
    parser.add_argument(
        "--rank",
        "-r",
        type=str,
        required=True,
        help="Required: taxonomic rank used to summarise"
        "Kraken2 classification results.",
    )
    parser.add_argument(
        "--tobacco_mapping",
        "-t",
        type=str,
        required=True,
        help="Required: path tobacco mapping outputs",
    )
    parser.add_argument(
        "--contigs",
        "-c",
        type=str,
        required=True,
        help="Required: path to assembly summary output",
    )
    parser.add_argument(
        "--kraken2_read_list",
        metavar="N",
        type=str,
        nargs="+",
        required=True,
        help="Required: path to kraken2 read summary outputs",
    )
    parser.add_argument(
        "--kraken2_contig_list",
        metavar="N",
        type=str,
        nargs="+",
        required=True,
        help="Required: path to kraken2 contig summary outputs",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        required=False,
        default=".",
        help="Optional: directory to save outputs.",
    )

    return parser.parse_args()


def main():
    args = commandline()
    sample = args.sample
    rank = args.rank
    tobacco_mapping = args.tobacco_mapping
    contig_mapping = args.contigs
    kraken2_read_summary = args.kraken2_read_list
    kraken2_contig_summary = args.kraken2_contig_list
    save_path = args.output_dir

    ### read level sankey ###
    tobacco_mapped_reads = pd.read_csv(tobacco_mapping)
    try:
        tobacco_mapped_reads = tobacco_mapped_reads[tobacco_mapped_reads["mapq"] > 50]
    except:
        tobacco_mapped_reads = pd.DataFrame(columns=["name", "reference", "mapq"])

    all_read_results = parse_all_kraken_results(kraken2_read_summary)

    # return all_read_results

    # dictionary to add hover data (taxid and rank)
    read_taxid_dict = dict(
        zip(all_read_results["result_grouped"], all_read_results["taxid_grouped"])
    )
    read_taxid_dict["Mapped to tobacco genome"] = None
    read_taxid_dict["Not mapped to tobacco genome"] = None

    read_rank_dict = dict(
        zip(all_read_results["result_grouped"], all_read_results["group_rank"])
    )
    read_rank_dict["Mapped to tobacco genome"] = None
    read_rank_dict["Not mapped to tobacco genome"] = None

    # dictionary to add colours
    read_colours = colours(all_read_results)

    # format data
    read_results = process_read_class(tobacco_mapped_reads, all_read_results)
    # do the plotting
    read_sankey = make_sankey(
        read_results,
        "Read Classifications",
        read_taxid_dict,
        read_rank_dict,
        read_colours,
    )

    read_sankey_text = sankey_text("read", rank)

    # save the outputs
    read_sankey_save_name = os.path.join(save_path, f"{sample}_read_sankey.html")
    read_sankey.write_html(read_sankey_save_name)

    read_sankey_save_name = os.path.join(save_path, f"{sample}_read_sankey.json")
    read_sankey.write_json(read_sankey_save_name)

    read_text_save_name = os.path.join(save_path, f"{sample}_read_sankey_text.txt")
    with open(read_text_save_name, "w") as read_text_file:
        read_text_file.write(read_sankey_text)

    ### contig level sankey ###
    # do it all again just a lil bit different
    contig_mapping = pd.read_csv(contig_mapping)

    all_contig_results = parse_all_kraken_results(kraken2_contig_summary)

    # dictionaries for hover data
    contig_taxid_dict = dict(
        zip(all_contig_results["result_grouped"], all_contig_results["taxid_grouped"])
    )
    taxid_contig_only = dict(
        zip(contig_mapping["reference"], [None] * len(contig_mapping))
    )
    contig_taxid_dict = contig_taxid_dict | taxid_contig_only
    contig_taxid_dict["Not assembled"] = None

    contig_rank_dict = dict(
        zip(all_contig_results["result_grouped"], all_contig_results["group_rank"])
    )
    rank_contig_only = dict(
        zip(contig_mapping["reference"], [None] * len(contig_mapping))
    )
    contig_rank_dict = contig_rank_dict | rank_contig_only
    contig_rank_dict["Not assembled"] = None

    # colour per category
    contig_colours = colours(all_contig_results)

    # need read results to show which went to what contigs
    pluspf23_read_results = all_read_results[
        all_read_results["database"] == "pluspf2023"
    ]

    contig_results = process_contig_class(
        pluspf23_read_results, contig_mapping, all_contig_results
    )

    contig_sankey = make_sankey(
        contig_results,
        "Contig Classifications",
        contig_taxid_dict,
        contig_rank_dict,
        contig_colours,
    )

    contig_sankey_text = sankey_text("contig", rank)

    # save the outputs
    contig_sankey_save_name = os.path.join(save_path, f"{sample}_contig_sankey.html")
    contig_sankey.write_html(contig_sankey_save_name)

    contig_sankey_save_name = os.path.join(save_path, f"{sample}_contig_sankey.json")
    contig_sankey.write_json(contig_sankey_save_name)

    contig_text_save_name = os.path.join(save_path, f"{sample}_contig_sankey_text.txt")
    with open(contig_text_save_name, "w") as contig_text_file:
        contig_text_file.write(contig_sankey_text)


if __name__ == "__main__":
    main()
