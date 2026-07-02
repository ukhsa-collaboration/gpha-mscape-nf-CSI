#!/usr/bin/env python

import argparse
import pandas as pd
import os
from taxaplease import TaxaPlease

tp = TaxaPlease()


def parse_kraken_assignments(path: str) -> pd.DataFrame:
    """Read kraken2 std out

    Parameters:
    path (str): path to kraken2 std out in txt/tsv format

    Returns:
    pd.DataFrame: df of kraken2 per read classification

    """
    # read in std out
    kraken_result = pd.read_table(path, header=None)
    kraken_result.columns = ["classification", "read", "taxid", "length", "kmers"]

    sample_id = os.path.basename(path)
    sample_id = sample_id.split("_")[0]
    kraken_result["sample"] = sample_id

    kraken_result["result_dict"] = kraken_result["taxid"].apply(
        lambda x: tp.get_record(x)
    )

    kraken_result["is_dict"] = kraken_result["result_dict"].apply(
        isinstance, args=(dict,)
    )
    kraken_result["result"] = kraken_result["result_dict"].str.get("name")[
        kraken_result["is_dict"] == True
    ]
    kraken_result["result_rank"] = kraken_result["result_dict"].str.get("rank")[
        kraken_result["is_dict"] == True
    ]

    kraken_result.loc[kraken_result["classification"] == "U", "result"] = "Unclassified"
    kraken_result["result"] = kraken_result["result"].fillna("")

    kraken_result.drop(columns=["result_dict", "is_dict"], axis=1, inplace=True)

    return kraken_result


def group_taxids(df: pd.DataFrame, rank: str) -> pd.DataFrame:
    """Get name for specified rank for each result

    Parameters:
    df (pd.DataFrame): kraken2 classification per read
    rank (str): taxonomic rank to collapse taxid to

    Returns:
    pd.DataFrame: df of kraken2 per read classification with
    rank column added
    """
    # use taxa please to get higher rank taxid and name
    df["result_dict"] = df["taxid"].apply(
        lambda x: tp.get_specified_rank_record(x, rank)
    )
    # df['result_grouped'] = df['result_grouped'].apply(lambda x:tp.get_record(x))

    df["is_dict"] = df["result_dict"].apply(isinstance, args=(dict,))
    df["taxid_grouped"] = df["result_dict"].str.get("taxid")[df["is_dict"] == True]
    df["taxid_grouped"] = df["taxid_grouped"].fillna(df["taxid"])

    df["result_grouped"] = df["result_dict"].str.get("name")[df["is_dict"] == True]
    df["result_grouped"] = df["result_grouped"].fillna(df["result"])

    # if phage just call as phage
    df["result_grouped"] = df.apply(
        lambda x: "Phage" if tp.isPhage(x["taxid"]) == True else x["result_grouped"],
        axis=1,
    )

    df["group_rank"] = df["result_dict"].str.get("rank")[df["is_dict"] == True]
    df["group_rank"] = df["group_rank"].fillna(df["result_rank"])
    df["group_rank"] = df["group_rank"].fillna("")

    df.drop(columns=["is_dict", "result_dict"], axis=1, inplace=True)

    return df


def add_category(df: pd.DataFrame) -> pd.DataFrame:
    # get unclassified
    df.loc[df["taxid_grouped"] == 0, "category"] = "Unclassified"

    # sort classifications into groups
    df["category"] = df.apply(
        lambda x: "Phage" if tp.isPhage(x["taxid_grouped"]) == True else x["category"],
        axis=1,
    )

    df["category"] = df.apply(
        lambda x: "Virus" if tp.isVirus(x["taxid_grouped"]) == True else x["category"],
        axis=1,
    )

    df["category"] = df.apply(
        lambda x: "Bacteria"
        if tp.isBacteria(x["taxid_grouped"]) == True
        else x["category"],
        axis=1,
    )

    df["category"] = df.apply(
        lambda x: "Archaea"
        if tp.isArchaea(x["taxid_grouped"]) == True
        else x["category"],
        axis=1,
    )

    df["category"] = df.apply(
        lambda x: "Eukaryote"
        if tp.isEukaryote(x["taxid_grouped"]) == True
        else x["category"],
        axis=1,
    )

    df.loc[df["category"].isna(), "category"] = "Other"

    return df


def summarise_results(df: pd.DataFrame) -> pd.DataFrame:
    """Summarise results and get seq length stats

    Parameters:
    df (pd.DataFrame): kraken2 classification per read

    Returns:
    pd.DataFrame: df of kraken2 results grouped by
    classification with seq length stats per group.
    """
    df_grouped = (
        df.groupby(["taxid", "result", "result_grouped"])
        .size()
        .reset_index(name="count")
    )
    df_grouped_stats = (
        df.groupby(["result", "result_grouped"])[["length"]]
        .agg(["mean", "std"])
        .reset_index()
    )

    df_grouped_stats.columns = list(map("_".join, df_grouped_stats.columns.values))
    df_grouped_stats = df_grouped_stats.rename(
        columns={
            "result_": "result",
            "result_grouped_": "result_grouped",
            "taxid_": "result_grouped",
        }
    )

    df_combined = df_grouped.merge(
        df_grouped_stats, on=["result", "result_grouped"], how="left"
    )

    df_combined["length_mean"] = round(df_combined["length_mean"])
    df_combined["length_std"] = round(df_combined["length_std"])

    return df_combined


def kraken_text(database: str, rank: str, input_type: str) -> str:
    """Format title and text to go with Kraken2 summaries
    Parameters
        database (str): name of database used for classification
        rank (str): taxonomic rank classifications have been grouped
        to
        input_type (str): read or contigs classified
    Returns:
        text (str): string summarising Kraken2 results
    """

    text = f"""Kraken2 with {database} {input_type} classifications

{input_type.capitalize()} classified using Kraken2 with {database}. Classifications
are grouped to {rank} in result_grouped column. Where a classification
is above the taxonomic rank {rank} the result has not been grouped. Where
there are >10 {input_type}s classifications are grouped with summary data
provided."""

    return text


def commandline():
    parser = argparse.ArgumentParser(
        prog="Kraken2 Summary",
        epilog="""
            Summarise Kraken2 stdout for a given sample
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
        "--input_type",
        "-i",
        type=str,
        required=True,
        help="Required: what type of sample - read or contigs?",
    )
    parser.add_argument(
        "--kraken_results",
        "-k",
        type=str,
        required=True,
        help="Required: path to Kraken2 stdout.",
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
        "--database",
        "-d",
        type=str,
        required=True,
        help="Required: name of the database used when running Kraken2.",
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
    input_type = args.input_type
    kraken_results = args.kraken_results
    rank = args.rank
    database = args.database
    save_path = args.output_dir

    kraken_formatted = parse_kraken_assignments(kraken_results)
    kraken_formatted = group_taxids(kraken_formatted, rank)
    kraken_formatted = add_category(kraken_formatted)

    group_save_name = os.path.join(
        save_path, f"{sample}_{input_type}_kraken2_{database}_{rank}.csv"
    )
    kraken_formatted.to_csv(group_save_name, index=False)

    kraken_formatted_text = kraken_text(database, rank, input_type)

    kraken_text_save_name = os.path.join(
        save_path, f"{sample}_{input_type}_kraken2_{database}_text.txt"
    )
    with open(kraken_text_save_name, "w") as kraken_text_file:
        kraken_text_file.write(kraken_formatted_text)

    total = kraken_formatted["read"].nunique()
    if total > 10:
        kraken_table = summarise_results(kraken_formatted)
    else:
        kraken_table = kraken_formatted[
            ["read", "taxid", "result", "result_grouped", "length"]
        ]

    table_save_name = os.path.join(
        save_path, f"{sample}_{input_type}_kraken2_{database}_summary.csv"
    )
    kraken_table.to_csv(table_save_name, index=False)


if __name__ == "__main__":
    main()
