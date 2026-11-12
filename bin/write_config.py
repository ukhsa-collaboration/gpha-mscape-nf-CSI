#!/usr/bin/env python

import argparse
import os
from datetime import datetime

def commandline():
    parser = argparse.ArgumentParser(
        prog="Write config file for report generation",
        epilog="""
            Takes template config file and populates with sample ID. Writes output config file to use for report generation.
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
        "--config",
        "-c",
        type=str,
        required=True,
        help="path to template config file",
            
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

def write_config(config, search_text, replace_text):
    """Edits template config file to generate sample-specific config
    Parameters
        config (str): path to template config file
        search_text (str): placeholder text in the template to be replaced
        replace_text (str): text to insert into the template (sample ID)
    Returns:
        config_text: edited text to be written to the output config file
    """
  
    with open(config, 'r') as file:
        config_text = file.read()
        config_text = config_text.replace(search_text, replace_text)
        now = datetime.now()
        today = now.strftime("%Y-%m-%d")
        config_text = config_text.replace("today", today)
    return config_text


def main():
    args = commandline()
    sample = args.sample
    config = args.config
    save_path = args.output_dir
    config_text = write_config(config, "climb-id", sample)
    save_name = os.path.join(save_path, f"{sample}_config.json")
    with open(save_name, 'w') as file:
        file.write(config_text)

if __name__ == "__main__":
    main()