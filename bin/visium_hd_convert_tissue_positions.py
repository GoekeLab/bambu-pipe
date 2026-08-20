#!/usr/bin/env python3
"""10x ships the Visium HD spatial metadata as tissue_positions.parquet, one file per bin
resolution. This script converts the parquet file to CSV so that every downstream step can
read it. The 2um tissue_positions file is also used to generate a list of in-tissue
barcodes, which filters the BAM so that it only carries barcodes that are in the tissue.
"""
import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Convert a Visium HD tissue_positions.parquet to CSV; for the 2um resolution also write the in-tissue barcode list used to filter the BAM")
    parser.add_argument("tissue_positions", help="Path to tissue_positions.parquet")
    parser.add_argument("resolution", help="Bin resolution token, e.g. 002um, 008um, 016um")
    parser.add_argument("csv", help="Output path for the converted CSV")
    parser.add_argument("barcodes", help="Output path for the in-tissue barcode list (written only for the 002um resolution)")
    args = parser.parse_args()

    tissue_positions = pd.read_parquet(args.tissue_positions)
    tissue_positions.to_csv(args.csv, index=False)

    # the BAM carries 2um barcodes, so only the 2um resolution feeds the samtools tissue filter
    if args.resolution == "002um":
        in_tissue = tissue_positions.loc[tissue_positions["in_tissue"] == 1, "barcode"]
        in_tissue.to_csv(args.barcodes, index=False, header=False)


if __name__ == "__main__":
    main()
