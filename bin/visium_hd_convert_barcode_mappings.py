#!/usr/bin/env python3
"""10x ships the Visium HD barcode mapping as barcode_mappings.parquet, with one column
per bin resolution (2um -> 8um, 16um, 32um, etc). This script extracts the mapping from
2um to one of those resolutions (e.g. 2um -> 8um) and exports it into a CSV file.
"""
import argparse
import sys

import pandas as pd
import pyarrow.parquet as pq


def main():
    parser = argparse.ArgumentParser(
        description="Subset a Visium HD barcode_mappings.parquet to a single bin resolution, "
        "writing a CSV of the 2um-barcode -> bin assignment with columns 'barcode,bin'"
    )
    parser.add_argument("barcode_mappings", help="Path to barcode_mappings.parquet")
    parser.add_argument("resolution", help="Bin resolution token, e.g. 008um, 016um")
    parser.add_argument("output", help="Output path for the converted CSV")
    args = parser.parse_args()

    # read the schema first so only the two needed columns are loaded from the parquet;
    # the file has one row per 2um barcode and one column per resolution
    schema = pq.read_schema(args.barcode_mappings)
    barcode_column = schema.names[0]
    bin_column = f"square_{args.resolution}"
    if bin_column not in schema.names:
        sys.exit(
            f"error: column '{bin_column}' not found in {args.barcode_mappings} "
            f"(available: {', '.join(schema.names)})"
        )

    barcode_mappings = pd.read_parquet(args.barcode_mappings, columns=[barcode_column, bin_column])
    barcode_mappings = barcode_mappings.rename(columns={barcode_column: "barcode", bin_column: "bin"})
    barcode_mappings.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
