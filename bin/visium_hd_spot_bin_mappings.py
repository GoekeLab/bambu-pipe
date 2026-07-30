#!/usr/bin/env python3
import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Map every 2um spot bambu quantified to its bin at one resolution")
    parser.add_argument("barcode_mappings", help="Path to the resolution's barcode mappings CSV (columns: barcode,bin)")
    parser.add_argument("barcodes", help="Path to barcodes.tsv.gz from the 2um count directory; one bambu id per column")
    parser.add_argument("sample", help="Sample name that bambu prefixes onto each barcode to form the id")
    parser.add_argument("output", help="Output path for the spot to bin mapping CSV")
    args = parser.parse_args()

    # bambu ids are sampleName_barcode, so strip the prefix to get back to the mapping's barcodes
    prefix = f"{args.sample}_"
    spot_barcodes = pd.read_csv(args.barcodes, header=None)[0].str.removeprefix(prefix)

    # keep the bins of the spots bambu quantified, then add the prefixed forms so
    # downstream modules never have to rebuild a bambu id
    barcode_mappings = pd.read_csv(args.barcode_mappings)
    spot_mappings = barcode_mappings[barcode_mappings["barcode"].isin(spot_barcodes)].copy()
    spot_mappings["sample_barcode"] = prefix + spot_mappings["barcode"]
    spot_mappings["sample_bin"] = prefix + spot_mappings["bin"]

    spot_mappings[["sample_barcode", "barcode", "sample_bin", "bin"]].to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
