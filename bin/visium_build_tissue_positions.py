#!/usr/bin/env python3
"""Build the tissue positions CSV (Space Ranger format) for visium-v* samples from the
Loupe manual alignment JSON and the chemistry's coordinates file. Also writes the
in-tissue barcode list used by samtools to filter out-of-tissue barcodes from the BAM.
"""
import argparse
import csv
import json
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("loupe_alignment", help="Loupe manual alignment .json")
    parser.add_argument("spatial_coordinates", help="chemistry coordinates file (barcode, x, y; 1-based)")
    parser.add_argument("tissue_positions", help="output tissue positions CSV")
    parser.add_argument("tissue_barcodes", help="output in-tissue barcode list")
    args = parser.parse_args()

    with open(args.loupe_alignment) as f:
        oligos = json.load(f)["oligo"]
    spots = {(o["row"], o["col"]): o for o in oligos}

    coords = {}
    with open(args.spatial_coordinates) as f:
        for line in f:
            barcode, x, y = line.split()
            coords[(int(y) - 1, int(x) - 1)] = barcode

    # the JSON grid and the (shifted) coordinates grid must coincide exactly
    if set(spots) != set(coords):
        sys.exit(
            f"error: the {len(spots)} spots in {args.loupe_alignment} do not match the "
            f"{len(coords)}-spot layout in {args.spatial_coordinates} -- check that the "
            "alignment JSON was made for this slide chemistry"
        )

    n_tissue = 0
    with open(args.tissue_positions, "w", newline="") as fc, open(args.tissue_barcodes, "w") as fb:
        writer = csv.writer(fc)
        writer.writerow(["barcode", "in_tissue", "array_row", "array_col",
                         "pxl_row_in_fullres", "pxl_col_in_fullres"])
        for (row, col), oligo in sorted(spots.items(), key=lambda item: item[0]):
            in_tissue = int(bool(oligo.get("tissue")))
            writer.writerow([coords[(row, col)], in_tissue, row, col,
                             oligo["imageY"], oligo["imageX"]])
            if in_tissue:
                fb.write(coords[(row, col)] + "\n")
                fb.write(coords[(row, col)] + "-1\n")
                n_tissue += 1

    if n_tissue == 0:
        sys.exit("error: no spot in the alignment JSON is marked as tissue")
    print(f"{n_tissue} of {len(spots)} spots are in tissue")


if __name__ == "__main__":
    main()
