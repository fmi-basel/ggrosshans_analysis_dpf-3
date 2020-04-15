#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------


import sys
import os
import HTSeq
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from argparse import ArgumentParser, RawTextHelpFormatter

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------


def main():
    """ Main function """

    __doc__ = "Filter count sequences."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--counts",
        dest="counts",
        help="Counts table",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--category",
        dest="category",
        help="Category (options: 22G, 21T, 26G)",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--out",
        dest="out",
        help="Output counts file",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        required=False,
        help="Verbose"
    )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # get the arguments
    # -------------------------------------------------------------------------
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if options.verbose:
        sys.stdout.write(f"Reading counts table: {options.counts} {os.linesep}")

    df = pd.read_csv(options.counts, header=0, sep="\t")

    if options.verbose:
        sys.stdout.write(f"Filtering counts table with option: {options.category} {os.linesep}")

    if options.category == "22G":
        df_filter = df[(df["Name"].str.len() == 22) & (df["Name"].str.startswith("G"))].copy()
    elif options.category == "21T":
        df_filter = df[(df["Name"].str.len() == 21) & (df["Name"].str.startswith("T"))].copy()
    elif options.category == "26G":
        df_filter = df[(df["Name"].str.len() == 26) & (df["Name"].str.startswith("G"))].copy()

    if options.verbose:
        sys.stdout.write(f"Writing filtered counts table: {options.out} {os.linesep}")

    df_filter.to_csv(options.out, header=True, sep="\t", index=False)

    if options.verbose:
        sys.stdout.write(f"Done {os.linesep}")

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(0)
