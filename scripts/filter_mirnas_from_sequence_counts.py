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

    __doc__ = "Count occurences of sequences."

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
        "--mirnas",
        dest="mirnas",
        help="Fasta file of annotated miRNAs",
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
        sys.stdout.write(f"Reading miRNAs fasta file: {options.mirnas} {os.linesep}")

    fa = HTSeq.HTSeq.FastaReader(options.mirnas)

    mirnas = {}
    for i in fa:
        mirna = i.seq.decode('UTF8')
        if mirna not in mirnas:
            mirnas[mirna] = mirna

    if options.verbose:
        sys.stdout.write(f"Reading counts table: {options.counts} {os.linesep}")

    df = pd.read_csv(options.counts, header=0, sep="\t")

    if options.verbose:
        sys.stdout.write(f"Filtering mirnas from counts table {os.linesep}")

    df = df[df["Name"].isin(list(mirnas.keys())) == False]

    if options.verbose:
        sys.stdout.write(f"Writing filtered counts table: {options.out} {os.linesep}")

    df.to_csv(options.out, header=True, sep="\t", index=False)

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
