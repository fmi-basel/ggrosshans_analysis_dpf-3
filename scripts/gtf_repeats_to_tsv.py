#!/usr/bin/env python

import sys
import os
import pandas as pd
import HTSeq
from argparse import ArgumentParser, RawTextHelpFormatter

def main():
    """ Main function """

    __doc__ = "Filter repeats"

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--gtf",
        dest="gtf",
        help="Input gtf file with repeats",
        required=True,
        metavar="FILE"
    )
     
    parser.add_argument(
        "--tsv",
        dest="tsv",
        help="Output tsv file",
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
        sys.stdout.write(f"Parsing file: {options.gtf} and writing file: {options.tsv} {os.linesep}")

    gtf = HTSeq.GFF_Reader(options.gtf)
    with open(options.tsv, 'w') as w:
        w.write("\t".join(["id", "repName", "repClass", "repFamily" + os.linesep]))
        for gtf_line in gtf:
            w.write("\t".join([gtf_line.attr["gene_id"], 
                            gtf_line.attr["repName"], 
                            gtf_line.attr["repClass"], 
                            gtf_line.attr["repFamily"] + os.linesep]))

    if options.verbose:
        sys.stdout.write(f"Done{os.linesep}")


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
