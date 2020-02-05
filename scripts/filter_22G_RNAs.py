#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------


import sys
import os
import HTSeq
import pandas as pd
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
        "--bam",
        dest="bam",
        help="Alignment file in BAM (sorted and indexed) format. Output of bowtie.",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--bam_out",
        dest="bam_out",
        help="Output filtered bam file",
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

    bam = HTSeq.BAM_Reader(options.bam)
    bam_writer = HTSeq.BAM_Writer.from_BAM_Reader(options.bam_out, bam)

    # parse bam file based on chromosomes
    if options.verbose:
        sys.stdout.write(f"Fetching reads from bam file: {options.bam}, filtering 22G RNAs and writing to bam file: {options.bam_out} {os.linesep}")

    for aln in bam.fetch():
        seq = aln.read.seq.decode('utf-8')
        
        if (len(seq) == 22) and (seq[0]=="G") and (aln.aligned):
            bam_writer.write(aln)
    bam_writer.close()

    if options.verbose:
        sys.stdout.write(
            "Done {}".format(
                os.linesep
            )
        )

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
