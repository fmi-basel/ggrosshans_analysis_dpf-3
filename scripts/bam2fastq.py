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

    __doc__ = "Convert bowtie alignment file (BAM) to fastq. Multimappers are kept only once."

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
        "--fastq",
        dest="fastq",
        help="Output file in fastq format.",
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

    # parse alignment file and count
    if options.verbose:
        sys.stdout.write(
            "Parsing bam file: {} and writing fastq file: {} {}".format(
                str(options.bam),
                str(options.fastq),
                os.linesep
            )
        )

    # observed_reads
    observed_reads = {}

    w = open(options.fastq, 'w')
    bam = HTSeq.BAM_Reader(options.bam)
    for read in bam:
        if read.aligned:
            read_name = read.read.name
            read_sequence = read.read.seq.decode('utf-8')
            read_quality = read.read.qualstr.decode('utf-8')            
            if read_name not in observed_reads:
                observed_reads[read_name] = read_name
                w.write("@" + read_name + os.linesep)
                w.write(read_sequence + os.linesep)
                w.write("+" + os.linesep)
                w.write(read_quality + os.linesep)
    w.close()

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