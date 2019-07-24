#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------


import sys
import os
import HTSeq
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
        "--gff",
        dest="gff",
        help="Annotation file (transposons) in gff format",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--bed",
        dest="bed",
        help="Output file in bed format",
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

    # parse gtf and bam file
    if options.verbose:
        sys.stdout.write(f"Parsing gtf file: {options.gff} {os.linesep}")
        sys.stdout.write(f"and writiing bed file: {options.bed} {os.linesep}")

    gff_file= HTSeq.GFF_Reader(options.gff)

    w = open(options.bed, "w")
    for gff_line in gff_file:
        chrom = gff_line.iv.chrom
        start = str(gff_line.iv.start)
        end = str(gff_line.iv.end)
        strand = gff_line.iv.strand
        score = str(0)
        name = gff_line.attr['Name'] + "$" +gff_line.attr['Family']
        w.write("\t".join([chrom, start, end, name, score, strand + os.linesep]))
    w.close()

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
