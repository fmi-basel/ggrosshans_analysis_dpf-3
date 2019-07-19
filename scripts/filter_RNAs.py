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
        "--gtf",
        dest="gtf",
        help="Annotation file in gtf format",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--bam",
        dest="bam",
        help="Alignment file in BAM (sorted and indexed) format. Output of bowtie.",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--filter_RNAs",
        dest="filter_RNAs",
        help="List of RNAs to filter",
        required=True,
        metavar="STRING"
    )

    parser.add_argument(
        "--filter_chromosomes",
        dest="filter_chromosomes",
        help="List of chromosomes to filter",
        required=True,
        metavar="STRING"
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

    filter_RNAs = str(options.filter_RNAs).split()
    filter_chromosomes = str(options.filter_chromosomes).split()

    # parse gtf and bam file
    if options.verbose:
        sys.stdout.write(
            "Parsing gtf file: {} and fetching reads from bam file to filter RNAs by biotype: {} {}".format(
                str(options.gtf),
                str(options.bam),
                os.linesep
            )
        )

    # read names that fall in the "filtered RNAs" or
    reads_names_to_exclude = {}

    gtf = HTSeq.GFF_Reader(options.gtf)
    bam = HTSeq.BAM_Reader(options.bam)

    for line in gtf:
        chrom = str(line.iv.chrom)
        start = str(line.iv.start + 1)
        end = str(line.iv.end)
        strand = str(line.iv.strand)
        transcript_biotype = line.attr["transcript_biotype"]

        if transcript_biotype in filter_RNAs:
            region_string = chrom + ":" + start + "-" + end
            for aln in bam.fetch(region=region_string):
                if aln.iv.strand == strand:
                    read_name = aln.read.name
                    if read_name not in reads_names_to_exclude:
                        reads_names_to_exclude[read_name] = read_name

    # parse bam file based on chromosomes
    if options.verbose:
        sys.stdout.write(
            "Fetching reads from bam file to filter RNAs based on chromosomes: {} {}".format(
                str(options.gtf),
                str(options.bam),
                os.linesep
            )
        )

    for filter_chromosome in filter_chromosomes:
        for aln in bam.fetch(filter_chromosome):
            chrom = str(line.iv.chrom)
            read_name = aln.read.name
            if read_name not in reads_names_to_exclude:
                reads_names_to_exclude[read_name] = read_name

    # parse gtf and bam file
    if options.verbose:
        sys.stdout.write(
            "Filtering reads from bam file: {} and writing output to: {} {}".format(
                str(options.bam),
                str(options.bam_out),
                os.linesep
            )
        )

    bam_writer = HTSeq.BAM_Writer.from_BAM_Reader(options.bam_out, bam)
    for aln in bam:
        if (aln.aligned) & (aln.read.name not in reads_names_to_exclude):
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
