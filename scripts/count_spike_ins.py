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

    __doc__ = "Count spike in occurences"

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )
    
    parser.add_argument(
        "--fasta",
        dest="fasta",
        help="Fasta file with the spike-ins.",
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
        "--out",
        dest="out",
        help="Output directory",
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
        sys.stdout.write(f"Creating output directory: {str(options.out)} {os.linesep}")
        sys.stdout.flush()

    if not os.path.exists(options.out):
        os.makedirs(options.out)
        
    # parse fasta file and keep header names
    if options.verbose:
        sys.stdout.write(f"Parsing fasta file: {str(options.fasta)} {os.linesep}")
        sys.stdout.flush()
        
    # counts dictionary
    # key sequence name
    # value number of reads
    counts = {}
    
    fasta = HTSeq.FastaReader(options.fasta)
    
    for entry in fasta:
        name = entry.name
        counts[name] = 0

    # parse alignment file and count
    if options.verbose:
        sys.stdout.write(f"Parsing bam file: {str(options.bam)} {os.linesep}")
        sys.stdout.flush()

    bam = HTSeq.BAM_Reader(options.bam)
    for read in bam:
        if read.aligned:
            # ignore reads mapping to the - strand
            if read.iv.strand == "-":
               continue
            # find spike in id
            spike_in_id = str(read.iv.chrom)
            # find to how many loci it maps
            mapped_loci = read.optional_field('XM') - 1
    	    # ignore muntimappers
            #if mapped_loci>1:
            #    continue

            if spike_in_id in counts:
                counts[spike_in_id] += 1/mapped_loci
    
    print(counts)

    if options.verbose:
        sys.stdout.write(f"Creating pandas dataframe from counts dictionary {os.linesep}")
        sys.stdout.flush()

    df = pd.DataFrame.from_dict(counts, orient='index')
    df.reset_index(inplace=True)
    df.columns = [['Name', 'counts']]

    if options.verbose:
        sys.stdout.write(f"Writing files in output directory: {str(options.out)} {os.linesep}")

    df.to_csv(os.path.join(options.out, "counts.tsv"),
              sep="\t",
              header=True,
              index=False
    )

    if options.verbose:
        sys.stdout.write(f"Done {os.linesep}")
        sys.stdout.flush()

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
