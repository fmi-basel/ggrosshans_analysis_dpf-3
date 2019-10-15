#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------


import sys
import os
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
        "--ucsc_repeats",
        dest="ucsc_repeats",
        help="Ucsc repeats file (tsv)",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--gtf",
        dest="gtf",
        help="Output gtf file",
        required=True,
        metavar="FILE"
    )
    
    parser.add_argument(
        "--bed",
        dest="bed",
        help="Output bed file",
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
        sys.stdout.write(f"Parsing file: {options.ucsc_repeats} and writing file: {options.gtf} {os.linesep}")

    df = pd.read_csv(options.ucsc_repeats, header=0, sep="\t")
    
    w = open(options.gtf, 'w')
    for index, row in df.iterrows():
        
        chrom = str(row['genoName'])
        start = str(row['genoStart']+1)
        end = str(row['genoEnd'])
        strand = str(row["strand"])
        
        gene_id = f"{chrom}:{start}-{end}:{strand};"
        repName = f"repName \"{row['repName']}\";"
        repClass = f"repClass \"{row['repClass']}\";"
        repFamily = f"repFamily \"{row['repFamily']}\";"
        
        attributes = " ".join([gene_id, repName, repClass, repFamily + os.linesep])
        
        w.write("\t".join([chrom,
                          "uscs_repbase",
                          "exon",
                          str(start),
                          str(end),
                          ".",
                          strand,
                          ".",
                          attributes]))
    w.close()
    
    w = open(options.bed, 'w')
    for index, row in df.iterrows():
        
        chrom = str(row['genoName'])
        start = str(row['genoStart'])
        end = str(row['genoEnd'])
        strand = str(row["strand"])
        
        gene_id = f"{chrom}:{start}-{end}:{strand}"
        # repName = f"repName \"{row['repName']}\";"
        # repClass = f"repClass \"{row['repClass']}\";"
        # repFamily = f"repFamily \"{row['repFamily']}\";"
        
        attributes = "|".join([gene_id, row['repName'], row['repClass'], row['repFamily']])
        
        w.write("\t".join([chrom,
                           str(start),
                           str(end),
                           attributes,
                           str(0),
                           strand + os.linesep]))
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
