#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------


import sys
import os
import pandas as pd
from Bio import SeqIO
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
        "--fastq",
        dest="fastq",
        help="Input fastq file",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--adapter",
        dest="adapter",
        help="3' adapter",
        required=True,
        metavar="STR"
    )
    
    parser.add_argument(
        "--extend-after-adapter",
        dest="extend_after_adapter",
        help="Extend after 3' adapter",
        required=True,
        metavar="INT"
    )

    parser.add_argument(
        "--out",
        dest="out",
        help="Output file",
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
        sys.stdout.write(f"Searcing fastq file for potential UMI sequences: {options.fastq} {os.linesep}")
    
    potential_umis = {}
    
    adapter_len = len(str(options.adapter))
    for rec in SeqIO.parse(options.fastq, "fastq"):
        
        # detect match (-1 in case of not match)
        first_position = (str(rec.seq).find(str(options.adapter)))
        if first_position != -1:
            
            potential_umi = str(rec.seq[first_position+adapter_len:first_position+adapter_len+int(options.extend_after_adapter)])
            
            if len(potential_umi) != int(options.extend_after_adapter):
                continue
        
            if potential_umi not in potential_umis:
                potential_umis[potential_umi] = 1
            else:
                potential_umis[potential_umi] += 1
    
    df = pd.DataFrame.from_dict(data=potential_umis, orient="index")
    df.reset_index(inplace=True)
    df.columns = ["potential_umi", "count"]
    df.sort_values("count", inplace=True, ascending=False)
    df.to_csv(options.out, header=True, sep="\t", index=False)
         
    if options.verbose:
        sys.stdout.write("Done {os.linesep}")

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
