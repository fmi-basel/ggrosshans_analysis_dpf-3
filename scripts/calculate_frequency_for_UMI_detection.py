#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------


import sys
import os
import pandas as pd
import numpy as np
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
        sys.stdout.write(f"Extracting sequences after adapter {options.adapter} from fastq file {options.fastq} and calculating A/C/G/T frequency {os.linesep}")
    
    A = np.zeros(100)
    C = np.zeros(100)
    G = np.zeros(100)
    T = np.zeros(100)
    
    adapter_len = len(str(options.adapter))

    for rec in SeqIO.parse(options.fastq, "fastq"):
        
        # detect match (-1 in case of not match)
        first_position = (str(rec.seq).find(str(options.adapter)))
        if first_position != -1:
            
            sequence_after_adapter = str(rec.seq[first_position+adapter_len:])
            for i in range(0, len(sequence_after_adapter)):
                if sequence_after_adapter[i] == "A":
                    A[i]+=1
                elif sequence_after_adapter[i] == "C":
                    C[i]+=1
                elif sequence_after_adapter[i] == "G":
                    G[i]+=1
                elif sequence_after_adapter[i] == "T":
                    T[i]+=1
                elif sequence_after_adapter[i] == "N":
                    pass
                else:
                    sys.stderr.write(f"Weird sequence: {sequence_after_adapter[i]} detected {os.linesep}")
                    sys.exit(-1)

    if options.verbose:
        sys.stdout.write(f"Write out frequencies in: {options.out} {os.linesep}")

    df = pd.DataFrame(data=[A, C, G, T]).T
    df.columns = ["A", "C", "G", "T"]
    df.to_csv(options.out, header=True, sep="\t")
            
         
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
