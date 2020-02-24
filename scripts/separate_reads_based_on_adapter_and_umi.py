#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------


import sys
import os
import gzip
import pandas as pd
import numpy as np
from Bio.SeqIO.QualityIO import FastqGeneralIterator
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
        "--umi-length",
        dest="umi_length",
        help="UMI length",
        required=True,
        metavar="STR"
    )

    parser.add_argument(
        "--out",
        dest="out",
        help="Output directory",
        required=True,
        metavar="FILE"
    )
    
    parser.add_argument(
        "--compress",
        action="store_true",
        dest="compress",
        default=False,
        required=False,
        help="gzip compress output fastq files [Default: False]"
    )
    
    parser.add_argument(
        "--attach-umis-to-header",
        action="store_true",
        dest="attach_umis_to_header",
        default=False,
        required=False,
        help="Attach umis to header [Default: False]"
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        required=False,
        help="Verbose [Default: False]"
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
        sys.stdout.write(f"Creating output directory: {options.out} {os.linesep}")
   
    if not os.path.exists(options.out):
        os.makedirs(options.out)
        
    # output files
    adapter_missing_fastq = os.path.join(options.out, "adapter_missing.fastq")
    umi_truncated_fastq = os.path.join(options.out, "umi_truncated.fastq")
    adapter_and_umi_present_fastq = os.path.join(options.out, "adapter_and_umi_present.fastq")
    stats_out = os.path.join(options.out, "stats.tsv")
    adapter_missing_fastq_gz = adapter_missing_fastq + ".gz"
    umi_truncated_fastq_gz = umi_truncated_fastq + ".gz"
    adapter_and_umi_present_fastq_gz = adapter_and_umi_present_fastq + ".gz"
        
    if options.verbose:
        sys.stdout.write(f"Parsing fastq file: {options.fastq}{os.linesep}and separating reads in different output files: {os.linesep} - {umi_truncated_fastq}{os.linesep} - {adapter_and_umi_present_fastq}{os.linesep} - {adapter_and_umi_present_fastq}{os.linesep}")
    
    adapter_len = len(str(options.adapter))

    stats = {'adapter_missing': 0, 'umi_truncated': 0, 'adapter_and_umi_present':0}

    with open(adapter_missing_fastq, 'w') as adapter_missing, \
         open(umi_truncated_fastq, 'w') as umi_truncated, \
         open(adapter_and_umi_present_fastq, 'w') as adapter_and_umi_present:

        for rec_id, rec_seq, rec_qual in FastqGeneralIterator(options.fastq):
            # detect match (-1 in case of not match)
            first_position = (str(rec_seq).find(str(options.adapter)))
            # adapter missing
            if first_position == -1:
                stats['adapter_missing']+=1
                adapter_missing.write(f"@{rec_id}{os.linesep}")
                adapter_missing.write(f"{rec_seq}{os.linesep}")
                adapter_missing.write(f"+{os.linesep}")
                adapter_missing.write(f"{rec_qual}{os.linesep}")
            # adapter present
            else:
                sequence_after_adapter = str(rec_seq[first_position+adapter_len:])
                # umi (or part) missing
                if len(sequence_after_adapter)<int(options.umi_length):
                    stats['umi_truncated']+=1
                    umi_truncated.write(f"@{rec_id}{os.linesep}")
                    umi_truncated.write(f"{rec_seq}{os.linesep}")
                    umi_truncated.write(f"+{os.linesep}")
                    umi_truncated.write(f"{rec_qual}{os.linesep}")
                # umi present
                else:
                    stats['adapter_and_umi_present']+=1
                    if options.attach_umis_to_header:
                        umi = sequence_after_adapter[0:int(options.umi_length)]
                        rec_id_sp = rec_id.split(" ")
                        rec_id_sp[0] = ":".join([rec_id_sp[0], umi])
                        rec_id = " ".join(rec_id_sp)

                    adapter_and_umi_present.write(f"@{rec_id}{os.linesep}")
                    adapter_and_umi_present.write(f"{rec_seq}{os.linesep}")
                    adapter_and_umi_present.write(f"+{os.linesep}")
                    adapter_and_umi_present.write(f"{rec_qual}{os.linesep}")

    if options.verbose:
        sys.stdout.write(f"Writing stats file: {stats_out}  {os.linesep}")    

    df = pd.DataFrame.from_dict(data=stats, orient="index")
    df.reset_index(inplace=True)
    df.columns = ["category", "number_of_reads"]
    df.to_csv(stats_out, header=True, sep="\t", index=False)
    
    if options.compress:
        if options.verbose:
            sys.stdout.write(f"Compressing {adapter_missing_fastq} -> {adapter_missing_fastq_gz}{os.linesep}")
        with open(adapter_missing_fastq, 'rb') as f_in, gzip.open(adapter_missing_fastq_gz, 'wb') as f_out:
            f_out.writelines(f_in)

        if options.verbose:
            sys.stdout.write(f"Compressing {umi_truncated_fastq} -> {umi_truncated_fastq_gz}{os.linesep}")
        with open(umi_truncated_fastq, 'rb') as f_in, gzip.open(umi_truncated_fastq_gz, 'wb') as f_out:
            f_out.writelines(f_in)
            
        if options.verbose:
            sys.stdout.write(f"Compressing {adapter_and_umi_present_fastq} -> {adapter_and_umi_present_fastq_gz}{os.linesep}")
        with open(adapter_and_umi_present_fastq, 'rb') as f_in, gzip.open(adapter_and_umi_present_fastq_gz, 'wb') as f_out:
            f_out.writelines(f_in)
        
        if options.verbose:
            sys.stdout.write(f"Removing: {adapter_missing_fastq}, {umi_truncated_fastq}, {adapter_and_umi_present_fastq}  {os.linesep}")
        os.remove(adapter_missing_fastq)
        os.remove(umi_truncated_fastq)
        os.remove(adapter_and_umi_present_fastq)

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
