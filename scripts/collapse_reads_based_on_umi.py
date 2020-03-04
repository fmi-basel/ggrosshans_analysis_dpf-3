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

    __doc__ = "Collapse reads based on UMIs"

    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--fastq", dest="fastq", help="Input fastq file", required=True, metavar="FILE"
    )

    parser.add_argument(
        "--out",
        dest="out",
        help="Output directory",
        required=True,
        metavar="FILE",
    )

    parser.add_argument(
        "--compress",
        action="store_true",
        dest="compress",
        default=False,
        required=False,
        help="gzip compress output fastq files [Default: False]",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        required=False,
        help="Verbose [Default: False]",
    )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # get the arguments
    # -------------------------------------------------------------------------
    try:
        options = parser.parse_args()
    except (Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    if options.verbose:
        sys.stdout.write(f"Creating output directory: {options.out} {os.linesep}")
   
    if not os.path.exists(options.out):
        os.makedirs(options.out)
        
    # output files
    output_fastq = os.path.join(options.out, "collapsed.fastq")
    output_fastq_gz = output_fastq + ".gz"
    stats_out = os.path.join(options.out, "stats.tsv")
    stats_out_ignore_reads = os.path.join(options.out, "stats.ignore_reads.tsv")

    if options.verbose:
        sys.stdout.write(f"Parsing fastq file: {options.fastq} and writing fastq file: {output_fastq} {os.linesep}")

    # dictionary
    # key: sequence, value: list of umis []
    seq_to_list_of_umis = {}
    ignored_reads = 0

    with open(output_fastq, "w") as wp:

        for rec_id, rec_seq, rec_qual in FastqGeneralIterator(options.fastq):
            # determine umi
            rec_id_sp_space = rec_id.split(" ")
            rec_id_sp = rec_id_sp_space[0].split(":")
            rec_id = ":".join(rec_id_sp[0:-1]) + " " + rec_id_sp_space[1]
            umi_seq = rec_id_sp[-1]
            # If the read sequence was not observed before
            # write out the read entry and keep the observed
            # UMI in the list of UMIs for this sequence
            if rec_seq not in seq_to_list_of_umis:
                seq_to_list_of_umis[rec_seq] = [umi_seq]
                wp.write(f"@{rec_id}{os.linesep}")
                wp.write(f"{rec_seq}{os.linesep}")
                wp.write(f"+{os.linesep}")
                wp.write(f"{rec_qual}{os.linesep}")
            else:
                # If the read sequence was observed before,
                # but the UMI is not present in the list of
                # UMIs for this sequence then write out the read
                # entry and add the UMI in the list of UMIs
                if umi_seq not in seq_to_list_of_umis[rec_seq]:
                    seq_to_list_of_umis[rec_seq].append(umi_seq)
                    wp.write(f"@{rec_id}{os.linesep}")
                    wp.write(f"{rec_seq}{os.linesep}")
                    wp.write(f"+{os.linesep}")
                    wp.write(f"{rec_qual}{os.linesep}")
                # If the read sequence was observed before
                # and the UMI for this sequence is already in 
                # the list of UMIs for this sequence
                # then ignore this read
                else:
                    ignored_reads+=1

    if options.verbose:
        sys.stdout.write(f"Writing out stat files: {stats_out} and {stats_out_ignore_reads}{os.linesep}")
        
    stats = {}   
    for key, value in seq_to_list_of_umis.items():
        stats[key] = len(value)
    
    df = pd.DataFrame.from_dict(data=stats, orient="index")
    df.reset_index(inplace=True)
    df.columns = ["sequence", "number_of_umis"]
    df.to_csv(stats_out, header=True, sep="\t", index=False)
    
    with open(stats_out_ignore_reads, 'w') as wp:
        wp.write(f"Number of reads ignored: {ignored_reads}{os.linesep}")

    if options.compress:
        if options.verbose:
            sys.stdout.write(f"Compressing {output_fastq} -> {output_fastq_gz}{os.linesep}")
        with open(output_fastq, 'rb') as f_in, gzip.open(output_fastq_gz, 'wb') as f_out:
            f_out.writelines(f_in)

        if options.verbose:
            sys.stdout.write(f"Removing: {output_fastq}{os.linesep}")
        os.remove(output_fastq)

    if options.verbose:
        sys.stdout.write(f"Done {os.linesep}")


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(0)
