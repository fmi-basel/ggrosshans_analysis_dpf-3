""" merge reads """

import argparse
import sys
import math
import os
import re
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter, FileType

###################

def main():
    """ Main """

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument('--input_condition1',
        dest='input_condition1',
        help='quantification files for condition1')

    parser.add_argument('--condition1_name',
        dest='condition1_name',
        help='condition1 name')

    parser.add_argument('--input_condition2',
        dest='input_condition2',
        help='quantification files for condition2')

    parser.add_argument('--condition2_name',
        dest='condition2_name',
        help='condition2 name')

    parser.add_argument('--outfile_counts',
        dest='outfile_counts',
        help='output_file')

    parser.add_argument('--outfile_conditions',
        dest='outfile_conditions',
        help='output_file')

    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    condition1_files = str(options.input_condition1).split()
    condition2_files = str(options.input_condition2).split()
    condition1 = options.condition1_name
    condition2 = options.condition2_name

    conditions =[]
    counts = pd.DataFrame()

    counter = 0
    for condition1_file in condition1_files:
        condition1_df = pd.read_csv(condition1_file,
            header=0,
            sep='\t',
            index_col=0,
            comment='#',
            engine='python')
        condition1_df['counts'] = condition1_df['counts'].round()
        condition1_df.counts = condition1_df.counts.astype(int)
        condition1_df.columns = [str(condition1) + "_" + str(counter)]
        conditions.append(condition1)
        if counts.empty:
            counts = condition1_df
        else:
            counts = counts.merge(condition1_df, how='outer', right_index=True, left_index=True)
        counter += 1
    counter = 0
    for condition2_file in condition2_files:
        condition2_df = pd.read_csv(condition2_file,
            header=0,
            sep='\t',
            index_col=0,
            comment='#',
            engine='python')
        condition2_df['counts'] = condition2_df['counts'].round()
        condition2_df.counts = condition2_df.counts.astype(int)
        condition2_df.columns = [str(condition2) + "_" + str(counter)]
        conditions.append(condition2)
        if counts.empty:
            counts = condition2_df
        else:
            counts = counts.merge(condition2_df, how='outer', right_index=True, left_index=True)
        counter += 1

    counts.fillna(0, inplace=True)
    counts.to_csv(options.outfile_counts, index=True)
    #conditions_file = open(options.outfile_conditions,'w')
    conditions=pd.DataFrame(conditions)
    conditions=conditions.T
    conditions.to_csv(options.outfile_conditions, sep=',', index=None, header=None)
    #conditions_file.write(",".join(conditions))
    #conditions_file.close()
    return
###################

if __name__ == '__main__':
    main()
