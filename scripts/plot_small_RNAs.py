#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------


import sys
import os
import HTSeq
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
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
        "--DE_edgeR_table",
        dest="DE_edgeR_table",
        help="Provide DE edgeR table",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--comparison_name",
        dest="comparison_name",
        help="Comparison name to add in the title of the figure",
        required=True,
        metavar="STRING"
    )

    parser.add_argument(
        "--mirnas",
        dest="mirnas",
        help="Fasta file of annotated miRNAs",
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
        sys.stdout.write(
            "Creating output directory: {} {}".format(
                str(options.out),
                os.linesep
            )
        )

    if not os.path.exists(options.out):
        os.makedirs(options.out)

    if options.verbose:
        sys.stdout.write(
            "Reading edgeR table: {} {}".format(
                str(options.DE_edgeR_table),
                os.linesep
            )
        )

    if options.verbose:
        sys.stdout.write(
            "Reading miRNAs fasta file: {} {}".format(
                str(options.mirnas),
                os.linesep
            )
        )

    fa = HTSeq.HTSeq.FastaReader(options.mirnas)

    mirnas = {}
    for i in fa:
        mirna = i.seq.decode('UTF8')
        if mirna not in mirnas:
            mirnas[mirna] = mirna

    df = pd.read_csv(options.DE_edgeR_table,
                     header=0,
                     sep="\t")

    if options.verbose:
        sys.stdout.write(
            "Subseting tables for small RNAs (22G, 21U and 26G) and writing the to files {}".format(
                os.linesep
            )
        )

    df_22G = df[(df["id"].str.len() == 22) & (df["id"].str.startswith("G")) & (df["id"].isin(list(mirnas.keys())) == False)].copy()
    df_22G["Significant"] = "No"
    df_22G.loc[df_22G["FDR"]<0.05, "Significant"] = "FDR<0.05"
    df_22G.to_csv(os.path.join(options.out, "22G.tsv"), header=True, index=False, sep="\t")

    df_21U = df[(df["id"].str.len() == 21) & (df["id"].str.startswith("T")) & (df["id"].isin(list(mirnas.keys())) == False)].copy()
    df_21U["Significant"] = "No"
    df_21U.loc[df_21U["FDR"]<0.05, "Significant"] = "FDR<0.05"
    df_21U.to_csv(os.path.join(options.out, "21U.tsv"), header=True, index=False, sep="\t")

    df_26G = df[(df["id"].str.len() == 26) & (df["id"].str.startswith("G")) & (df["id"].isin(list(mirnas.keys())) == False)].copy()
    df_26G["Significant"] = "No"
    df_26G.loc[df_26G["FDR"]<0.05, "Significant"] = "FDR<0.05"
    df_26G.to_csv(os.path.join(options.out, "26G.tsv"), header=True, index=False, sep="\t")

    df_mirnas = df[df["id"].isin(list(mirnas.keys()))].copy()
    df_mirnas["Significant"] = "No"
    df_mirnas.loc[df_mirnas["FDR"]<0.05, "Significant"] = "FDR<0.05"
    df_mirnas.to_csv(os.path.join(options.out, "miRNAs.tsv"), header=True, index=False, sep="\t")

    if options.verbose:
        sys.stdout.write(
            "Genetating MA plots {}".format(
                os.linesep
            )
        )

    sns.set(font_scale=2)
    sns_plot = sns.lmplot(x="logCPM",
                          y="logFC",
                          data=df_22G,
                          fit_reg=False,
                          hue="Significant",
                          hue_order = ["No", "FDR<0.05"],
                          palette=["black", "red"],
                          height=15,
                          aspect=1)
    ax = plt.gca()
    ax.set_title("22G " + options.comparison_name)
    sns_plot.savefig(os.path.join(options.out, "22G.png"))
    sns_plot.savefig(os.path.join(options.out, "22G.pdf"))

    sns_plot = sns.lmplot(x="logCPM",
                          y="logFC",
                          data=df_21U,
                          fit_reg=False,
                          hue="Significant",
                          hue_order = ["No", "FDR<0.05"],
                          palette=["black", "red"],
                          height=15,
                          aspect=1)
    ax = plt.gca()
    ax.set_title("21U " + options.comparison_name)
    sns_plot.savefig(os.path.join(options.out, "21U.png"))
    sns_plot.savefig(os.path.join(options.out, "21U.pdf"))

    sns_plot = sns.lmplot(x="logCPM",
                          y="logFC",
                          data=df_26G,
                          fit_reg=False,
                          hue="Significant",
                          hue_order = ["No", "FDR<0.05"],
                          palette=["black", "red"],
                          height=15,
                          aspect=1)
    ax = plt.gca()
    ax.set_title("26G " + options.comparison_name)
    sns_plot.savefig(os.path.join(options.out, "26G.png"))
    sns_plot.savefig(os.path.join(options.out, "26G.pdf"))

    sns_plot = sns.lmplot(x="logCPM",
                          y="logFC",
                          data=df_mirnas,
                          fit_reg=False,
                          hue="Significant",
                          hue_order = ["No", "FDR<0.05"],
                          palette=["black", "red"],
                          height=15,
                          aspect=1)
    ax = plt.gca()
    ax.set_title("miRNAs " + options.comparison_name)
    sns_plot.savefig(os.path.join(options.out, "miRNAs.png"))
    sns_plot.savefig(os.path.join(options.out, "miRNAs.pdf"))

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
