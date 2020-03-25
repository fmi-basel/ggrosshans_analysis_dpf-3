#!/urs/bin/bash

set -e

mkdir -p logs/cluster_log
mkdir -p logs/local_log

snakemake \
--use-conda \
--conda-prefix ../conda \
--configfile config.yaml \
-p \
--local-cores 50 \
--cores 50 \
--rerun-incomplete
