#!/urs/bin/bash

set -e

mkdir -p logs/cluster_log
mkdir -p logs/local_log

snakemake \
--use-conda \
--configfile config.yaml \
-p \
--local-cores 40 \
--cores 40 \
--rerun-incomplete
