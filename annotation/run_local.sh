#!/urs/bin/bash

set -e

mkdir -p logs/cluster_log
mkdir -p logs/local_log

snakemake \
--use-conda \
--configfile config.yaml \
-p \
--local-cores 60 \
--cores 60 \
--rerun-incomplete
