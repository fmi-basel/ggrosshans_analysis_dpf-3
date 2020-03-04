#!/urs/bin/bash

snakemake --use-conda --configfile config.yaml -np --dag | dot -Tpng > dag.png
