#!/urs/bin/bash

snakemake --rulegraph --configfile config.yaml | dot -Tpng > rulegraph.png
