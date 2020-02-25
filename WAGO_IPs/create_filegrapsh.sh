#!/urs/bin/bash

snakemake --filegraph --configfile config.yaml | dot -Tpng > filegraph.png
