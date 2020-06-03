# Analysis dpf-3

In the current repository you can find the pipelines used to process the high throughput sequencing libraries (small RNA-seq, total RNA-seq, small RNA-seq IPs) and scripts to generate the final plots for the manuscript:

```
Protease-mediated processing of Argonaute proteins controls small RNA association

Rajani kanth Gudipati, Kathrin Braun, Foivos Gypas, Daniel Hess, Jan Schreier, Sarah H. Carl, René F. Ketting, Helge Großhans
```

# Installation and execution

The repo consists of two parts:
- Snakemake pipelines that preprocess the data
- Jupyter notebooks (and saannome R script) to generate the final plots

## Install snakemake and conda

In order to re-run the pipelines you need to install snakemake and conda/miniconda.

#### Step 1: Download miniconda 3 installation file (if not already installed)

You can do this with one of the following options:

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

#### Step 2: Install miniconda 3

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

### Step 3: Create a new conda environment

Create a new conda environment with snakemake
```bash
conda create --name snakemake --channel bioconda --channel conda-forge snakemake=5.11.1
```

## Run the pipelines

Download the necessary datasets from GEO/SRA and update the samples.tsv files in the following directories:
- 00_annotation
- 01_small_RNA_seq_15_C
- 02_total_RNA_seq_15_C
- 03_WAGO_IPs
- 04_csr_1_IP

Activate the conda environment containing snakemake
```bash
conda activate snakemake
```

Run the script
```bash
bash run_local.sh
```

## Install dependencies and generate plots 

Once you have the results from the pipelines you can generate the plots

Create a new conda environment to run the plotting scripts
```bash
conda env create -f envs/environment.yml
```

Activate the environment
```bash
conda activate dpf_3_python
```

Then you can open the jupyter notebooks (e.g. with Jupyter Lab) and execute them.