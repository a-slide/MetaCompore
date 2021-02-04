# MetaCompore v0.1.0

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.30.1-brightgreen.svg)](https://snakemake.bitbucket.io)

---

**Metacompore is a snakemake pipeline Running multiple RNA modifications detection tools for RNA modification detection**

## Authors

* Adrien Leger (@a-slide)

## Usage

### Step 1: Obtain a copy of this workflow

Clone the last tarball archive of the pipeline to your local system, into the place where you want to perform the data analysis

```
wget https://github.com/a-slide/MetaCompore/releases/download/0.0.3/MetaCompore.tar.gz
```

### Step 2: Install dependencies

#### Singularity

If required, install singularity following the official documentation: https://sylabs.io/guides/3.7/user-guide/quick_start.html

#### Conda / Mamba

Install miniconda following the official documentation: https://docs.conda.io/en/latest/miniconda.html

you can also install mamba to speed up snakemake installation: https://github.com/mamba-org/mamba

#### Snakemake

Create a virtual environment containing snakemake with [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or
```
conda env create -f environment.yaml
```

You can also use [mamba](https://github.com/mamba-org/mamba) which will give you the same result, but much faster

```
mamba env create -f environment.yaml
```

### Step 3: Configure the workflow

Configure the workflow according to your needs by editing the files `config.yaml` to configure the workflow execution,

```
nano config.yaml
```

Edit the `samples.tsv` to specify your sample setup and fast5 source files

```
nano samples.tsv
```

### Step 4: Execute workflow

### Local Mode

Activate the conda environment:

```
conda activate snakemake
snakemake --use-singularity -j 4
```

### LSF cluster Mode

Set an LSF cluster profile https://github.com/Snakemake-Profiles/lsf

Edit the lsf rule specific config file `lsf.yaml`



## Disclaimer

Please be aware that MetaCompore is a research package that is still under development.

It was tested under Linux Ubuntu 16.04 and in an HPC environment running under Red Hat Enterprise 7.1.

Thank you

## citation

Adrien Leger. a-slide/MetaCompore. Zenodo. [Not yet available]

## licence

MIT (https://mit-license.org/)

Copyright © 2020 Adrien Leger

## Authors

* Adrien Leger / contact@adrienleger.com / https://adrienleger.com
