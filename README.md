# MetaCompore

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.30.1-brightgreen.svg)](https://snakemake.bitbucket.io)

**Metacompore is a snakemake pipeline Running multiple RNA modifications detection tools for RNA modification detection**

## Authors

* Adrien Leger (@a-slide)

## Usage

### Step 1: Obtain a copy of this workflow

Clone this github repository to your local system, into the place where you want to perform the data analysis

```
git clone git@github.com:a-slide/MetaCompore.git {local_directory}
cd {local_directory}
```

### Step 2: Install dependencies

#### Singularity

If required, install singularity following the official documentation: https://sylabs.io/guides/3.7/user-guide/quick_start.html

#### Conda / Mamba

**This is not required if you choose option 2 to install snakemake**

Install miniconda following the official documentation: https://docs.conda.io/en/latest/miniconda.html

you can also install mamba to speed up snakemake installation: https://github.com/mamba-org/mamba

#### Snakemake

**1) Conda/Mamba option**

Create a virtual environment containing snakemake with [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or
```
conda env create -f environment.yaml
```

You can also use [mamba](https://github.com/mamba-org/mamba) which will give you the same result, but much faster

```
mamba env create -f environment.yaml
```

**2) Singularity option**

Alternatively, you can also pull a snakemake container with singularity

```
singularity pull --arch amd64 library://aleg/default/snakemake:5.30.1
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

#### With Singularity

```
singularity run snakemake_5.30.1.sif --use-singularity -j 4
```

#### With conda

Activate the conda environment:

```
conda activate snakemake
snakemake --use-singularity -j 4
```

### Step 5: Investigate results
