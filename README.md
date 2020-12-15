# MetaCompore

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.30.1-brightgreen.svg)](https://snakemake.bitbucket.io)

**Metacompore is a snakemake pipeline Running multiple RNA modifications detection tools for RNA modification detection**

## Authors

* Adrien Leger (@a-slide)

## Usage

### Step 1: Install dependencies

#### Singularity

If required, install singularity following the official documentation: https://sylabs.io/guides/3.7/user-guide/quick_start.html

#### Conda / Mamba

**This is optional if you chose option 3 to install snakemake**

Install miniconda following the official documentation: https://docs.conda.io/en/latest/miniconda.html

you can also install mamba to speed up snakemake installation: https://github.com/mamba-org/mamba

#### Install snakemake

**1) Conda option**

Create a virtual environment containing snakemake with [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or
```
conda env create -f environment.yaml
```

**2) Mamba option**

You can also use [mamba](https://github.com/mamba-org/mamba) which will give you the same result, but much faster

```
mamba env create -f environment.yaml
```

**3) Singularity option**

Alternatively, you can also pull a snakemake container with singularity

```
singularity pull --arch amd64 library://aleg/default/snakemake:5.30.1
```

### Step 2: Obtain a copy of this workflow

Clone this github repository to your local system, into the place where you want to perform the data analysis

```
git clone git@github.com:a-slide/MetaCompore.git {local_directory}
cd {local_directory}
```

### Step 3: Configure workflow

Configure the workflow according to your needs by editing the files `config.yaml` to configure the workflow execution,

```
nano config.yaml
```

Edit the `samples.tsv` to specify your sample setup and fast5 source files

```
nano samples.tsv
```

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

### Step 5: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.
An example (using some trivial test data) can be seen [here](https://cdn.rawgit.com/snakemake-workflows/rna-seq-kallisto-sleuth/master/.test/report.html).

### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 7: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/MetaCompore.git` or `git remote add -f upstream https://github.com/snakemake-workflows/MetaCompore.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


### Step 8: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.

## Testing

Test cases are in the subfolder `.test`. They are automatically executed via continuous integration with [Github Actions](https://github.com/features/actions).
