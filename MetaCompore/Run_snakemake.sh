#!/bin/bash
set -euo pipefail

snakemake $* --use-singularity -j 4 --keep-going
