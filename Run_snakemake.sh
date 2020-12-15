#!/bin/bash
set -euo pipefail

snakemake $* --use-singularity -j 4 --wms-monitor http://127.0.0.1:5000
