# -*- coding: utf-8 -*-

##### Imports #####

import logging
import numpy as np
import h5py
import pandas as pd
from collections import Counter
import pyfaidx
from scipy.signal import find_peaks

##### DEFINE SCRIPT FUNCTION #####
def tombo_postprocess (stats_h5, fasta, stats_tsv, log, min_cov=50, p_val_lim=0.01, quantile_lim=0.5, min_distance=9):
    """
    Extract data from tombo stat database, cleanup and denoise
    """
    logging.basicConfig(filename=log, filemode="w", level=logging.DEBUG)

    # Define variables
    sig_lim = -np.log10(p_val_lim)
    min_pval = np.nextafter(float(0), float(1))
    max_pval = 1

    # Init collections
    tx_id_list=[]
    c = Counter()
    first = True

    # Get transcript lengths in dict for convenience
    logging.debug(f'Load Fasta reference lengths\n')
    with pyfaidx.Fasta(fasta) as fa:
        tx_len_dict = {i.name:len(i) for i in fa}

    logging.debug(f'Extract data from hdf5 database\n')
    with h5py.File(stats_h5,'r') as h5, open(stats_tsv, "w") as tsv_fp:
        for block_id, block_data in h5["Statistic_Blocks"].items():
            tx_id = block_data.attrs['chrm']

            if tx_id in tx_id_list:
                c["Duplicated transcript"]+=1
                continue

            start = block_data.attrs['start']
            assert start == 0
            strand = block_data.attrs['strand']
            assert strand == "+"

            tx_len = tx_len_dict[tx_id]

            block_df = pd.DataFrame(block_data.get("block_stats")[()])
            block_df = block_df.dropna()
            block_df = block_df[(block_df["cov"]>=min_cov) & (block_df["control_cov"]>=min_cov)]
            if block_df.empty:
                c["Low coverage transcripts discarded"]+=1
                continue

            block_df.rename(columns={"stat":"pvalue", "cov":"test_cov"}, inplace=True)
            block_df["pvalue"] = np.clip(block_df["pvalue"], min_pval, max_pval)
            block_df["ref_id"] = tx_id
            block_df["log_pvalue"] = -np.log10(block_df["pvalue"])
            c["Valid transcripts"]+=1

            # Peak calling in -log10 space
            log_pvalue = pd.Series(data=block_df["log_pvalue"].values, index=block_df["pos"].values)
            log_pvalue = log_pvalue.reindex(range(tx_len))
            log_pvalue = log_pvalue.fillna(0)
            sig_val = log_pvalue[log_pvalue>sig_lim]
            c["Significant pvalues"]+=len(sig_val)

            if sig_val.empty:
                peaks=[]
                block_df["peak"] = False
            else:
                threshold = np.quantile(sig_val, quantile_lim)
                peaks, extra = find_peaks(log_pvalue, height=threshold, distance=min_distance)
                peaks_status = [i in peaks for i in block_df["pos"]]
                block_df["peak"] = peaks_status
                c["Peaks detected"] += len(peaks)

            tx_id_list.append(tx_id)
            block_df = block_df.reindex(columns=["ref_id", "pos", "pvalue", "log_pvalue", "peak", "test_cov","control_cov"])
            block_df.to_csv(tsv_fp, sep="\t", header=first, index=False)
            first=False

    logging.debug(f'Counts\n')
    for i, j in c.items():
        logging.debug(f'\t{i}:{j}\n')

##### RUN SCRIPT FUNCTION #####
tombo_postprocess (
    stats_h5=snakemake.input.stats_h5,
    fasta=snakemake.input.fasta,
    stats_tsv=snakemake.output.stats_tsv,
    log = snakemake.log[0],
    min_cov=snakemake.params.opt.get("min_cov", 60),
    p_val_lim=snakemake.params.opt.get("p_val_lim", 0.01),
    quantile_lim=snakemake.params.opt.get("quantile_lim", 0.5),
    min_distance=snakemake.params.opt.get("min_distance", 9))
