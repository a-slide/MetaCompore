# -*- coding: utf-8 -*-

##### Imports #####

import logging
import numpy as np
import h5py
import pandas as pd
from collections import Counter
import pyfaidx
from scipy.signal import find_peaks
import datetime

##### DEFINE SCRIPT FUNCTION #####
def tombo_postprocess (
    res_h5,
    fasta,
    res_tsv,
    log,
    min_cov=50,
    p_val_lim=0.01,
    quantile_lim=0.5,
    min_distance=9):

    logging.basicConfig(filename=log, filemode="w", level=logging.INFO, format='%(message)s')
    logging.info("timestamp: {}".format(str(datetime.datetime.now())))
    for i, j in locals().items():
        logging.info("\t{}: {}\n".format(i,j))

    # Define variables
    sig_lim = -np.log10(p_val_lim)
    min_pval = np.nextafter(float(0), float(1))
    max_pval = 1

    # Init collections
    tx_id_set=set()
    c = Counter()
    # first = True

    # Get transcript lengths in dict for convenience
    logging.info(f'Load Fasta reference lengths\n')
    with pyfaidx.Fasta(fasta) as fa:
        tx_len_dict = {i.name:len(i) for i in fa}

    logging.info(f'Extract data from hdf5 database\n')
    with h5py.File(res_h5,'r') as h5, open(res_tsv, "w") as res_fp:
        res_fp.write("ref_id\tpos\tpvalue\tpeak\n")

        for block_id, block_data in h5["Statistic_Blocks"].items():

            # Extract attrs
            tx_id = block_data.attrs['chrm']
            start = block_data.attrs['start']
            strand = block_data.attrs['strand']

            if tx_id in tx_id_set:
                c["Duplicated transcript"]+=1
            elif start > 0:
                c["Transcript with invalid start"]+=1
            elif strand != "+":
                c["Transcript with invalid strand"]+=1
            else:
                tx_df = pd.DataFrame(block_data.get("block_stats")[()])
                tx_df = tx_df.dropna()
                tx_df = tx_df[(tx_df["cov"]>=min_cov) & (tx_df["control_cov"]>=min_cov)]
                if tx_df.empty:
                    c["Low coverage transcripts discarded"]+=1
                    continue

                tx_df.rename(columns={"stat":"pvalue"}, inplace=True)
                tx_df["pvalue"] = tx_df["pvalue"].fillna(1)
                tx_df["pvalue"] = np.clip(tx_df["pvalue"], min_pval, max_pval)

                # Peak calling in -log10 space
                c["All transcripts"]+=1
                x = pd.Series(data=-np.log10(tx_df["pvalue"]).values, index=tx_df["pos"].values)
                x = x.reindex(range(tx_len_dict[tx_id]))
                x = x.fillna(0)
                sig_val = x[x>sig_lim]

                if sig_val.empty:
                    c["Transcripts without significant pvalues"]+=1
                else:
                    c["Transcripts with significant pvalues"]+=1
                    threshold = np.quantile(sig_val, quantile_lim)
                    peaks = find_peaks(x, height=threshold, distance=min_distance)[0]

                    # Write significant hits +- peaks
                    for i in tx_df.itertuples():
                        if i.pvalue <=p_val_lim:
                            c["Significant pvalues"]+=1
                            if i.pos in peaks:
                                c["Peaks detected"]+=1
                                peak = True
                            else:
                                peak = False
                            res_fp.write(f"{tx_id}\t{i.pos}\t{i.pvalue}\t{peak}\n")

                tx_id_set.add(tx_id)

    logging.info(f'Counts\n')
    for i, j in c.items():
        logging.info(f'\t{i}:{j}\n')

##### RUN SCRIPT FUNCTION #####
tombo_postprocess (
    res_h5=snakemake.input.res_h5,
    fasta=snakemake.input.fasta,
    res_tsv=snakemake.output.res_tsv,
    log=snakemake.log[0],
    min_cov=snakemake.params.opt.get("min_cov", 30),
    p_val_lim=snakemake.params.opt.get("p_val_lim", 0.01),
    quantile_lim=snakemake.params.opt.get("quantile_lim", 0.5),
    min_distance=snakemake.params.opt.get("min_distance", 9))
