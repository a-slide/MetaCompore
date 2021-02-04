# -*- coding: utf-8 -*-

##### Imports #####

import logging
import numpy as np
import pandas as pd
from collections import Counter, OrderedDict
import pyfaidx
from scipy.signal import find_peaks
import datetime
import os

##### DEFINE SCRIPT FUNCTION #####
def nanocompore_postprocess (
    input_tsv,
    fasta,
    outdir,
    log,
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

    tests_d = {
        "GMM_logit_pvalue": {"peak":"GMM_logit_peak" ,"label":"GMM_context_0"},
        "GMM_logit_pvalue_context_2": {"peak":"GMM_logit_peak_context_2", "label":"GMM_context_2"},
        "KS_intensity_pvalue": {"peak":"KS_intensity_peak", "label":"KS_intensity_context_0"},
        "KS_intensity_pvalue_context_2": {"peak":"KS_intensity_peak_context_2", "label":"KS_intensity_context_2"},
        "KS_dwell_pvalue": {"peak":"KS_dwell_peak", "label":"KS_dwell_context_0"},
        "KS_dwell_pvalue_context_2": {"peak":"KS_dwell_peakontext_2", "label":"KS_dwell_context_2"}}

    # Get transcript lengths in dict for convenience
    logging.info('Load Fasta reference lengths\n')
    with pyfaidx.Fasta(fasta) as fa:
        tx_len_dict = {i.name:len(i) for i in fa}

    # Get data and cleanup
    logging.info('Load and cleanup data\n')
    df = pd.read_csv(input_tsv, sep="\t", dtype={'chr':str})

    test_sig_d = OrderedDict()
    logging.info('Iterate over transcripts and call peaks\n')

    for test in list(tests_d.keys()):
        c = Counter()
        test_label = tests_d[test]["label"]

        with open (os.path.join(outdir, f"nanocompore_results_{test_label}.tsv"), "w") as res_fp:
            res_fp.write("ref_id\tpos\tpvalue\tpeak\n")
            # Extract data for current test and cleanup
            test_df = df[["ref_id", "pos", test]].copy()
            test_df = test_df.rename(columns={test:"pvalue"})
            test_df["pvalue"]= test_df["pvalue"].fillna(1)
            test_df["pvalue"] = test_df["pvalue"].clip(min_pval, max_pval)

            # Iterate over tx for peak calling
            for tx, tx_df in test_df.groupby("ref_id"):
                c["All transcripts"]+=1
                x = pd.Series(data=-np.log10(tx_df["pvalue"]).values, index=tx_df["pos"].values)
                x = x.reindex(range(tx_len_dict[tx]))
                x = x.fillna(0)
                sig_val = x[x>=sig_lim]

                if sig_val.empty:
                    c["Transcripts without significant pvalues"]+=1
                else:
                    c["Transcripts with significant pvalues"]+=1
                    threshold = np.quantile(sig_val, quantile_lim)
                    peaks = find_peaks(x, height=threshold, distance=min_distance)[0]
                    if peaks.size == 0:
                        c["Transcripts without peaks called"]+=1
                    else:
                        c["Transcripts with peaks called"]+=1

                    # Write significant hits +- peaks
                    for i in tx_df.itertuples():
                        if i.pvalue <=p_val_lim:
                            c["Significant pvalues"]+=1
                            if i.pos in peaks:
                                c["Peaks detected"]+=1
                                peak = True
                            else:
                                peak = False
                            res_fp.write(f"{i.ref_id}\t{i.pos}\t{i.pvalue}\t{peak}\n")

        logging.info(f'{test_label} counts\n')
        for i, j in c.items():
            logging.info(f'\t{i}:{j}\n')

##### RUN SCRIPT FUNCTION #####

nanocompore_postprocess (
    input_tsv=snakemake.input.res_tsv,
    fasta=snakemake.input.fasta,
    outdir=os.path.dirname(snakemake.output[0]),
    log = snakemake.log[0],
    p_val_lim=snakemake.params.opt.get("p_val_lim", 0.01),
    quantile_lim=snakemake.params.opt.get("quantile_lim", 0.5),
    min_distance=snakemake.params.opt.get("min_distance", 9))
