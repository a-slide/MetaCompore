# -*- coding: utf-8 -*-

##### Imports #####

import logging
import pandas as pd
import numpy as np
import datetime

##### DEFINE SCRIPT FUNCTION #####
def eligos2_postprocess (
    input_tsv,
    output_tsv,
    log,
    min_oddR,
    min_esb,
    min_cov,
    max_adj_pval,
    discard_homopolymers,
    ref_base):

    logging.basicConfig(filename=log, filemode="w", level=logging.INFO, format='%(message)s')
    logging.info("timestamp: {}".format(str(datetime.datetime.now())))
    for i, j in locals().items():
        logging.info("\t{}: {}\n".format(i,j))

    min_pval = np.nextafter(float(0), float(1))
    max_pval = 1

    logging.info("Loading data")
    df = pd.read_csv(input_tsv, sep='\t')
    df = df.dropna(subset=["adjPval"])
    df["adjPval"] = df["adjPval"].clip(min_pval, max_pval)

    logging.info("Cleanup data and calculate pvalue")
    if min_oddR is not None:
        df = df[df["oddR"] >= min_oddR]
    if min_esb is not None:
        df = df[df["ESB_test"] >= min_esb]
    if min_cov is not None:
        df = df[df["total_reads"] >= min_cov]
    if max_adj_pval is not None:
        df = df[df["adjPval"] <= max_adj_pval]
    if discard_homopolymers is True:
        df = df[df["homo_seq"] == "--"]
    if ref_base in ["A","T","C","G"]:
        df = df[df["ref"] == ref_base]

    df = df[['chrom', 'start_loc', 'pval', 'adjPval', 'oddR']]
    df = df.rename(columns={'chrom':"refid", 'start_loc':"pos", 'pval':"pvalue", 'adjPval':"adj_pvalue", 'oddR':"odds_ratio"})

    logging.info("Write output file")
    df.to_csv(output_tsv, index=False, sep="\t")

##### RUN SCRIPT FUNCTION #####
eligos2_postprocess (
    input_tsv=snakemake.input.res_tsv,
    output_tsv=snakemake.output.res_tsv,
    log=snakemake.log[0],
    min_oddR=snakemake.params.opt.get("min_oddR", 1.2),
    min_esb=snakemake.params.opt.get("min_esb", 0),
    min_cov=snakemake.params.opt.get("min_cov", 30),
    max_adj_pval=snakemake.params.opt.get("max_adj_pval", 0.01),
    discard_homopolymers=snakemake.params.opt.get("discard_homopolymers", True),
    ref_base=snakemake.params.opt.get("ref_base", "A"))
