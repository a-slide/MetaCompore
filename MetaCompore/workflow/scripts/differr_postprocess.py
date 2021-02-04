# -*- coding: utf-8 -*-

##### Imports #####

import logging
import pandas as pd
import numpy as np
import datetime

##### DEFINE SCRIPT FUNCTION #####
def differr_postprocess (res_bed, res_tsv, log):

    logging.basicConfig(filename=log, filemode="w", level=logging.INFO, format='%(message)s')
    logging.info("timestamp: {}".format(str(datetime.datetime.now())))
    for i, j in locals().items():
        logging.info("\t{}: {}\n".format(i,j))

    min_pval = np.nextafter(float(0), float(1))
    max_pval = 1

    logging.info("Loading data")
    df = pd.read_csv(res_bed, sep="\t", usecols=[0,1,4,6,7,8,9],
                     names=[
                        "refid",
                        "pos",
                        "end",
                        "name",
                        "score",
                        "strand",
                        "odds_ratio",
                        "G_stat",
                        "-log10_pval",
                        "-log10_FDR",
                        "G_stat_control",
                        "-log10_pval_control",
                        "G_stat_test",
                        "-log10_pval_test"])

    logging.info("Cleanup data and calculate pvalue")
    df = df.replace(np.inf, np.nan)
    df = df.dropna(subset=["score", "-log10_pval"])
    df["pvalue"]=np.power(10, -df["-log10_pval"]).clip(min_pval, max_pval)
    df["adj_pvalue"]=np.power(10, -df["-log10_FDR"]).clip(min_pval, max_pval)
    df=df[["refid", "pos", "pvalue", "adj_pvalue", "odds_ratio", "G_stat"]]

    logging.info("Write output file")
    df.to_csv(res_tsv, index=False, sep="\t")

##### RUN SCRIPT FUNCTION #####
differr_postprocess (
    res_bed=snakemake.input.res_bed,
    res_tsv=snakemake.output.res_tsv,
    log=snakemake.log[0])
