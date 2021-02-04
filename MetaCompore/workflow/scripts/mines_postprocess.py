# -*- coding: utf-8 -*-

##### Imports #####

import logging
import pandas as pd
import numpy as np
import datetime

##### DEFINE SCRIPT FUNCTION #####
def mines_postprocess (res_bed, res_tsv, log, min_cov):

    logging.basicConfig(filename=log, filemode="w", level=logging.INFO, format='%(message)s')
    logging.info("timestamp: {}".format(str(datetime.datetime.now())))
    for i, j in locals().items():
        logging.info("\t{}: {}\n".format(i,j))

    logging.info("Loading data")
    df = pd.read_csv(res_bed, sep="\t", names=["refid", "pos", "end", "kmer", "unique key", "strand", "fraction modified", "coverage"], usecols=[0,1,3,6,7])

    logging.info("Filtering out low coverage position")
    df = df[df["coverage"]>=min_cov]

    logging.info("Write output file")
    df.to_csv(res_tsv, index=False, sep="\t")

##### RUN SCRIPT FUNCTION #####
mines_postprocess (
    res_bed=snakemake.input.res_bed,
    res_tsv=snakemake.output.res_tsv,
    log = snakemake.log[0],
    min_cov=snakemake.params.opt.get("min_cov", 30))
