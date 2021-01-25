# -*- coding: utf-8 -*-

##### Imports #####

import os
import gzip
import logging
from collections import Counter, defaultdict

##### DEFINE SCRIPT FUNCTION #####

def min_ref_coverage (reads_index_list, ref_list, min_cov, log):

    logging.basicConfig(filename=log, filemode="w", level=logging.DEBUG)

    try:
        logging.debug("Parse index files")
        d = defaultdict(Counter)
        all_ref = set()
        for fn in reads_index_list:
            sample_id = os.path.basename(fn)
            with gzip.open(fn, "rt") as fp:
                header = next(fp)
                for l in fp:
                    ref_id = l.split("\t")[1]
                    d[sample_id][ref_id]+=1
                    all_ref.add(ref_id)
        logging.debug("Total references: {}".format(len(all_ref)))

        logging.debug("Count valid references per sample")
        valid_ref_count = Counter()
        for sample_id, ref_d in d.items():
            for ref_id, i in ref_d.items():
                if i >= min_cov:
                    valid_ref_count[ref_id]+=1

        logging.debug("Select valid references for all samples")
        valid_ref_list = []
        for ref_id, i in valid_ref_count.items():
            if i == len(reads_index_list):
                valid_ref_list.append(ref_id)
        logging.debug("Valid references: {}".format(len(valid_ref_list)))

        logging.debug("Writing valid references to file")
        with open(ref_list, "w") as fp:
            for ref in valid_ref_list:
                fp.write("{}\n".format(ref))

    except:
        logging.exception('Error while running min_ref_coverage')
        raise

##### RUN SCRIPT FUNCTION #####

min_ref_coverage(
    reads_index_list = snakemake.input.reads_index_list,
    ref_list = snakemake.output.ref_list,
    min_cov = snakemake.params.opt.get("min_cov", 0),
    log = snakemake.log[0])
