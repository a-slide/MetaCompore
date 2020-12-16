# -*- coding: utf-8 -*-

##### Imports #####

import gzip
import logging

##### DEFINE SCRIPT FUNCTION #####
def nanopolish_eventalign_gather(input_tsv_list, output_tsv, log):

    logging.basicConfig(filename=log, filemode="w", level=logging.DEBUG)

    try:
        # Concatenate and manage gzip compression
        write_open_fun, write_open_mode = (gzip.open, "wt") if output_tsv.endswith(".gz") else (open, "w")

        with write_open_fun(output_tsv, write_open_mode) as output_fp:
            first = True
            for input_tsv in input_tsv_list:
                logging.debug(f'Reading file {input_tsv}\n')
                read_open_fun, read_open_mode = (gzip.open, "rt") if input_tsv.endswith(".gz") else (open, "r")
                with read_open_fun(input_tsv, read_open_mode) as input_fp:
                    if not first:
                        # flush header if not first file
                        _ = input_fp.readline()
                    else:
                        first=False
                    output_fp.write(input_fp.read())
            logging.debug('All files done\n')
    except:
        logging.exception('Error while running nanopolish_eventalign_gather')
        raise

##### RUN SCRIPT FUNCTION #####

nanopolish_eventalign_gather(
    input_tsv_list = snakemake.input.tsv_list,
    output_tsv = snakemake.output.tsv,
    log = snakemake.log[0])
