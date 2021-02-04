# -*- coding: utf-8 -*-

##### Imports #####

from collections import OrderedDict
import logging
from pyBioTools import Fasta
from pyfaidx import Faidx
import datetime

##### DEFINE SCRIPT FUNCTION #####

def get_transcriptome(fa_input, fa_output, fai_output, log):

    logging.basicConfig(filename=log, filemode="w", level=logging.INFO, format='%(message)s')
    logging.info("timestamp: {}".format(str(datetime.datetime.now())))
    for i, j in locals().items():
        logging.info("\t{}: {}\n".format(i,j))

    try:
        # Parse fasta file uncompress and simplify transcript ids
        logging.info("Read input transcriptome fasta file")
        with open(fa_output, "w") as fa_out:
            for rec in Fasta.Reader(fa_input):
                fa_out.write(">{}\n{}\n".format(rec.short_name, rec.seq))

        logging.info("Index fasta file")
        with Faidx(fa_output) as fa_out:
            fa_out.build_index()

    except:
        logging.exception('Error while running get_transcriptome')
        raise

##### RUN SCRIPT FUNCTION #####

get_transcriptome(
    fa_input = str(snakemake.input.fasta), # str required for FTP or HTTP sources
    fa_output = snakemake.output.fasta,
    fai_output = snakemake.output.fai,
    log = snakemake.log[0])
