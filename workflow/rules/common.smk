# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
import sys
from collections import *
import shutil

# Third party lib
import pandas as pd
from snakemake.logging import logger
from snakemake.io import load_configfile
from snakemake.utils import validate

##### Input files validation functions #####

def config_load_validate(configfile, schema):
    config = load_configfile(configfile)
    validate(config, schema=schema)
    logger.debug(config)
    return config

def samples_load_validate(samplefile, schema):
    samples_df = pd.read_csv (samplefile, comment="#", skip_blank_lines=True, sep="\t")
    validate(samples_df, schema=schema)
    samples_df.set_index("sample_id", drop=True, inplace=True)
    logger.debug(samples_df)
    return samples_df

##### Getter functions #####

def get_fast5 (wildcards):
    return samples_df.loc[wildcards.sample, "fast5_dir"]

def get_threads (config, rule_name, default=1):
    try:
        return config[rule_name]["threads"]
    except (KeyError, TypeError):
        return default

def get_opt (config, rule_name, default=""):
    try:
        return config[rule_name]["opt"]
    except KeyError:
        return default

def get_mem (config, rule_name, default=1000):
    try:
        return config[rule_name]["mem"]
    except KeyError:
        return default

##### Helper functions #####

def mkdir (fn, exist_ok=False):
    """ Create directory recursivelly. Raise IO error if path exist or if error at creation """
    try:
        os.makedirs (fn, exist_ok=exist_ok)
    except:
        raise NanocomporeError ("Error creating output folder `{}`".format(fn))

def access_file (fn, **kwargs):
    """ Check if the file is readable """
    return os.path.isfile (fn) and os.access (fn, os.R_OK)
