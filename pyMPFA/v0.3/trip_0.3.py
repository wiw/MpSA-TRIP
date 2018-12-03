#!/usr/bin/env python
# encoding: utf-8
"""
    trip_0.3.py
    programm software version 0.3
    Search for barcodes and associated mutations.
    Report on the work done and statistics output.\n
    author: Anton V. Ivankin
    e-mail: anton.ivankin@gmail.com
    source: https://github.com/wiw/MpSA-TRIP/
"""

import os
import pickle
import subprocess
import argparse
import json
import logging
import logging.config
import gzip
import regex
import pdb
import string
import Levenshtein
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

Logger = logging.getLogger(__name__)
CONFIG_FILE_PATH = "config.json"


def load_main_config(config_path):
    try:
        if not os.path.exists(config_path):
            config_path = CONFIG_FILE_PATH
        with open(config_path, "rt") as handle:
            config = json.load(handle)
            return config
    except:
        Logger.exception("Don't load config file from {}".format(config_path))


def load_logging_config(path):
    if os.path.exists(path):
        with open(path, 'rt') as f:
            config = json.load(f)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(filename=None, level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y.%m.%d %H:%M:%S')


def parse_arguments():
    try:
        p = argparse.ArgumentParser(
            description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
        p.add_argument("--config", "-c", help="Path to config file")
        args = p.parse_args()
        return args
    except:
        Logger.exception("Don't parse arguments!")


def open_input(INPUT):
    try:
        if INPUT.endswith('gz'):
            return gzip.open(INPUT, 'rb')
        else:
            return open(INPUT, 'rb')
    except:
        Logger.exception("Opener error!")


def load_pickle(data):
    try:
        with open_input(data) as handle:
            unserialized_data = pickle.load(handle)
            return unserialized_data
    except:
        Logger.exception("Don't load your data {}".format(data))


def select_file(config, experiment):
    if config['input_file'].keys() == ["single"]:
        return [config['input_file']['single']['path']]
    elif config['input_file'].keys() == ["multiple"]:
        return config['input_file']['muliple'][experiment]


def reverseComplement(x):
    try:
        complement = string.maketrans('ACGTNSRYMKWHBVD', 'TGCANSRYMKWHBVD')
        if type(x) == unicode:
            complement = complement.decode("latin-1")
        revx = x.translate(complement)[::-1]
        return revx
    except:
        Logger.exception("Don't reverse string. Return initial value")
        return x


def regexp_maker(option, config, mod_read_structure=None):
    try:
        regexp = "^"
        experiment = option["experiment"]
        read_structure = config["content"][experiment]["read_structure"]
        if mod_read_structure:
            read_structure = mod_read_structure
    except:
        Logger.exception("Don't explain 'read_structure'")
    try:
        for key in read_structure:
            try:
                if type(key) is int:
                    regexp += "([ATGC]{{{}}})".format(str(key))
                else:
                    if key == "index":
                        index_error = str(config["core"]["error"].get("index_error", "0"))
                        index = option["index_seq"].upper()
                        regexp += "({}){{s<={}}}".format(index, index_error)
                    if regex.search("constant", key) is not None:
                        const_seq, const_err = config["content"][experiment][key]
                        const_seq, const_err = const_seq.upper(), str(const_err)
                        regexp += "({}){{s<={}}}".format(const_seq, const_err)
                    if key == "sub_index":
                        sub_index_str = [k.upper() for k in config["content"][experiment][key]]
                        if config["content"][experiment]["barcode_reversed"]:
                            sub_index_str = [reverseComplement(x) for x in sub_index_str]
                        sub_index_str = "|".join(sub_index_str)
                        sub_index_error = str(config["core"]["error"].get("sub_index_error", "0"))
                        regexp += "(?P<sub_index>{}){{s<={}}}".format(sub_index_str, sub_index_error)
                    if key == "barcode":
                        barcode_length = ",".join(
                            [str(int(round(config["core"]["length"]["barcode_length"] * x, 0))) for x in [0.9, 1.1]])
                        regexp += "(?P<barcode>[ATGC]{{{}}})".format(barcode_length)
            except:
                Logger.exception("Error with prepare key {]".format(key))
            try:
                if key == "sequence":
                    if config["project"] == "TRIP":
                        key_position_number = [i for i, x in enumerate(
                            config["content"][experiment]["read_structure"]) if x == key][0]
                        if key_position_number == len(config["content"][experiment]["read_structure"]) - 1:
                            regexp += "(?P<sequence_1>.*)"
                            return regexp
                        else:
                            mod_read_structure = config["content"][experiment]["read_structure"][key_position_number + 1:]
                            after_seq = regexp_maker(index, experiment, config, mod_read_structure)
                            regexp += "(((?P<sequence_0>.*){}.*)|(?P<sequence_1>.*))".format(after_seq)
                            return regexp
                    else:
                        sequence_length = str(config["core"]["length"].get("mutation_length", "8"))
                        regexp += "(?P<sequence>.{{{}}})".format(sequence_length)
            except:
                Logger.exception("Error with make regexp for 'sequence' key")
        return regexp
    except:
        Logger.exception("Undefined error in 'regexp_maker'")


def check_sub_index(dubbio_idx, config_sub_index, config_sub_index_error):
    try:
        chkResult = []
        for i in config_sub_index:
            if Levenshtein.distance(str(dubbio_idx), str(i)) <= int(config_sub_index_error):
                chkResult.append(str(i))
        if len(chkResult) == 1:
            return chkResult[0]
    except:
        return False


stat_dict = {
    "counter": {
        "all_lines": 0,
        "matched_barcodes": 0,
    }
}


def collect_grains(option, config, regexp, stat_dict):
    MAIN_DICT = {
        option["experiment"]: {
            option["index_name" + label]: {}
        }
    }
    with open_input(option["input"]) as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            match = regexp.match(seq)
            stat_dict["counter"]["all_lines"] += 1
            if match is not None:
                mgroup = {k: v for k, v in match.groupdict().items() if v is not None}


for k in mgroup:
    if regex.search("sequence", k) is None:
        mgroup[k] = reverseComplement(
            mgroup[k]) if config["content"][experiment]["barcode_reversed"] else mgroup[k]
    if k == "sub_index":
        mgroup[k] = check_sub_index(mgroup[k], config["content"][experiment]
                                    ["sub_index"], config["core"]["error"]["sub_index_error"])


# def main(args):
try:
    if os.path.exists(args.config):
        config = load_main_config(args.config)
        load_logging_config(config['logging'])
except:
    Logger.exception("There is no config on the specified path!")
try:
    for experiment in config['content']:
        INPUT = select_file(config, experiment)
        for f in INPUT:
            if f.endswith("pickle"):
                THIS_VAR_NAME_NEED_TO_BE_CHANGE = load_pickle(f)
                Logger.info("Data '{}' has been loaded from pickle file".format(f))
                continue
            for idx_name, idx_seq in config["content"][experiment]["index"].items():
                cur_option = {
                    "input": f,
                    "experiment": experiment,
                    "index_name": idx_name,
                    "index_seq": idx_seq
                }
                regexp = regex.compile(regexp_maker(cur_option, config))
            # MpSA: collect barcodes; combinations of barcode-mutation; check probability bc-mut; write result
            # TRIP: collect barcodes; combinations of: barcode-name_read, name_read-genome, name_read-genome_mapping and barcode-genome_mapping; check probability bc-gm; write result
except:
    Logger.exception("Undefined IOError!")


if __name__ == "__main__":
    try:
        args = parse_arguments()
        main(args)
    except:
        Logger.exception("Undefined error in <main>!")
