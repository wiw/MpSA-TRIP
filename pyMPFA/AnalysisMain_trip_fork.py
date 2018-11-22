#!/usr/bin/env python
# encoding: utf-8

import AnalysisMain as alm
import os
import pickle
import re
import copy
import json
import subprocess
import gzip
import datetime
import logging
from glob import glob
from string import Template
from pprint import pprint as view
from json2html import *
from collections import OrderedDict

CONFIG = {
    "experiment_dir": "/home/anton/backup/input/trip/RUN_2018-06-07/results/trip_6_2",
    "content": ["expression", "normalization"],
    "exception": {
        "experiment_dir": "/home/anton/backup/input/trip/RUN_2018-06-07/results/trip_6_2__prob_80",
        "content": ["mapping"],
    },
    "statistics_output": "/home/anton/backup/input/trip/RUN_2018-06-07/results/statistics/trip_6_2__prob_80",
    "rscript": "/usr/bin/Rscript",
    "output_control": "control.json",
    "output_data": "data.json",
    "output_rpl_count": "rpl_count.json",
    "html_template": "/home/anton/data/TRIP/pyMPFA/report.html.tpl",
    "pympfa_src": "/home/anton/data/TRIP/pyMPFA",
}

if not os.path.exists(CONFIG["statistics_output"]):
    os.makedirs(CONFIG["statistics_output"])
logging.basicConfig(filename=os.path.join(CONFIG["statistics_output"], "analysis.log"),
                    level=logging.INFO,
                    format='\n%(asctime)s - %(levelname)s\n%(message)s',
                    datefmt='%Y.%m.%d %H:%M:%S')
Logger = logging.getLogger(__name__)


def reformatting_data(data):
    def get_promotor_indexes(data):
        indexes = []
        for key, value in data.items():
            indexes.extend(value.keys())
        indexes = set(indexes)
        return indexes
    pmi = get_promotor_indexes(data)
    experiment_set = data.keys()
    new_data = {}
    for pmi_item in pmi:
        new_data.setdefault(pmi_item, {})
        for experiment_set_item in experiment_set:
            new_data[pmi_item].setdefault(
                experiment_set_item, data[experiment_set_item][pmi_item])
    return new_data, pmi


def make_filename(CONFIG, item):
    now = datetime.datetime.now().strftime("%Y_%m_%d")
    fn = os.path.join(CONFIG["statistics_output"],
                      "{}_report_{}.html".format(now, item))
    return fn


def main(CONFIG):
    # prepare data to analysis
    data, control, output_rpl_count = alm.load_pickle(CONFIG)
    del control
    output_rpl_count = list(set(output_rpl_count))
    data, indexes = reformatting_data(data)
    # aligning and reads counting
    for pmi_item in indexes:
        map_norm_1_2_data, expr_1_2_data = alm.align_map_norm_count_expr(
            data[pmi_item], CONFIG)
        mapped_ratio_data = alm.align_map_vs_norm_replicates(map_norm_1_2_data)
        mapped_norm_expression_data = alm.align_mapped_ratio_vs_expr_replicates(
            mapped_ratio_data, expr_1_2_data)
        REPORT = {
            "output_rpl_count": OrderedDict(sorted({k: v for k, v in output_rpl_count}.items())),
            "data": OrderedDict(sorted(data[pmi_item].items())),
            "map_norm_1_2_data": OrderedDict(sorted(map_norm_1_2_data.items())),
            "expr_1_2_data": OrderedDict(sorted(expr_1_2_data.items())),
            "mapped_ratio_data": mapped_ratio_data,
            "mapped_norm_expression_data": mapped_norm_expression_data
        }
        OUTPUT = {
            "output_data": mapped_norm_expression_data,
            "output_rpl_count": output_rpl_count
        }
        alm.make_report(REPORT, CONFIG, label=pmi_item)
        try:
            for _out in OUTPUT:
                json_output = os.path.join(
                    CONFIG["statistics_output"], pmi_item, CONFIG[_out])
                if not os.path.exists(os.path.dirname(json_output)):
                    os.makedirs(os.path.dirname(json_output))
                alm.dump_to_json(OUTPUT[_out], json_output)
                Logger.info(
                    "{} {} write succesful! ^_^".format(pmi_item, _out))
        except:
            Logger.exception("{} {} write failed! >_<".format(pmi_item, _out))
        # Run biological analysis in R implementation
        alm.biological_sense(CONFIG, use_method="trip", label=pmi_item)


if __name__ == '__main__':
    main(CONFIG)
