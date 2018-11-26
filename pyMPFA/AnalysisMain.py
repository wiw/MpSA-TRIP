#!/usr/bin/env python
# encoding: utf-8

# from pprint import pprint as view
import os
import pickle
import re
import copy
import json
import subprocess
import gzip
import datetime
import logging
# from progressbar import ProgressBar, Bar
from glob import glob
from string import Template
# from pprint import pprint as view
from json2html import *
from collections import OrderedDict
# Configuration dictionary for lib 29-36 with bcread=3, bcmut_probability=0.8
# CONFIG = {
#     "experiment_dir": "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_3__bcmutProb_80/Lib_29-36",
#     "content": ["control_e", "control_n", "control_m", "expression", "normalization"],
#     "exception": {
#         "experiment_dir": "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_3__bcmutProb_80/Lib_29-36",
#         "content": ["mapping"]
#     },
#     "control": {
#         "wt-bc1": "TTCCAAGTGCAGGTTAGGCG",
#         "wt-bc2": "TGTGTACGGCTTGCTCTCAA",
#         "deltaC-bc3": "GAGCCCGGATCCACTCCAAG",
#         "deltaC-bc4": "TGTCACGTCAGCTAACCCAC"
#     },
#     "statistics_output": "/home/anton/backup/input/trip/RUN_2018-05-10/results/statistics/Lib_29-36_bcRead_3__bcmutProb_80_bmc_genomics",
#     "rscript": "/usr/bin/Rscript",
#     "output_control": "control.json",
#     "output_data": "data.json",
#     "output_rpl_count": "rpl_count.json",
#     "html_template": "/home/anton/data/TRIP/pyMPFA/report.html.tpl",
#     "pympfa_src": "/home/anton/data/TRIP/pyMPFA"
# }
# Config for lib 33-40 Run from 2018-10-19
# CONFIG = {
#     "experiment_dir": "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_3__bcmutProb_80/Lib_33-40",
#     "content": ["control_e", "control_n", "control_m", "expression", "normalization"],
#     "exception": {
#       "experiment_dir": "/home/anton/backup/input/trip/RUN_2018-05-10/results/article_bmc_genomics/Lib_33-40",
#       "content": ["mapping"]
#     },
#     "control": {
#         "wt-bc1": "TTCCAAGTGCAGGTTAGGCG",
#         "wt-bc2": "TGTGTACGGCTTGCTCTCAA",
#         "deltaC-bc3": "GAGCCCGGATCCACTCCAAG",
#         "deltaC-bc4": "TGTCACGTCAGCTAACCCAC"
#     },
#     "statistics_output": "/home/anton/backup/input/trip/RUN_2018-05-10/results/statistics/Lib_33-40_bmc_genomics",
#     "rscript": "/usr/bin/Rscript",
#     "output_control": "control.json",
#     "output_data": "data.json",
#     "output_rpl_count": "rpl_count.json",
#     "html_template": "/home/anton/data/TRIP/pyMPFA/report.html.tpl",
#     "pympfa_src": "/home/anton/data/TRIP/pyMPFA"
# }

# indexList = {"m1": "CAAGATAA", "m2": "GGACAACG"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_BC_Mut_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_BC_Mut_mapping"

# Configuration dictionary for lib 29-36 with bcread=3, bcmut_probability=0.8
CONFIG = {
    "experiment_dir": "/home/anton/backup/input/trip/RUN_2018-11-20/results/bcRead_3__bcMutProb_80__bcError_2/Lib_37-44",
    "content": ["control_e", "control_n", "control_m", "expression", "normalization"],
    "exception": {
        "experiment_dir": "/home/anton/backup/input/trip/RUN_2018-11-20/results/bcRead_2__bcMutProb_50__bcError_2/Lib_37-44",
        "content": ["mapping"]
    },
    "control": {
        "wt-bc1": "TTCCAAGTGCAGGTTAGGCG",
        "wt-bc2": "TGTGTACGGCTTGCTCTCAA",
        "deltaC-bc3": "GAGCCCGGATCCACTCCAAG",
        "deltaC-bc4": "TGTCACGTCAGCTAACCCAC"
    },
    "statistics_output": "/home/anton/backup/input/trip/RUN_2018-11-20/results/statistics/Lib_37-44__bcRead_2__bcMutProb_50__bcError_2",
    "rscript": "/usr/bin/Rscript",
    "output_control": "control.json",
    "output_data": "data.json",
    "output_rpl_count": "rpl_count.json",
    "html_template": "/home/anton/data/TRIP/pyMPFA/report.html.tpl",
    "pympfa_src": "/home/anton/data/TRIP/pyMPFA"
}

if not os.path.exists(CONFIG["statistics_output"]):
    os.makedirs(CONFIG["statistics_output"])
logging.basicConfig(filename=os.path.join(CONFIG["statistics_output"], "analysis.log"),
                    level=logging.INFO,
                    format='\n%(asctime)s - %(levelname)s\n%(message)s',
                    datefmt='%Y.%m.%d %H:%M:%S')
Logger = logging.getLogger(__name__)


def head(list_or_dict, start=0, n=10):
    start, n = int(start), int(int(start) + int(n))
    if type(list_or_dict) == dict:
        for i in dict(list_or_dict.items()[start:n]):
            print i, dict(list_or_dict.items()[start:n])[i]
    elif type(list_or_dict) == list:
        for i in list_or_dict[start:n]:
            print i


def dump_to_json(obj, output_file):
    try:
        with open(output_file, "w") as handle:
            json.dump(obj, handle)
    except:
        Logger.exception("Doesn't write to {}".format(output_file))


def SaveDictToPy(dictVar, filename):
    with open(filename + ".py", "wb") as handle:
        for k, v in dictVar.items():
            if type(v) == str:
                handle.write(str(k) + " = '" + str(v) + "'\n")
            else:
                handle.write(str(k) + " = " + str(v) + "\n")


def load_pickle(CONFIG):
    all_data = {}
    control_data = {}
    rpl_count_data = []

    def pickle_opener(i):
        with open(i, "rb") as handle:
            unserialized_data = pickle.load(handle)
        return unserialized_data

    def open_input(INPUT):
        if INPUT.endswith('gz'):
            return gzip.open(INPUT, 'rb')
        else:
            return open(INPUT, 'rb')

    def get_records(input_file):
        with open_input(input_file) as f:
            TotalSeqRecords = int(sum(1 for _ in f)) / 4
        return TotalSeqRecords

    def search_counting_replicates(item, path_to_rpl):
        out = []
        for fastq in glob(path_to_rpl):
            count = get_records(fastq)
            match = re.match(".*_([a-z])(?P<repl>[0-9])",
                             os.path.basename(os.path.dirname(fastq)))
            if match is not None:
                replicate_id = item + "_" + match.group("repl")
            else:
                replicate_id = item + "_Und"
            out.append((replicate_id, count))
        return out
    for item in CONFIG["content"]:
        regex0 = re.compile(".*" + item)
        replicate_dir = filter(
            regex0.match, os.listdir(CONFIG["experiment_dir"]))[0]
        path_to = os.path.join(CONFIG["experiment_dir"], replicate_dir, "Dump")
        path_to_rpl = os.path.join(
            CONFIG["experiment_dir"], replicate_dir, "**/index*.fastq")
        if re.search("control.*", item) is None:
            rpl_count_data.extend(
                search_counting_replicates(item, path_to_rpl))
        if item in ("control_m", "mapping"):
            what_is_pickle = "resultDict"
        else:
            what_is_pickle = "bcDict"
        regex1 = re.compile(".*" + what_is_pickle + ".*")
        pickles = filter(regex1.match, os.listdir(path_to))
        for f in pickles:
            filename, ext = os.path.splitext(f)
            filename = re.sub("-", "_", filename)
            fst, snd, repl, dct = filename.split("_")
            del(fst, snd, ext, dct)
            filename = item + "-" + repl[1:]
            f = os.path.join(path_to, f)
            if re.search("control", item) is not None:
                control_data[filename] = pickle_opener(f)
            else:
                all_data[filename] = pickle_opener(f)
    if CONFIG.get("exception") is not None and bool(CONFIG["exception"]):
        exception_all_data, exception_control_data, exception_rpl_count_data = load_pickle(
            CONFIG["exception"])
        all_data.update(exception_all_data)
        control_data.update(exception_control_data)
        rpl_count_data.extend(exception_rpl_count_data)
    return all_data, control_data, rpl_count_data


def get_control_count(CONFIG, control):
    compiled_control = {}

    def convert_and_count_map(mapping_data):
        mut, count, mutated_variants = mapping_data
        for mut_bc, mut_list in mutated_variants.items():
            for sets in mut_list:
                mutm, countm = sets
                count += countm
        return count
    for item in control:
        compiled_control.setdefault(item, copy.deepcopy(CONFIG["control"]))
        for alias, bc in CONFIG["control"].items():
            try:
                if re.search("control_m.*", item) is None:
                    counter = sum(
                        [count for mutated_bc, count in control[item][bc]])
                else:
                    counter = convert_and_count_map(control[item][bc])
            except KeyError:
                counter = 0
                Logger.exception(
                    "Barcode {} for alias {} is not found".format(bc, alias))
            compiled_control[item][alias] = counter
    return compiled_control


def align_map_norm_count_expr(data, CONFIG):

    def _get_count(doubled_sets_in_list):
        return sum([count for item, count in doubled_sets_in_list])

    def _get_sets(CONFIG):
        sets = [i for i in CONFIG["content"] if not re.search("control", i)]
        if CONFIG.get("exception") is not None and bool(CONFIG["exception"]):
            sets_add = [i for i in CONFIG["exception"]
                        ["content"] if not re.search("control", i)]
            sets.extend(sets_add)
        return sets
    sets = _get_sets(CONFIG)
    map_norm_1_2_data, expr_1_2_data = {}, {}
    fmt_map_data = {}
    for experiment in sets:
        get_replicates = [i for i in data.keys() if re.search(experiment, i)]
        # For mapping - compare replicates with each other and return common combinations in set (bc, seq)
        if experiment == "mapping":
            for repl in get_replicates:
                tmp = {bc: mut[0] for bc, mut in data[repl].items()}
                fmt_map_data[repl] = tmp
            if len(get_replicates) == 2:  # intersect two replicates
                first, second = get_replicates
                unique_mapping = [comb for comb in fmt_map_data[first].viewitems(
                ) if comb in fmt_map_data[second].viewitems()]
            elif len(get_replicates) == 3:  # intersect three replicates
                first, second, third = get_replicates
                unique_mapping = [comb for comb in fmt_map_data[first].viewitems(
                ) if comb in fmt_map_data[second].viewitems() and comb in fmt_map_data[third].viewitems()]
            elif len(get_replicates) == 4:  # intersect four replicates
                first, second, third, fourth = get_replicates
                unique_mapping = [comb for comb in fmt_map_data[first].viewitems(
                ) if comb in fmt_map_data[second].viewitems() and comb in fmt_map_data[third].viewitems() and comb in fmt_map_data[fourth].viewitems()]
            unique_mapping_data = {bc: mut for bc, mut in unique_mapping}
            map_norm_1_2_data[experiment] = unique_mapping_data
        # For normalization - compare replicates with each other and return common barcodes, but with its own value in each replica
        if experiment == "normalization":
            for repl in get_replicates:
                map_norm_1_2_data.setdefault(repl, {})
            if len(get_replicates) == 2:
                first, second = get_replicates
                for bc, value in data[first].viewitems():
                    if bc in data[second]:
                        map_norm_1_2_data[first][bc] = _get_count(value)
                        map_norm_1_2_data[second][bc] = _get_count(
                            data[second][bc])
            elif len(get_replicates) == 3:
                first, second, third = get_replicates
                for bc, value in data[first].viewitems():
                    if bc in data[second] and bc in data[third]:
                        map_norm_1_2_data[first][bc] = _get_count(value)
                        map_norm_1_2_data[second][bc] = _get_count(
                            data[second][bc])
                        map_norm_1_2_data[third][bc] = _get_count(
                            data[third][bc])
        # For expression - don't compare replicates, instead, we count the reads for each barcode considering mutated barcodes
        if experiment == "expression":
            for repl in get_replicates:
                expr_1_2_data.setdefault(repl, {})
                for bc, value in data[repl].items():
                    expr_1_2_data[repl][bc] = _get_count(value)
    return map_norm_1_2_data, expr_1_2_data


def align_map_vs_norm_replicates(map_norm_1_2_data):
    mapped_ratio_data = {}
    _tmp_norm_keys, union_norm_data = [], {}
    normalization_sets = [
        i for i in map_norm_1_2_data.keys() if re.search("normalization", i)]
    # Get union keys from both normalization replicates
    for repl in normalization_sets:
        _tmp_norm_keys.extend(map_norm_1_2_data[repl].keys())
    union_norm_keys = set().union(_tmp_norm_keys)
    # Making a dictionary for each barcode with the values of both normalizing replicates
    for bc in union_norm_keys:
        union_norm_data.setdefault(bc, {})
        for repl in normalization_sets:
            union_norm_data[bc].setdefault(repl, map_norm_1_2_data[repl][bc])
    # Intersect union normalization dictionary with mapping data
    if map_norm_1_2_data.get("mapping") is not None and bool(map_norm_1_2_data["mapping"]):
        for bc, mut in map_norm_1_2_data["mapping"].items():
            if bc in union_norm_data:
                mapped_ratio_data.setdefault(bc, {"mutation": mut})
                mapped_ratio_data[bc].update(union_norm_data[bc])
    return mapped_ratio_data


def align_mapped_ratio_vs_expr_replicates(mapped_ratio_data, expr_1_2_data):
    mapped_norm_expression_data = {}
    # Intersect "mapping-normalization" data with expression data
    for bc, value in mapped_ratio_data.items():
        mapped_norm_expression_data.setdefault(bc, value)
        for repl in expr_1_2_data:
            # If barcode exists in expression data then associated with him expression data, but ...
            if bc in expr_1_2_data[repl]:
                mapped_norm_expression_data[bc].update(
                    {repl: expr_1_2_data[repl][bc]})
            else:
                # if not, then barcode assigned ZERO.
                mapped_norm_expression_data[bc].update({repl: 0})
    return mapped_norm_expression_data


def make_report(REPORT, CONFIG, label=False):
    def count_effective_reads(data, item_name):
        eff_count = []
        if re.search("mapping", item_name) is None:
            for key, value in data[item_name].items():
                eff_count.append(sum([_count for _string, _count in value]))
        else:
            for key, value in data[item_name].items():
                eff_count.append(sum([x for x in value if type(x) == int]))
                mutation_dict_existance = [x for x in value if type(x) == dict]
                if len(mutation_dict_existance) == 1:
                    mutation_dict = mutation_dict_existance[0]
                    for _key, _value in mutation_dict.items():
                        eff_count.append(
                            sum([_count for _string, _count in _value]))
        return sum(eff_count)
    # path to template and set filename
    name, ext_one, ext_two = os.path.basename(
        CONFIG["html_template"]).split(".")
    current = datetime.datetime.now().strftime("%Y_%m_%d")
    if label is not False:
        path_to_html_report = os.path.join(
            CONFIG["statistics_output"], label, current + "_" + label + "_" + name + "." + ext_one)
        write_to_report = dict(library=" ".join(
            [label, os.path.basename(CONFIG["experiment_dir"])]))
    else:
        path_to_html_report = os.path.join(
            CONFIG["statistics_output"], current + "_" + name + "." + ext_one)
        write_to_report = dict(
            library=os.path.basename(CONFIG["experiment_dir"]))
    if not os.path.exists(os.path.dirname(path_to_html_report)):
        os.makedirs(os.path.dirname(path_to_html_report))
    # Make tables for report
    write_to_report["current"] = re.sub("_", ".", current)
    if os.path.exists(CONFIG["html_template"]) and os.path.isfile(CONFIG["html_template"]):
        for element in REPORT:
            if element in ["data", "map_norm_1_2_data", "expr_1_2_data"]:
                _tmp = {}
                for experiment in REPORT[element]:
                    _tmp.setdefault(experiment, len(
                        REPORT[element][experiment]))
                source_for_json = OrderedDict(sorted(_tmp.items()))
            if element in ["output_rpl_count", "output_control_data"]:
                source_for_json = REPORT[element]
            if element in ["mapped_ratio_data", "mapped_norm_expression_data"]:
                source_for_json = len(REPORT[element])
            write_to_report[element] = json2html.convert(
                json=source_for_json, table_attributes="class=\"table table-bordered table-hover\"")
        _t_tmp = {}
        for item in REPORT["data"]:
            _t_tmp.setdefault(
                item, count_effective_reads(REPORT["data"], item))
        _source_for_json = OrderedDict(sorted(_t_tmp.items()))
        write_to_report["effective_reads_count"] = json2html.convert(
            json=_source_for_json, table_attributes="class=\"table table-bordered table-hover\"")
        with open(CONFIG["html_template"], "rb") as handle:
            _for_substitute = Template(
                handle.read()).safe_substitute(write_to_report)
            with open(path_to_html_report, "wb") as out_file:
                out_file.write(_for_substitute)


def biological_sense(CONFIG, use_method="mpsa", label=False):
    pathToScript = os.path.join(
        CONFIG["pympfa_src"], "StatAnalysisMain_fork.R")
    option = [use_method]
    for item in ["statistics_output", "output_control", "output_data", "output_rpl_count"]:
        ins = label if label is not False else ""
        if item == "statistics_output":
            wd_path = os.path.join(CONFIG["statistics_output"], ins)
            if os.path.exists(wd_path):
                option.append(wd_path)
            continue
        path_to = os.path.join(CONFIG["statistics_output"], ins, CONFIG[item])
        if os.path.exists(path_to):
            option.append(path_to)
    cmd = [CONFIG["rscript"], "--verbose", pathToScript] + option
    Logger.info(cmd)
    subprocess.call(cmd)


def simple_dict_writer(CONFIG, dictionary, filename):
    filepath = os.path.join(CONFIG["statistics_output"], filename + ".csv")
    with open(filepath, "wb") as handle:
        fieldnames = ['key', 'value']
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for k, v in dictionary.items():
            writer.writerow({
                'key': k,
                'value': v})


def main(CONFIG):
    # prepare data to analysis
    data, control, output_rpl_count = load_pickle(CONFIG)
    output_control_data = get_control_count(CONFIG, control)
    # aligning and reads counting
    map_norm_1_2_data, expr_1_2_data = align_map_norm_count_expr(data, CONFIG)
    mapped_ratio_data = align_map_vs_norm_replicates(map_norm_1_2_data)
    mapped_norm_expression_data = align_mapped_ratio_vs_expr_replicates(
        mapped_ratio_data, expr_1_2_data)
    # prepare to report & make report
    REPORT = {
        "output_rpl_count": OrderedDict(sorted({k: v for k, v in output_rpl_count}.items())),
        "output_control_data": OrderedDict(sorted(output_control_data.items())),
        "data": OrderedDict(sorted(data.items())),
        "map_norm_1_2_data": OrderedDict(sorted(map_norm_1_2_data.items())),
        "expr_1_2_data": OrderedDict(sorted(expr_1_2_data.items())),
        "mapped_ratio_data": mapped_ratio_data,
        "mapped_norm_expression_data": mapped_norm_expression_data
    }
    OUTPUT = {
        "output_data": mapped_norm_expression_data,
        "output_control": output_control_data,
        "output_rpl_count": output_rpl_count
    }
    make_report(REPORT, CONFIG)
    # write output data to json
    try:
        for _out in OUTPUT:
            dump_to_json(OUTPUT[_out], os.path.join(
                CONFIG["statistics_output"], CONFIG[_out]))
            Logger.info("{} write succesful! ^_^".format(_out))
    except:
        Logger.exception("{} write failed! >_<".format(_out))
    # Run biological analysis in R implementation
    # biological_sense(CONFIG)      # disabled


if __name__ == '__main__':
    main(CONFIG)
