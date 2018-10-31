#!/usr/bin/env python
# encoding: utf-8
import os
import csv
import regex
import datetime
import json
import logging
import logging.config
import pickle
import math
import gzip
import copy
import pandas as pd
import LegacyFuncHD as lf
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from matplotlib import pyplot as plt
from collections import OrderedDict
"""
**** SUPPORT FUNCTIONS
"""


# def GetTotalSeqRecords(input_file):
#     import subprocess
#     salt = "zcat -f | " if input_file.endswith("gz") else ""
#     cmd = '{}wc -l {}'.format(salt, input_file)
#     result = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
#     out, err = result.communicate()
#     count, path = out.split()
#     TotalSeqRecords = int(count) / 4
#     return TotalSeqRecords


def GetTotalSeqRecords(input_file):
    '''
    This function count of number strings in fastq file and return sequences number
    '''
    def open_input(INPUT):
        if INPUT.endswith('gz'):
            return gzip.open(INPUT, 'rb')
        else:
            return open(INPUT, 'rb')
    with open_input(input_file) as f:
        TotalSeqRecords = int(sum(1 for _ in f)) / 4
    return TotalSeqRecords


def head(list_or_dict, start=0, n=10):
    start, n = int(start), int(int(start) + int(n))
    if type(list_or_dict) == dict:
        for i in dict(list_or_dict.items()[start:n]):
            print i, dict(list_or_dict.items()[start:n])[i]
    elif type(list_or_dict) == list:
        for i in list_or_dict[start:n]:
            print i


def simpleWrite(list_or_dict, path, name, delim="\t"):
    filePath = os.path.join(path, name)
    with open(filePath, "wb") as handle:
        if type(list_or_dict) == dict:
            header = ["id", "value"]
            writer = csv.DictWriter(handle, fieldnames=header, delimiter=delim)
            writer.writeheader()
            for i in list_or_dict:
                writer.writerow({
                    "id": i,
                    "value": list_or_dict[i]})
        elif type(list_or_dict) == list:
            header = ["id"]
            writer = csv.DictWriter(handle, fieldnames=header, delimiter=delim)
            writer.writeheader()
            for i in list_or_dict:
                writer.writerow({
                    "id": i})


def jsonWrite(var_dict, outpath, outfilename):
    path_to_out = os.path.join(outpath, outfilename)
    with open(path_to_out, 'wb') as handle:
        json.dump(var_dict, handle)


def EstimateCalculationTime(obj):
    if type(obj) == dict or type(obj) == list:
        l = float(len(obj))
        velocity = 5*10**5
        sec = (float(l**2)/2)/velocity
        now = datetime.datetime.now().strftime("%H:%M")
        if sec <= 60:
            t = str(float(sec)) + " sec. Current time: " + now
        elif sec <= 3600 and sec > 60:
            t = str(float(sec)/(60)) + " min. Current time: " + now
        elif sec <= 24*60**2 and sec > 3600:
            t = str(float(sec)/(60**2)) + "h."
        else:
            t = str(float(sec)/(24*60**2)) + "d. Current time: " + now
        return t
    return "I don't estimate calculation time for this object. Sorry..."


def makeStatFromBowtieAlign(bwtAlignerDict):
    LogInfo('Some statistics from bowtie aligner')
    expr = regex.compile('[ ]*(?P<count>\d.*) ((\((?P<pct>[0-9\.\%].*)\) .*aligned (?P<times>.*time.*|))|reads.*)')
    for key in bwtAlignerDict:
        LogInfo("Bowtie align report for {}".format(key))
        for line in bwtAlignerDict[key]:
            m = expr.match(line)
            if m is not None:
                if m.group('times') is not None:
                    LogInfo("     Aligned {}: {} ({})".format(m.group('times'), m.group('count'), m.group('pct')))
                else:
                    LogInfo("     {} input reads".format(m.group('count')))


def make_histogramm_plot(pkl):
    def pickle_opener(i):
        with open(i, "rb") as handle:
            unserialized_data = pickle.load(handle)
        return unserialized_data

    def _get_count(doubled_sets_in_list):
        return sum([count for item, count in doubled_sets_in_list])
    data = pickle_opener(pkl)
    length_data = len(data)
    _bins = math.log(length_data, 2) + 1
    dumpdir = os.path.dirname(pkl)
    filename, ext = os.path.splitext(os.path.basename(pkl))
    hset = [_get_count(value) for key, value in data.items()]
    plt.figure(figsize=(10, 10))
    plt.title("Histogramm plot for " + filename)
    n, bins, patches = plt.hist(hset, _bins)
    plt.savefig(os.path.join(dumpdir, filename + ".png"), fmt='png')
    return length_data


def setup_logging(default_path='logging.json', default_level=logging.INFO, env_key='LOG_CFG'):
    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
    if os.path.exists(path):
        with open(path, 'rt') as f:
            config = json.load(f)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=default_level)


def LogInfo(msg):
    logger = logging.getLogger(__name__)
    logger.info(msg)


def LogErr(msg):
    logger = logging.getLogger(__name__)
    logger.error(msg)

##################################
# Human report generator functions
##################################


def open_input(INPUT):
    try:
        if INPUT.endswith('gz'):
            return gzip.open(INPUT, 'rb')
        else:
            return open(INPUT, 'rb')
    except:
        LogErr("Opener error!")

# Function generate unordered html list (<ul>) from simple python dict with fixed css class 'long_string'.
# And replace signs more&less in text string


def dict2html_list(any_dict):
    try:
        if isinstance(any_dict, dict):
            def replace_brackets(s):
                if type(s) == str:
                    pattern = OrderedDict([("<", "&lt;"), (">", "&gt;")])
                    for old, new in pattern.items():
                        s = s.replace(old, new)
                    return s
                return s

            def unordered(fn):
                def wrapper(any_dict):
                    return "<ul class='long_string'>" + fn(any_dict) + "</ul>"
                return wrapper

            @unordered
            def make_dict_to_html_list(any_dict):
                html_list = ["<li><strong>{}</strong> - {}</li>".format(k, replace_brackets(v)) for k, v in any_dict.items()]
                html_list = "".join(html_list)
                return html_list
            return make_dict_to_html_list(any_dict)
    except:
        raise TypeError


def get_source_file(options):
    if options.get("r1") is False or options.get("r2") is False:
        return options['input_file']
    keys = ["r1", "r2"]
    foo = {x: options[x] for x in keys}
    return dict2html_list(foo)


def extract_experiment_and_lib(options):
    # Regular expression splits text string by example: "trip-6_2_mapping" -> "trip_6_2", "mapping"
    foo = regex.split('([a-z\-_0-9]*)_([a-z]*)', options['experiment_label'])
    foo = [x for x in foo if len(x) != 0]
    return foo


def make_all_index_stat(options):
    def make_regexp(idx, idx_error):
        regexp = "^(?P<index>{}){{s<={}}}.*".format(idx,
                                                    str(idx_error))
        return regex.compile(regexp)
    output_stat, reads_count = {}, 0
    with open_input(options["forward"]) as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            reads_count += 1
            for idx_error in range(3):
                output_stat.setdefault(idx_error, {})
                for idx in options["all_indexes"]:
                    regexp = make_regexp(idx, idx_error)
                    match = regexp.match(seq)
                    if match is not None:
                        output_stat[idx_error][idx] = output_stat[idx_error].setdefault(idx, 0) + 1
    return output_stat, reads_count


def format_all_index_stat_to_html(output_stat, reads_count, to_html=True):
    output_data = pd.DataFrame()
    pd.set_option('display.max_colwidth', -1)
    for idx_error, stat in output_stat.items():
        pd_data = pd.DataFrame(stat.items(), columns=['index', str(idx_error)]).sort_values(by=['index'])
        if len(output_data) == 0:
            output_data = copy.deepcopy(pd_data)
        else:
            output_data = output_data.merge(pd_data, left_on='index', right_on='index', how='outer')
    output_data.loc['total'] = pd.Series(output_data.sum(), index=[str(x) for x in output_stat.keys()])
    output_data.ix['total', 'index'] = reads_count
    if to_html:
        return output_data.to_html().replace("\n", "")
    return output_data


def count_non_barcoded_reads(options):
    import PairedEndFunc as pend

    def reverse_idx_promotor(promotor):
        return pend.reverseComplement(promotor.upper())
    non_barcoded_data, c = {}, 0
    DEFAULT = {"const1": "CGCCAGGGTTTTCCCAGTCACAAGGGCCGGCCACAACTCGTCGAC",
               "const2": "CTCGATCTCTAGACCCTCCG",
               "idx_illumina": ["AGGCAGCA", "AGCTTTCT"],
               "idx_promotor": ["ttgag", "gatgg", "agctc", "tgtgt", "agtca",
                                "tggcc", "acgta", "ctagt", "ctgct", "caatt",
                                "tcaaa", "actgc", "tcgct", "ccgag"],
               }
    with open_input(options['forward']) as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            for idx_illumina in DEFAULT['idx_illumina']:
                non_barcoded_data.setdefault(idx_illumina, {})
                for idx_promotor in DEFAULT['idx_promotor']:
                    regex_string = regex.compile("{}({}){{s<=4}}{}{{s<=1}}{}{{s<=2}}".format(idx_illumina,
                                                                                             DEFAULT['const1'],
                                                                                             reverse_idx_promotor(idx_promotor),
                                                                                             DEFAULT['const2']))
                    match = regex_string.match(seq)
                    if match is not None:
                        non_barcoded_data[idx_illumina][idx_promotor] = non_barcoded_data[idx_illumina].setdefault(idx_promotor, 0) + 1
            c += 1
            if c % 10**4 == 0:
                print("Current loaded {} reads".format(c))
    return non_barcoded_data


def collect_report_data(options):
    report_data = {}
    # Header section
    report_data['experiment'], report_data['lib'] = extract_experiment_and_lib(options)
    report_data['date'] = datetime.datetime.now().strftime("%Y-%m-%d")
    # Options section
    report_data['options'] = dict2html_list(options)
    # De-multiplexing section
    # Generate html table with statistical data about splitting by ALL indexes in view of variable error count in indexes ${all_index_stat}
    output_stat, reads_count = make_all_index_stat(options)
    report_data["all_index_stat"] = format_all_index_stat_to_html(output_stat, reads_count)
    # Generate html table with stat_data by using at the moment indexes.
    tmp = format_all_index_stat_to_html(output_stat, reads_count, False)
    report_data['current_index_stat'] = tmp[tmp['index'].isin(options['indexList'].values())][["index", "0"]]
    # Reads count before/after comm by promotor index. Collect data in 'colb.main_paired' Return results as html table
    # #################### smth code
    # Bowtie statistics in table view $alignment_stat_table +  filtration statistics
    # #################### smth code
    # Finally statistics in $output_stat_table
    # #################### smth code


def parse_collection_data(collection_of_output_data):
    # Collect all the found four-letter sequences and make a frequency analysis of the occurrence of letters ATGC.
    # Necessary intervention in the process of creating a dictionary bcDict. Return html table
    # #################### smth code
    # Further use only files from 'paired_indexes'. Need to load path to this files to report
    # #################### smth code
    # Count length of genome from the fwd and rev reads and use the shortest
    # #################### smth code
    pass


def main_report(collection_of_output_data, options):
    report_data = collect_report_data(options)
    # Safe substitute in template by the simple options. When this is used built-in dict - "options"
    pass
