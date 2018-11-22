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
from collections import OrderedDict, Counter
from string import Template
from functools import reduce  # forward compatibility for Python 3
import operator

"""
**** SUPPORT FUNCTIONS
"""


# def get_sequence_count(input_file):
#     import subprocess
#     salt = "zcat -f | " if input_file.endswith("gz") else ""
#     cmd = '{}wc -l {}'.format(salt, input_file)
#     result = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
#     out, err = result.communicate()
#     count, path = out.split()
#     TotalSeqRecords = int(count) / 4
#     return TotalSeqRecords


def get_sequence_count(input_file):
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


def get_sequence_count_2(input_file):
    '''
    This function count of number strings in fastq file and return sequences number
    '''
    count = 0
    if not input_file.endswith('gz'):
        with open(input_file) as f:
            for line in f.xreadlines():
                count += 1
        return count / 4
    return 0


def head(list_or_dict, start=0, n=10):
    start, n = int(start), int(int(start) + int(n))
    if type(list_or_dict) == dict:
        for i in dict(list_or_dict.items()[start:n]):
            print i, dict(list_or_dict.items()[start:n])[i]
    elif type(list_or_dict) == list:
        for i in list_or_dict[start:n]:
            print i


def simple_write(list_or_dict, path, name, delim="\t"):
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


def json_write(var_dict, outpath, outfilename):
    path_to_out = os.path.join(outpath, outfilename)
    with open(path_to_out, 'wb') as handle:
        json.dump(var_dict, handle)


def estimate_calculation_time(obj):
    if type(obj) == dict or type(obj) == list:
        l = float(len(obj))
        velocity = 5 * 10 ** 5
        sec = (float(l ** 2) / 2) / velocity
        now = datetime.datetime.now().strftime("%H:%M")
        if sec <= 60:
            t = str(float(sec)) + " sec. Current time: " + now
        elif sec <= 3600 and sec > 60:
            t = str(float(sec) / (60)) + " min. Current time: " + now
        elif sec <= 24 * 60 ** 2 and sec > 3600:
            t = str(float(sec) / (60 ** 2)) + "h."
        else:
            t = str(float(sec) / (24 * 60 ** 2)) + "d. Current time: " + now
        return t
    return "I don't estimate calculation time for this object. Sorry..."


def make_stat_from_bowtie(bwtAlignerDict):
    log_info('Some statistics from bowtie aligner')
    expr = regex.compile(
        '[ ]*(?P<count>\\d.*) ((\\((?P<pct>[0-9\\.\\%].*)\\) .*aligned (?P<times>.*time.*|))|reads.*)')
    for key in bwtAlignerDict:
        log_info("Bowtie align report for {}".format(key))
        for line in bwtAlignerDict[key]:
            m = expr.match(line)
            if m is not None:
                if m.group('times') is not None:
                    log_info("     Aligned {}: {} ({})".format(
                        m.group('times'), m.group('count'), m.group('pct')))
                else:
                    log_info("     {} input reads".format(m.group('count')))


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


def log_info(msg):
    logger = logging.getLogger(__name__)
    logger.info(msg)


def log_error(msg):
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
        log_error("Opener error!")


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
                html_list = ["<li><strong>{}</strong> - {}</li>".format(
                    k, replace_brackets(v)) for k, v in any_dict.items()]
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
    # Regular expression splits text string by example: "trip-6_2_mapping" ->
    # "trip_6_2", "mapping"
    foo = regex.split('([a-z\\-_0-9]*)_([a-z]*)', options['experiment_label'])
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
                        output_stat[idx_error][idx] = output_stat[idx_error].setdefault(
                            idx, 0) + 1
            if reads_count % 10 ** 5 == 0:
                print("Prepaired {} reads".format(reads_count))
    return output_stat, reads_count


def format_all_index_stat(output_stat, reads_count):
    output_data = pd.DataFrame()
    pd.set_option('display.max_colwidth', -1)
    for idx_error, stat in output_stat.items():
        pd_data = pd.DataFrame(stat.items(), columns=[
            'index', str(idx_error)]).sort_values(by=['index'])
        if len(output_data) == 0:
            output_data = copy.deepcopy(pd_data)
        else:
            output_data = output_data.merge(
                pd_data, left_on='index', right_on='index', how='outer')
    output_data.loc['total'] = pd.Series(
        output_data.sum(), index=[str(x) for x in output_stat.keys()])
    output_data.ix['total', 'index'] = reads_count
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
                    regex_string = regex.compile("{}({}){{s<=4}}{}{{s<=1}}{}{{s<=2}}"
                                                 .format(idx_illumina,
                                                         DEFAULT['const1'],
                                                         reverse_idx_promotor(
                                                             idx_promotor),
                                                         DEFAULT['const2']))
                    match = regex_string.match(seq)
                    if match is not None:
                        non_barcoded_data[idx_illumina][idx_promotor] = non_barcoded_data[idx_illumina].setdefault(
                            idx_promotor, 0) + 1
            c += 1
            if c % 10 ** 4 == 0:
                print("Current loaded {} reads".format(c))
    return non_barcoded_data


def collect_report_data(options):
    report_data = {}
    # Header section
    report_data['experiment'], report_data['lib'] = extract_experiment_and_lib(
        options)
    report_data['date'] = datetime.datetime.now().strftime("%Y-%m-%d")
    # Options section
    report_data['options'] = dict2html_list(options)
    # De-multiplexing section
    # Generate html table with statistical data about splitting by ALL indexes
    # in view of variable error count in indexes ${all_index_stat}
    output_stat, reads_count = make_all_index_stat(options)
    report_data["all_index_stat"] = format_all_index_stat(
        output_stat, reads_count)
    # Generate html table with stat_data by using at the moment indexes.
    report_data['current_index_stat'] = report_data["all_index_stat"][report_data[
        "all_index_stat"]['index'].isin(options['indexList'].values())][["index", "0"]]
    for index_stat in ["all_index_stat", "current_index_stat"]:
        report_data[index_stat] = report_data[index_stat].to_html(
            classes=["table", "table-hover", "table-sm"], justify="center").replace("\n", "")
    # Count length of genome from the fwd and rev reads and use the shortest
    report_data['exp_fwd_genome_len'] = str(options['expected_min_genome_len'])
    encoded_report_data = {k: v.encode('ascii', 'ignore')
                           for k, v in report_data.items()}
    return encoded_report_data


def parse_collection_data(collection_of_output_data):
    structure = {'spacer_4freq': ['main_paired_stat', 'four_letters_seq_collection', 'fwd'],
                 'current_reads_count': ['paired_indexes'],
                 'comm_table': ['main_paired_stat', 'comm_stat'],
                 'bowtie_table': ['main_paired_stat', 'bwt_aligner_stat'],
                 'output_stat_table': {'genuine bc': ['bcDict'],
                                       'genuine bc & genome combinations': ['seqDict'],
                                       'genuine "bc & genome" combinations after filter': ['resultDict']},
                 }

    add_to_report = {}

    # Collect all the found four-letter sequences and make a frequency analysis of the occurrence of letters ATGC.
    # Return html table

    def stat_four_seq_data(idata):
        output= {}
        for promotor_index, four_seq_data in idata.items():
            output.setdefault(promotor_index, dict(Counter(four_seq_data)))
        return output

    def format_four_seq_data(sdata):
        def wrapper_head(promotor_index, pandas_df):
            return "<div><p>{}</p>{}</div>".format(promotor_index, pandas_df.replace("\n", ""))

        output = ""
        for promotor_index, stat_four_data in sdata.items():
            tmp_df = pd.DataFrame(stat_four_data.items(), columns=[
                'four_seq', 'count']).sort_values(by=['count'], ascending=False)
            # Magic number - 15, first 15 elements from statistics
            tmp_df = tmp_df.head(15)
            tmp_df_html = tmp_df.to_html(
                index=False, classes=["table", "table-hover", "table-sm"], justify="center")
            output += wrapper_head(promotor_index, tmp_df_html)
        return output

    def get_by_map(data_dict, map_list):
        return reduce(operator.getitem, map_list, data_dict)

    # Get reads count from current file
    def get_current_reads_data(cdata):
        output = {}
        for path in cdata.values():
            output.setdefault(path, str(get_sequence_count_2(path)))
        return output

    # Reads count before/after comm by promotor index. Collect data in 'colb.main_paired' Return results as html table
    # #################### smth code  comm_table  comm_difference
    def main_comm_stat(comm_tables):
        def format_comm_stat(comm_item, fst_column_name):
            tmp_df = pd.DataFrame(comm_item.items()).T
            tmp_df.columns = ['fwd', 'comm', 'rev']
            tmp_df = tmp_df[1:]
            tmp_df['pmi'] = str(fst_column_name)
            return tmp_df

        def merge_comm_stat(df_comm_list):
            tmp_df = pd.concat(df_comm_list)
            tmp_df = tmp_df[['pmi', 'fwd', 'rev', 'comm']]
            return tmp_df
        comm_list = []
        for promotor_index, tables in comm_tables.items():
            comm_list.append(format_comm_stat(tables, promotor_index))
        merged_comm = merge_comm_stat(comm_list)
        merged_comm_html = merged_comm.to_html(
            index=False, classes=["table", "table-hover", "table-sm"], justify="center")
        return merged_comm_html

    # Bowtie statistics in table view $alignment_stat_table +  filtration statistics
    # TODO: promotor index split by "_" -> bowtie stat parse "aligned
    # concordantly 0 times" & "aligned concordantly exactly 1 time" & "aligned
    # concordantly >1 times" -> to dict -> to dataFrame -> concat

    def main_bowtie_stat(bowtie_table):
        def parse_bowtie_log(log, regexp_compiled):
            output = []
            for line in log:
                match = regexp_compiled.match(line)
                if match is not None:
                    count, pct = match.group('count'), match.group('pct')
                    output.append([int(count), float(pct)])
            return output

        def collect_bowtie_log(bowtie_table):
            regexp = regex.compile(
                '^([ ].*?)(?P<count>[0-9].*) \\((?P<pct>[\\.0-9].*)%\\) .* concordantly (?P<header>[a-z0-9 >].*)$')
            collect_result = {}
            promotor_counter = 0
            for promotor_index, bowtie_output in bowtie_table.items():
                parse_result = parse_bowtie_log(bowtie_output, regexp)
                _sum = [sum(x) for x in zip(*parse_result)]
                promotor_index = promotor_index.split("_")[0]
                output_parse_result = [promotor_index] + parse_result
                output_parse_result.append(_sum)
                collect_result.setdefault(
                    promotor_counter, output_parse_result)
                promotor_counter += 1
            return collect_result

        def format_bowtie_log(collect_result):
            pd_table = pd.DataFrame(collect_result.values())
            pd_table.columns = ['promotor index', 'aligned conc. 0 times',
                                'aligned conc. exactly 1 time', 'aligned conc. >1 times', 'summary']
            bowtie_html = pd_table.to_html(
                index=False, classes=["table", "table-hover", "table-sm"], justify="center")
            return bowtie_html
        collect_result = collect_bowtie_log(bowtie_table)
        bowtie_html = format_bowtie_log(collect_result)
        return bowtie_html

    # Finally statistics in $output_stat_table
    # #################### smth code   output_stat_table
    def main_output_stat(collection_of_output_data, replicate):
        set_of_columns = ['genuine bc', 'genuine bc & genome combinations',
                          'genuine "bc & genome" combinations after filter']
        output_table = pd.DataFrame(columns=set_of_columns)
        for synonim, name_of_dict in structure['output_stat_table'].items():
            received_data = get_by_map(collection_of_output_data, [
                                       replicate] + name_of_dict)
            for promotor_index, actual_data in received_data.items():
                output_table.at[promotor_index, synonim] = len(actual_data)
        output_table['promotor_index'] = output_table.index
        output_table = output_table[['promotor_index'] + set_of_columns]
        output_to_html = output_table.to_html(
            index=False, classes=["table", "table-hover", "table-sm"], justify="center")
        return output_to_html

    replicate_count = 1
    for replicate in collection_of_output_data:
        # spacer_4freq
        four_seq_data_counting = stat_four_seq_data(get_by_map(
            collection_of_output_data, [replicate] + structure['spacer_4freq']))
        four_seq_data_counting_html = format_four_seq_data(
            four_seq_data_counting)
        add_to_report.setdefault(
            "{}_{}".format("spacer_4freq", replicate_count), four_seq_data_counting_html)
        # current_reads_data
        current_reads_data = get_current_reads_data(get_by_map(
            collection_of_output_data, [replicate] + structure['current_reads_count']))
        current_reads_data_html = dict2html_list(current_reads_data)
        add_to_report.setdefault("{}_{}".format(
            "current_reads_data", replicate_count), current_reads_data_html)
        # comm
        comm_table = main_comm_stat(get_by_map(collection_of_output_data, [
            replicate] + structure['comm_table']))
        add_to_report.setdefault("{}_{}".format(
            "comm_table", replicate_count), comm_table)
        # bowtie_table
        bowtie_table = main_bowtie_stat(get_by_map(collection_of_output_data, [
                                        replicate] + structure['bowtie_table']))
        add_to_report.setdefault("{}_{}".format(
            "bowtie_table", replicate_count), bowtie_table)
        # output_stat_table
        output_table = main_output_stat(collection_of_output_data, replicate)
        add_to_report.setdefault("{}_{}".format(
            "output_table", replicate_count), output_table)
        replicate_count += 1
    encoded_add_to_report_data = {k: v.encode('ascii', 'ignore')
                           for k, v in add_to_report.items()}
    return encoded_add_to_report_data


def write_report_to_template(report_data, options):
    # Safe substitute in template by the simple options. When this is used
    # built-in dict - "options"
    all_subst, all_subst_data = [report_data, options], {}
    for element in all_subst:
        all_subst_data.update(element)
    current = datetime.datetime.now().strftime("%Y_%m_%d")
    html_template_name = "mapping_report.html.tpl"
    html_output_name, tpl_ext = os.path.splitext(html_template_name)
    html_output_name = "{}_{}_{}".format(
        current, report_data['experiment'], html_output_name)
    html_output_path = os.path.join(options["workdir"], html_output_name)
    with open(html_template_name, "r") as tpl:
        substitution = Template(tpl.read()).safe_substitute(all_subst_data)
        with open(html_output_path, "w") as out:
            out.write(substitution)


def main_report(collection_of_output_data, options):
    report_data = collect_report_data(options)
    report_data.update(parse_collection_data(collection_of_output_data))
    write_report_to_template(report_data, options)
