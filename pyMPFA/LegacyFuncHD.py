#!/usr/bin/env python
# #!C:/Python27/python.exe
# encoding: utf-8


import os
import copy
import gzip
import regex
import pickle
import platform
import random
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import Counter, OrderedDict
import SupportFunc as supp
import PairedEndFunc as pend

CONFIG = {"Windows": "F:/Dropbox/1_S1_L001_R1_001.fastq.gz",
          "Linux": "/home/anton/backup/input/trip/RUN_2018-06-07/1_S1_L001_R1_001.fastq.gz",
          "idx_illumina": ["AGGCAGCA", "AGCTTTCT"],
          "idx_promotor": ["ttgag", "gatgg", "agctc", "tgtgt", "agtca",
                           "tggcc", "acgta", "ctagt", "ctgct", "caatt",
                           "tcaaa", "actgc", "tcgct", "ccgag"],
          "const": {1: ('cgccagggttttcccagtcacaagggccggccacaactc', 4),
                    2: ('ctc', 1)},
          "regexp_template": "^___idx_illumina___const__1___(?P<barcode>.*)([ATGC]{4})___idx_promotor___const__2___.*",
          }


def open_input(INPUT):
    try:
        if INPUT.endswith('gz'):
            return gzip.open(INPUT, 'rb')
        else:
            return open(INPUT, 'rb')
    except:
        print("Opener error!")


def reverse_idx_promotor(CONFIG):
    temp_idx = [pend.reverseComplement(idx_pr.upper()) for idx_pr in CONFIG["idx_promotor"]]
    return temp_idx


def save_pickle(obj, name, path):
    locationObj = os.path.join(path, name)
    filename = locationObj + ".pickle"
    with open(filename, 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return filename


def load_pickle(name, path):
    locationObj = os.path.join(path, name)
    with open(locationObj + '.pickle', 'rb') as handle:
        unserialized_data = pickle.load(handle)
        return unserialized_data


def convert_count2portion(count_data):
    try:
        portion_data = {}
        if isinstance(count_data, dict):
            all_count = [i for i in count_data.values()]
            if all([isinstance(value, int) for value in all_count]):
                all_count_sum = sum(all_count)
                for item, value in count_data.items():
                    portion = round(float(value) / all_count_sum, 2)
                    portion_data.setdefault(item, portion)
            return portion_data
    except:
        raise TypeError


def select_source_file(CONFIG):
    OS = platform.system()
    return CONFIG[OS]


def make_samples_for_testing(string, length, count):
    try:
        if isinstance(string, str) and isinstance(length + count, int):
            output = []
            for _ in range(0, count):
                i = 0
                word = ""
                while i < length:
                    word += random.choice(string)
                    i += 1
                output.append(word)
            return output
        else:
            raise TypeError
            print("Use correct input data. Arguments sequence: sample string, length of type 'int' of each sequence, count of all sequences in type 'int'")
    except:
        raise Exception

# To obtain the frequency of occurrence of a letter in a certain position, source data are used in this format:
# ['atgc', 'atgt', 'ggta', 'ggat', 'actg', 'cgct', 'gtca', 'attc', 'gcat', 'tatg']
# First, the rows are expanded into a vertical "matrix". Namely - for each letter
# of each item in the list - this letter is highlighted in a separate list in the dictionary with a key
# corresponding to the position of the letter. Then the function collections.Counter consider the frequency of occurrence
# of each letter for each key in the dictionary. Further, the frequency of occurrence is converted into fractions
# occurrences. The resulting dictionary is written to the csv file.


def get_positional_matrix(dataset):
    try:
        p_matrix = {}
        check = True if isinstance(dataset, list) else False
        check = True if isinstance(dataset[0], str) and all(dataset) else False
        strings_length = {len(s) for s in dataset}
        check = True if len(strings_length) == 1 else False
        if check:
            strings_length = int(list(strings_length)[0])
            for letter_pos in range(0, strings_length):
                for string in dataset:
                    p_matrix.setdefault(letter_pos, []).append(string[letter_pos])
                if bool(p_matrix.get(letter_pos)):
                    temp_var = dict(Counter(p_matrix[letter_pos]))
                    temp_var = convert_count2portion(temp_var)
                    p_matrix[letter_pos] = OrderedDict(sorted(temp_var.iteritems(), key=lambda x: x[0]))
            return p_matrix
    except:
        raise TypeError

# Formation of a regular expression from CONFIG at the stage of the search for barcodes to determine the Illumina index, promoter index and search for barcodes.


def make_regexp(CONFIG):
    def wrapper_appeal(name, match):
        return "(?P<{}>{})".format(name, match)
    template = CONFIG["regexp_template"].split("___")
    regexp = ""
    for element in template:
        if len(element.split("__")) == 1:
            # single element
            if CONFIG.get(element, False):
                # element of dict need to parse
                if isinstance(CONFIG[element], list):
                    match = "|".join(CONFIG[element]).upper()
                    regexp += wrapper_appeal(element, match)
            else:
                # simple text element
                regexp += element
        elif len(element.split("__")) == 2:
            # inherit element
            # correct template
            master_key, slave_key = element.split("__")
            slave_key = int(slave_key)
            # next parse dict and fill regexp
            if CONFIG.get(master_key, False):
                if CONFIG[master_key].get(slave_key, False):
                    if isinstance(CONFIG[master_key][slave_key], tuple):
                        sequence, error_value = CONFIG[master_key][slave_key]
                        regexp += "({}){{s<={}}}".format(sequence.upper(), str(error_value))
                else:
                    # Doesn't match slave key, skipping element
                    continue
            else:
                # Doesn't match master key, skipping element
                continue
        else:
            # bad template, skipping
            continue
    return regex.compile(regexp)


# Task I.
# Reading the source file, parsing lines in accordance with the regular expression "regex_template"
# and the recording of bar codes in memory with a mark of the corresponding Illumina index and promoter index

def fastq_parse(CONFIG):
    source_data = select_source_file(CONFIG)
    regexp = make_regexp(CONFIG)
    bc_data = {}
    if os.path.isfile(source_data):
        reads_count = supp.get_sequence_count(source_data)
        print("Fastq file: {}; reads count: {}".format(source_data, reads_count))
        with open_input(source_data) as fastq:
            read_number = 0
            for title, seq, qual in FastqGeneralIterator(fastq):
                read_number += 1
                match = regexp.match(seq)
                if match is not None:
                    bc_data.setdefault(match.group("idx_illumina"), {})
                    bc_data[match.group("idx_illumina")].setdefault(match.group("idx_promotor"), []).append(match.group("barcode"))
                if read_number % 10 ** 5 == 0:
                    print("Processed {} reads from {}".format(read_number, reads_count))
    return bc_data
# From the dictionary with a specific set of indexes ("idx_illumina" -> "idx_promotor") collect and count unique barcodes.
# For each set, calculate how many barcodes how long is in the sample


def analysis_of_data(bc_data):
    analyzed_data = {}
    if isinstance(bc_data, dict) and bc_data:
        for idx_illumina, idx_promotor_data in bc_data.items():
            analyzed_data.setdefault(idx_illumina, {})
            for idx_promotor, core_data in idx_promotor_data.items():
                unique_barcodes = dict(Counter(core_data))
                length_of_unq_barcodes = [len(unq_bc) for unq_bc, count in unique_barcodes.items()]
                counter_of_length_unq_barcodes = dict(Counter(length_of_unq_barcodes))
                counter_of_length_unq_barcodes = OrderedDict(sorted(counter_of_length_unq_barcodes.items(), key=lambda x: x[0]))
                analyzed_data[idx_illumina].setdefault(idx_promotor, counter_of_length_unq_barcodes)
        return analyzed_data


def align_analyzed_data(analyzed_data):
    max_value = 23
    for idx_illumina in analyzed_data:
        for idx_promotor in analyzed_data[idx_illumina]:
            for pos_key in range(max_value):
                analyzed_data[idx_illumina][idx_promotor].setdefault(pos_key, 0)
                analyzed_data[idx_illumina][idx_promotor] = OrderedDict(sorted(analyzed_data[idx_illumina][idx_promotor].items(), key=lambda x: x[0]))
    return analyzed_data
# Task II.
# Based on the previous data, select barcodes only 18 nk long.
# Conduct their frequency analysis


def eighteen_choice(bc_data):
    eighteen_data = {}
    if isinstance(bc_data, dict) and bc_data:
        for idx_illumina, idx_promotor_data in bc_data.items():
            eighteen_data.setdefault(idx_illumina, {})
            for idx_promotor, core_data in idx_promotor_data.items():
                eighteen_barcodes = [bc for bc in core_data if len(bc) == 18]
                if eighteen_barcodes:
                    eighteen_matrix = get_positional_matrix(eighteen_barcodes)
                    eighteen_data[idx_illumina].setdefault(idx_promotor, eighteen_matrix)
        return eighteen_data
# Write received data to file
# Task III
# Write to file barcode sequences with length 18bp


def w_eighteen_sequences(bc_data, path):
    file_path = os.path.join(path, 'Sequences_of_eighteen_bcs.xlsx')
    writer = pd.ExcelWriter(file_path, engine='xlsxwriter')
    for idx_illumina, idx_promotor_data in bc_data.items():
        start_col_number, delimiter = 0, 3
        for idx_promotor, bc_sequences in idx_promotor_data.items():
            print(idx_illumina, idx_promotor, start_col_number)
            eighteen_barcodes = [bc for bc in bc_sequences if len(bc) == 18]
            pd_data = pd.DataFrame(eighteen_barcodes)
            pd_data.to_excel(writer, sheet_name=idx_illumina, index_label=idx_promotor, startcol=start_col_number)
            start_col_number += len(pd_data.columns) + delimiter
    writer.save()


def main_write(CONFIG, bc_data, analyzed_data, eighteen_data):
    file_path = os.path.dirname(select_source_file(CONFIG))
    w_analysis(analyzed_data, file_path)
    w_eighteen(format_eighteen_data(eighteen_data), file_path)
    w_eighteen_sequences(bc_data, file_path)
# Preparation and recording of data on the distribution of lengths of barcodes in the found data


def w_analysis(analyzed_data, path):
    file_path = os.path.join(path, 'Length_distribution.xlsx')
    writer = pd.ExcelWriter(file_path, engine='xlsxwriter')
    for idx_illumina, idx_promotor_data in analyzed_data.items():
        start_col_number, delimiter = 0, 1
        for idx_promotor, core_data in idx_promotor_data.items():
            pd_data = pd.DataFrame(core_data.items(), columns=["length", "count"])
            pd_data.to_excel(writer, sheet_name=idx_illumina, startcol=start_col_number, index_label=idx_promotor)
            start_col_number += len(pd_data.columns) + delimiter
    writer.save()
# Preparation and recording of data on frequency analysis of 18-nk bar codes


def w_eighteen(eighteen_data, path):
    file_path = os.path.join(path, 'Position_eighteen_bc_distribution.xlsx')
    writer = pd.ExcelWriter(file_path, engine='xlsxwriter')
    for idx_illumina, idx_promotor_data in eighteen_data.items():
        start_row_number, delimiter = 0, 3
        for idx_promotor, letter_matrix in idx_promotor_data.items():
            pd_data = pd.DataFrame(letter_matrix.values(), columns=['A', 'C', 'G', 'T'])
            pd_data = pd_data.T
            pd_data.to_excel(writer, sheet_name=idx_illumina, index_label=idx_promotor, startrow=start_row_number)
            start_row_number += len(pd_data.index) + delimiter
    writer.save()
# Formatting data of frequency analysis of 18-nk bar codes. Remove undefined nucleotide.


def format_eighteen_data(eighteen_data):
    formatted_data = {}
    for idx_illumina, idx_promotor_data in eighteen_data.items():
        formatted_data.setdefault(idx_illumina, {})
        for idx_promotor, letter_matrix in idx_promotor_data.items():
            formatted_data[idx_illumina].setdefault(idx_promotor, OrderedDict())
            for position, ordered_values in letter_matrix.items():
                data_without_N = [x for x in ordered_values.values() if x != 0.0]
                formatted_data[idx_illumina][idx_promotor].update({position: data_without_N})
    return formatted_data


def main(CONFIG):
    CONFIG_MOD = copy.deepcopy(CONFIG)
    CONFIG_MOD["idx_promotor"] = reverse_idx_promotor(CONFIG)
    bc_data = fastq_parse(CONFIG_MOD)
    analyzed_data = align_analyzed_data(analysis_of_data(bc_data))
    eighteen_data = eighteen_choice(bc_data)
    main_write(CONFIG, bc_data, analyzed_data, eighteen_data)


if __name__ == '__main__':
    main(CONFIG)
    print("End of time ^_^")
