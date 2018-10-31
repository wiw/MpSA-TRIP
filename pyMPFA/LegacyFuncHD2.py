#!/usr/bin/env python
# #!C:/Python27/python.exe
# encoding: utf-8
# This script make pwm matrix from cleared txt files who contained only barcode sequences


import os
import pandas as pd
import regex
from weblogolib import *
from collections import Counter, OrderedDict
import PairedEndFunc as pend

CONFIG = {"input": "/home/anton/backup/input/trip/RUN_2018-06-07/Episomal_expression_results/bc_sequences",
          "output_sub_folder": "results",
          }


def reverse_idx_promotor(promotor):
        return pend.reverseComplement(promotor.upper())


def load_files(CONFIG):
    path = CONFIG["input"]
    files = [x for x in os.listdir(path) if x.endswith(".txt")]
    loaded_data = {}
    for item in files:
        abspath_item = os.path.join(path, item)
        with open(abspath_item, "r") as handle:
            for line in handle:
                match = regex.search("(N.*)", line)
                if match is None:
                    formatted_line = pend.reverseComplement(line.rstrip().upper())
                    loaded_data.setdefault(item, []).append(formatted_line)
    return loaded_data, files


def convert_count2portion(count_data):
    try:
        portion_data = {}
        if isinstance(count_data, dict):
            all_count = [i for i in count_data.values()]
            if all([isinstance(value, int) for value in all_count]):
                all_count_sum = sum(all_count)
                for item, value in count_data.items():
                    portion = float(value) / all_count_sum
                    portion_data.setdefault(item, portion)
            return portion_data
    except:
        raise TypeError


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


def convert_data_to_fastq(CONFIG, loaded_data):
    path = os.path.join(CONFIG["input"], CONFIG["output_sub_folder"])
    for file_name, bcs in loaded_data.items():
        file_core, ext = os.path.splitext(file_name)
        file_path = os.path.join(path, file_core + ".fa")
        with open(file_path, "w") as fout:
            c = 1
            for bc in bcs:
                fout.write(">BARCODE_{}\n{}\n".format(c, bc))
                c += 1


def make_positional_matrix(loaded_data):
    pwm_data = {}
    for files, bcs in loaded_data.items():
        pwm_data.setdefault(files, get_positional_matrix(bcs))
    return pwm_data


def w_pwm(pwm_data, file_path):
    with open(file_path, "w") as handle:
        pd_data = pd.DataFrame(pwm_data.values(), columns=['A', 'C', 'G', 'T'])
        pd_data = pd_data.T
        pd_data.to_csv(handle, sep="\t", )


def main_write(CONFIG, pwm_data):
    folder = os.path.join(CONFIG["input"], CONFIG["output_sub_folder"])
    if not os.path.exists(folder):
        os.makedirs(folder)
    for files, pwm in pwm_data.items():
        file_path = os.path.join(folder, "Position_weight_matrix_" + files)
        w_pwm(pwm, file_path)


def main(CONFIG):
    loaded_data, files = load_files(CONFIG)
    convert_data_to_fastq(CONFIG, loaded_data)
    pwm_data = make_positional_matrix(loaded_data)
    main_write(CONFIG, pwm_data)

if __name__ == '__main__':
    main(CONFIG)
