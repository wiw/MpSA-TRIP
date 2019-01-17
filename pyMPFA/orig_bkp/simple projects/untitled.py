#!/usr/bin/env python
# encoding: utf-8
import os
import regex
import csv


wd = "/home/anton/backup/input/trip/RUN_2018-06-07/philion/"
filename = os.path.join(wd, "GSE22069_kc167_r6.tsv")
output = os.path.join(wd, "GSE22069_kc167_r6_reformat.tsv")
regexp = regex.compile("(?P<seqname>.*):(?P<start>.*)\\.{2}(?P<end>.*)")


def format_kccell():
    with open(filename, "r") as handle, open(output, "w") as outfile:
        loadedcsv = csv.reader(handle, delimiter="\t")
        next(loadedcsv, None)
        header = ["seqname", "start", "end", "chromatin"]
        writer = csv.DictWriter(outfile, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in loadedcsv:
            sets_of_seqname, chromatin = row
            try:
                match = regexp.match(sets_of_seqname)
                result = match.groupdict()
            except AttributeError:
                print(sets_of_seqname)
                raise
            result['chromatin'] = chromatin
            writer.writerow(result)

def main():
    format_kccell()

if __name__ == '__main__':
    main()