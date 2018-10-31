#!/usr/bin/env python
# encoding: utf-8
import regex
import progressbar
import random
import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from toolshed import nopen
from SupportFunc import GetTotalSeqRecords
"""
**** READ INDEXES FUNCTIONS
"""


def SplitFastqByIndexes(input_file, indexFile, index, indexError, const_1, const_1Error, regExpIndex, no_trim):
    # print("    Processing index: '{}' ...".format(index))
    regIndex = regex.compile(regExpIndex[index])
    records = GetTotalSeqRecords(input_file)
    bar = progressbar.ProgressBar(maxval=records, widgets=[progressbar.Bar(left='<', marker='.', right='>')]).start()
    t = 0
    with open(indexFile, "w") as handle:
        for title, seq, qual in FastqGeneralIterator(nopen(input_file)):
            bar.update(t)
            t += 1
            match = regex.search(regIndex, seq)
            if match is not None:
                if no_trim:
                    seqStr = [title, seq, qual]
                else:
                    seqStr = [title, seq[len(match[0]):], qual[len(match[0]):]]
                handle.write("@{}\n{}\n+\n{}\n".format(*seqStr))
    bar.finish()


def RandomReadIndexes(indexFile, indexFileRand, probability):
    # print("    Random read indexFile: '{}' ...".format(os.path.basename(indexFile)))
    records = GetTotalSeqRecords(indexFile)
    bar = progressbar.ProgressBar(maxval=records, widgets=[progressbar.Bar(left='<', marker='.', right='>')]).start()
    t = 0
    randomSeq = [random.random() < probability for x in xrange(records)]
    with open(indexFileRand, "w") as handle:
        for title, seq, qual in FastqGeneralIterator(nopen(indexFile)):
            if randomSeq[t]:
                handle.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))
            bar.update(t)
            t += 1
    bar.finish()


def filterShadyReads(indexFile, reFilter, indexFiltFile):
    regFilter = regex.compile(reFilter)
    records = GetTotalSeqRecords(indexFile)
    bar = progressbar.ProgressBar(maxval=records, widgets=[progressbar.Bar(left='<', marker='.', right='>')]).start()
    t = 0
    with open(indexFiltFile, "w") as handle:
        for title, seq, qual in FastqGeneralIterator(nopen(indexFile)):
            bar.update(t)
            t += 1
            match = regex.search(regFilter, seq)
            if match is None:
                seqStr = [title, seq, qual]
                handle.write("@{}\n{}\n+\n{}\n".format(*seqStr))
    bar.finish()
    filt_records = GetTotalSeqRecords(indexFiltFile)
    output = (records, filt_records)
    return output


def SplitFastqByIndexesPaired(index, indexFile, options):
    _filename, _ext = os.path.splitext(os.path.basename(indexFile))
    indexFile_rev = os.path.join(os.path.dirname(indexFile), "{}_rev{}".format(_filename, _ext))
    regIndex = regex.compile(options["regExpIndex"][index])
    _index_titles = set()
    with open(indexFile, "w") as f_handle:
        for f_title, f_seq, f_qual in FastqGeneralIterator(nopen(options["r1"])):
            match = regex.search(regIndex, f_seq)
            if match is not None:
                _name, _idx = f_title.split()
                _index_titles.add(_name)
                if options["no_trim_index"]:
                    seqStr = [f_title, f_seq, f_qual]
                else:
                    seqStr = [f_title, f_seq[len(match[0]):], f_qual[len(match[0]):]]
                f_handle.write("@{}\n{}\n+\n{}\n".format(*seqStr))
    with open(indexFile_rev, "w") as r_handle:
        for r_title, r_seq, r_qual in FastqGeneralIterator(nopen(options["r2"])):
            _r_name, _r_idx = r_title.split()
            if _r_name in _index_titles:
                seqStr = [r_title, r_seq, r_qual]
                r_handle.write("@{}\n{}\n+\n{}\n".format(*seqStr))
    return {"fwd": indexFile, "rev": indexFile_rev}
