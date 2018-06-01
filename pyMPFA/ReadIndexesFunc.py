#C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8
import regex, progressbar, random
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
    t=0.0
    with open(indexFile, "w") as handle:
        for title, seq, qual in FastqGeneralIterator(nopen(input_file)):
            bar.update(t)
            t+=1
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
    t=0
    randomSeq = [random.random() < probability for x in xrange(records)]
    with open(indexFileRand, "w") as handle:
        for title, seq, qual in FastqGeneralIterator(nopen(indexFile)):
            if randomSeq[t]:
                handle.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))
            bar.update(t)
            t+=1
    bar.finish()

def filterShadyReads(indexFile, reFilter, indexFiltFile):
    regFilter = regex.compile(reFilter)
    records = GetTotalSeqRecords(indexFile)
    bar = progressbar.ProgressBar(maxval=records, widgets=[progressbar.Bar(left='<', marker='.', right='>')]).start()
    t=0.0
    with open(indexFiltFile, "w") as handle:
        for title, seq, qual in FastqGeneralIterator(nopen(indexFile)):
            bar.update(t)
            t+=1
            match = regex.search(regFilter, seq)
            if match is None:
                seqStr = [title, seq, qual]
                handle.write("@{}\n{}\n+\n{}\n".format(*seqStr))
    bar.finish()
    filt_records = GetTotalSeqRecords(indexFiltFile)
    output = (records, filt_records)
    return output