#C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8

import os, pickle, param, picks, itertools
from SupportFunc import GetTotalSeqRecords, simpleWrite, head
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
# import progressbar, csv, Levenshtein, regex, string, argparse, random, subprocess
# from collections import Counter
# from Bio import SeqIO
# from Bio.SeqIO.QualityIO import FastqGeneralIterator
# import itertools
# from itertools import islice
# from toolshed import nopen

"""
**** SAVE/LOAD FUNCTIONS
"""

def SaveDictToPy(dictVar, filename):
    with open(filename + ".py", "wb") as handle:
        for k, v in dictVar.items():
            if type(v) == str:
                handle.write(str(k) + " = '" + str(v) + "'\n")
            else:
                handle.write(str(k) + " = " + str(v) + "\n")

def Pdump(obj, name, folder):
    locationObj = os.path.join(folder, name)
    with open(locationObj + ".pickle", 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

def Pload(name, folder):
    locationObj = os.path.join(folder, name)
    with open(locationObj + '.pickle', 'rb') as handle:
        unserialized_data = pickle.load(handle)
        return unserialized_data

################
# MUTATION

workdir = picks.workdir
indexName = {k:os.path.join(workdir, k, "index_{}.fastq".format(v)) for k,v in param.indexList.items()}
dump = picks.PdumpDir
FRC = {k:GetTotalSeqRecords(v) for k,v in indexName.items()}

# load bcDict from pickle
simpleWrite(FRC, dump, "frc.txt")
n1 = {k:sum([i[1] for i in v]) for k, v in Pload(FRC.keys()[3] + "_bcDict", dump).items()}
n2 = {k:sum([i[1] for i in v]) for k, v in Pload(FRC.keys()[2] + "_bcDict", dump).items()}
e1 = {k:sum([i[1] for i in v]) for k, v in Pload(FRC.keys()[1] + "_bcDict", dump).items()}
e2 = {k:sum([i[1] for i in v]) for k, v in Pload(FRC.keys()[0] + "_bcDict", dump).items()}

# Aligning dict between them
data = {"n1":n1, "n2":n2, "e1":e1, "e2":e2}
variations = list(itertools.permutations(["n1", "n2", "e1", "e2"], 2))
for comb in variations:
    for i in data[comb[0]]:
        if i not in data[comb[1]].keys():
            data[comb[1]][i] = "NA"
for k, v in data.items():
    simpleWrite(v, dump, k + ".txt")

# load resultDict from pickle for bc_mut combinations 0.95/0.8 CutOff
# m1 = Pload(FRC.keys()[0] + "_resultDict", dump)
# m2 = Pload(FRC.keys()[1] + "_resultDict", dump)

m1 = {k:v[0] for k, v in Pload(FRC.keys()[0] + "_resultDict", dump).items()}
m2 = {k:v[0] for k, v in Pload(FRC.keys()[1] + "_resultDict", dump).items()}

# merge mapping reads
mp1 = [(k, v) for k, v in m1.items()]
mp2 = [(k, v) for k, v in m2.items()]

mdata = {"m1":mp1, "m2":mp2}
mvariations = list(itertools.permutations(["m1", "m2"], 2))
mapRes = {}
for comb in mvariations:
    for i in mdata[comb[0]]:
        if i in mdata[comb[1]]:
            if not mapRes.get(i[0]):
                mapRes[i[0]] = i[1]
simpleWrite(mapRes, dump, "mp.txt")

# make Venn diagramm for mapping replicates 0.95/0.8 CutOff
plt.figure(figsize=(10,10))
plt.title("Venn Diagram for mapping reads. 0.95 CutOff")
venn2([set(mp1), set(mp2)], set_labels = ("Mapping replicate #1", "Mapping replicate #2"))
plt.savefig(os.path.join(dump, "venn_mapping_R1_R2_95.pdf"), fmt='pdf')

################
# GENOME
normdir = "/home/anton/backup/input/trip/RUN_2017-11-27/results/repeat/sample_S1_L001_R1_001_norm"
exprdir = "/home/anton/backup/input/trip/RUN_2017-11-27/results/repeat/sample_S1_L001_R1_001_expr"
mapdir = "/home/anton/backup/input/trip/RUN_2017-11-27/results/repeat/sample_S1_L001_R1_001_map"
normDict = {"n18-1-1":"AGTCGCCG", "n18-1-2":"TAAACATC", "n18-2-1":"ACAATTCG", "n18-2-2":"TACTTGTC"}
exprDict = {"e18-1-1":"GGTATGTT", "e18-1-2":"GAGGGACC", "e18-2-1":"TAGCTCTA", "e18-2-2":"TAATTGCG"}
mapDict = {"m18-1-1":"TTCGGAGT", "m18-1-2":"ACTCATTT", "m18-2-1":"GGGATCCG", "m18-2-2":"TCAAGCAA"}
indexNameE = {k:os.path.join(exprdir, k, "index_{}.fastq".format(v)) for k,v in exprDict.items()}
indexNameN = {k:os.path.join(normdir, k, "index_{}.fastq".format(v)) for k,v in normDict.items()}
indexNameM = {k:os.path.join(mapdir, k, "index_{}.fastq".format(v)) for k,v in mapDict.items()}
indexNameMF = {k:os.path.join(mapdir, k, "filt_index_{}.fastq".format(v)) for k,v in mapDict.items()}

dumpE = os.path.join(exprdir, "Dump")
dumpN = os.path.join(normdir, "Dump")
dumpM = os.path.join(mapdir, "Dump")

FRCE = {k:GetTotalSeqRecords(v) for k,v in indexNameE.items()}
FRCN = {k:GetTotalSeqRecords(v) for k,v in indexNameN.items()}

FRCM = {k:GetTotalSeqRecords(v) for k,v in indexNameM.items()}
FRCMF = {"f"+k:GetTotalSeqRecords(v) for k,v in indexNameMF.items()}

# load bcDict from pickle
simpleWrite(FRCE, dumpE, "frce.txt")
simpleWrite(FRCN, dumpN, "frcn.txt")
simpleWrite(FRCM, dumpM, "frcm.txt")
simpleWrite(FRCMF, dumpM, "frcmf.txt")

popName = {"pop1":["e18-1-1", "e18-1-2", "n18-1-1", "n18-1-2"], "pop2":["e18-2-1", "e18-2-2", "n18-2-1", "n18-2-2"]}


data = {}
for p in popName:
    data[p] = {}
    for exp in popName[p]:
        data[p][exp] = {}
        if exp[:1] == "e":
            dump = dumpE
        else:
            dump = dumpN
        tDict = Pload(str(exp) + "_bcDictPI", dump)
        for pI in param.pmi:
            data[p][exp][pI] = {k:sum([i[1] for i in v]) for k, v in tDict[pI].items()}

# Aligning dict between them
for p in popName:
    variations = list(itertools.permutations(popName[p], 2))
    for comb in variations:
        for pi in param.pmi:
            for i in data[p][comb[0]][pi]:
                if i not in data[p][comb[1]][pi].keys():
                    data[p][comb[1]][pi][i] = 0

for p in popName:
    for i in data[p]:
        for y in data[p][i]:
            if i[:1] == "e":
                dump = dumpE
            else:
                dump = dumpN
            simpleWrite(data[p][i][y], dump, "{}_{}.txt".format(i, y))

# make Venn diagramm for technician replicate

# plt.figure(figsize=(10,10))
# plt.title("Venn Diagram for mapping reads. 0.95 CutOff")
# venn2([set(mp1), set(mp2)], set_labels = ("Mapping replicate #1", "Mapping replicate #2"))
# plt.savefig(os.path.join(dump, "venn_mapping_R1_R2_95.pdf"), fmt='pdf')
################
# MAPPING GENOME

mapName = {"pop1":["m18-1-1", "m18-1-2"], "pop2":["m18-2-1", "m18-2-2"]}
mapdata = {}
for p in mapName:
    mapdata[p] = {}
    for mpp in mapName[p]:
        mapdata[p][mpp] = {}
        tDict = Pload(str(mpp) + "_bcDictPI", dumpM)
        for pI in param.pmi:
            mapdata[p][mpp][pI] = {k:sum([i[1] for i in v]) for k, v in tDict[pI].items()}

for p in mapName:
    for i in mapdata[p]:
        for y in mapdata[p][i]:
            simpleWrite(mapdata[p][i][y], dumpM, "{}_{}.txt".format(i, y))

################
# GENOME v.2
normdir = "/home/anton/backup/input/trip/RUN_2017-11-27/results/repeat/sample_S1_L001_R1_001_norm"
exprdir = "/home/anton/backup/input/trip/RUN_2017-11-27/results/repeat/sample_S1_L001_R1_001_expr"
normDict = {"n18-1-1":"AGTCGCCG", "n18-1-2":"TAAACATC", "n18-2-1":"ACAATTCG", "n18-2-2":"TACTTGTC"}
exprDict = {"e18-1-1":"GGTATGTT", "e18-1-2":"GAGGGACC", "e18-2-1":"TAGCTCTA", "e18-2-2":"TAATTGCG"}
indexNameE = {k:os.path.join(exprdir, k, "index_{}.fastq".format(v)) for k,v in exprDict.items()}
indexNameN = {k:os.path.join(normdir, k, "index_{}.fastq".format(v)) for k,v in normDict.items()}

dumpE = os.path.join(exprdir, "Dump")
dumpN = os.path.join(normdir, "Dump")

FRCE = {k:GetTotalSeqRecords(v) for k,v in indexNameE.items()}
FRCN = {k:GetTotalSeqRecords(v) for k,v in indexNameN.items()}

# load bcDict from pickle
simpleWrite(FRCE, dumpE, "frce.txt")
simpleWrite(FRCN, dumpN, "frcn.txt")

popName = {"p18-1":["e18-1-1", "e18-1-2", "n18-1-1", "n18-1-2"], "p18-2":["e18-2-1", "e18-2-2", "n18-2-1", "n18-2-2"]}
expName = []
# pmi = ['AGCTC', 'ACGTA', 'CTGCT', 'AGTCA', 'TCAAA', 'TTGAG', 'TCGCT']

data = {}
for p in popName:
    data[p] = {}
    for exp in popName[p]:
        data[p][exp] = {}
        if exp[:1] == "e":
            dump = dumpE
        else:
            dump = dumpN
        tDict = Pload(str(exp) + "_bcDictPI", dump)
        for pI in param.pmi:
            data[p][exp][pI] = {k:sum([i[1] for i in v]) for k, v in tDict[pI].items()}

# merge and sum dict
dataMerged = {}
for p in popName:
    for i in ["e", "n"]:
        for exp in popName[p]:
            if exp[:1] == i:
                if exp[-1:] == "1":
                    expm = exp[0:len(exp)-2]
                    dataMerged[expm] = {}
                    techRepl = expm + "-2"
                    for pi in param.pmi:
                        x = data[p][exp][pi]
                        y = data[p][techRepl][pi]
                        dataMerged[expm][pi] = { k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y) }
# Count bc in every PI
countDict = {}
for exp in dataMerged:
    countDict[exp] = {}
    for pi in dataMerged[exp]:
        countDict[exp][pi] = len(dataMerged[exp][pi])

for exp in countDict:
    simpleWrite(countDict[exp], dumpE, "{}_count.txt".format(exp))


# Aligning dict between them
variations = list(itertools.permutations(dataMerged.keys(), 2))
for comb in variations:
    for pi in param.pmi:
        for i in dataMerged[comb[0]][pi]:
            if i not in dataMerged[comb[1]][pi].keys():
                dataMerged[comb[1]][pi][i] = 0

# Write to csv file and import to R
for exp in dataMerged:
    for pi in dataMerged[exp]:
        if exp[:1] == "e":
            dump = dumpE
        else:
            dump = dumpN
        simpleWrite(dataMerged[exp][pi], dump, "{}_{}.txt".format(exp, pi))