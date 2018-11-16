#C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8

import param, csv, os, picks, regex, inspect
import SupportFunc as supp
from TripMain_0_2 import Pload
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from toolshed import nopen
from collections import Counter
from PairedEndFunc import reverseComplement
# import ReadIndexesFunc as rind
# import CollectBcMutFunc as colb
# import ReliableCombBcMutFunc as relc
# import WriteFunc as wrt
import Levenshtein
import matplotlib.pyplot as plt
import seaborn as sns
# import progressbar, csv, Levenshtein, regex, string, argparse, random, subprocess
# from collections import Counter
# from Bio import SeqIO
# from Bio.SeqIO.QualityIO import FastqGeneralIterator
# import itertools
# from itertools import islice
# from toolshed import nopen
supp.setup_logging()
# logger = logging.getLogger(__name__)


def WriteResultsToFile(resultDict, bcDict, seqDict, statDict, workdir, indexFile, customTxt=''):
    indexString = os.path.basename(indexFile).split(".")[0].split("_")[1]
    csvFile = os.path.join(workdir, "{}_barcode-mutation_count_{}.csv".format(indexString, customTxt))
    with open(csvFile, "wb") as handle:
        fieldnames = ['Barcode',
                    'Mutation',
                    'MutationCount',
                    'BCSequence',
                    'BCSequenceCount',
                    'MutatedBCassociatedWith',
                    'MutationVariants',
                    'Frequency',
                    'LostBarcodes',
                    'ProbableMutations',
                    'theirFrequency',
                    'HybridsPortion']
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for barcodeID in resultDict:
            BCSequence = ""
            BCSequenceCount = ""
            BCSequence = '\n'.join(str(x[0]) for x in bcDict[barcodeID]) + "\n"
            BCSequenceCount = '\n'.join(str(x[1]) for x in bcDict[barcodeID]) + "\n"
            TotalBCSequence = sum([x[1] for x in bcDict[barcodeID]])
            ProbableMutations = '\n'.join(str(k) for k in statDict[barcodeID]['mutations'])
            theirFrequency = '\n'.join(str(v) for v in statDict[barcodeID]['mutations'].values())
            HybridsPortion = str(statDict[barcodeID]['hybrids'] * 100) + "%"
            if barcodeID in seqDict[barcodeID]:
                MutVarForUnqBc = '\n'.join(str(x[0]) for x in seqDict[barcodeID][barcodeID])
                CountMutForUnqBc = '\n'.join(str(x[1]) for x in seqDict[barcodeID][barcodeID])
                TotalMutForUncBc = len(seqDict[barcodeID][barcodeID])
            else:
                MutVarForUnqBc = "none"
                CountMutForUnqBc = "none"
                TotalMutForUncBc = 1
            if len(resultDict[barcodeID]) > 2:
                wtd = resultDict[barcodeID][2]
                if len(wtd) != 0:
                    tmplist = []
                    for i in wtd.values():
                        tmplist.extend(i)
                    AssociatedMutations = '\n'.join(str(x[0]) for x in tmplist)
                    Frequency = '\n'.join(str(x[1]) for x in tmplist)
                    MutatedBarcodes = ""
                    for bc in wtd.keys():
                        countMut = len(wtd[bc])
                        if bc == wtd.keys()[len(wtd.keys())-1]:
                            MutatedBarcodes += str(bc)+"\n"*(countMut-1)
                        else:
                            MutatedBarcodes += str(bc)+"\n"*countMut
                else:
                    MutatedBarcodes = Frequency = AssociatedMutations = ""
            else:
                MutatedBarcodes = Frequency = AssociatedMutations = ""
            writer.writerow({
                'Barcode': barcodeID,
                'Mutation': resultDict[barcodeID][0],
                'MutationCount': resultDict[barcodeID][1],
                'BCSequence': BCSequence,
                'BCSequenceCount': BCSequenceCount,
                'MutatedBCassociatedWith': barcodeID+"\n"*TotalMutForUncBc+MutatedBarcodes,
                'MutationVariants': MutVarForUnqBc+"\n"+AssociatedMutations,
                'Frequency': CountMutForUnqBc+"\n"+Frequency,
                'LostBarcodes': TotalBCSequence-resultDict[barcodeID][1],
                'ProbableMutations': ProbableMutations,
                'theirFrequency': theirFrequency,
                'HybridsPortion': HybridsPortion})
    return os.path.basename(csvFile)


def mean(numbers):
    try:
        return float(sum(numbers)) / max(len(numbers), 1)
    except:
        supp.log_error("Undefined error in function {}".format(inspect.stack()[0][3]))


mapd, normd, exprd = "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_Genome_Mapping_A1-4/Dump", "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_Genome_Norm_A10-13/Dump", "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_Genome_Expr_A20-23/Dump"

namesList = ["18-1-1", "18-1-2", "18-2-1", "18-2-2"]

mlist, nlist, elist = ["m" + i for i in namesList], ["n" + i for i in namesList], ["e" + i for i in namesList]

bc, res = "bcDictPI", "resultDictPI"
pmi_syn = ["HexA", "Hsp70", "MtnA", "PCNA", "Pyk", "Tbp", "Promoterless"]
pmiDict = dict(zip(param.pmi, pmi_syn))

mapBCDict, normBCDict, mapResDict = {}, {}, {}
for item in namesList:
    # mapBCDict["m"+item] = Pload("m"+item+"_"+bc, mapd)
    normBCDict["n"+item] = Pload("n"+item+"_"+bc, normd)
    mapResDict["m"+item] = Pload("m"+item+"_"+res, mapd)

mainDict = {"map": mapResDict, "norm": normBCDict}
mainBioDict, mainMergeDict = {}, {}

try:
    for element in mainDict:
        supp.log_info(element)
        for replicate in mainDict[element]:
            supp.log_info(replicate)
            label, bio, tech = replicate.split('-')[0], replicate.split('-')[1], replicate.split('-')[2]
            for bio_replicate in [1, 2]:
                for pmi in pmiDict:
                    if int(bio) == bio_replicate:
                        supp.log_info("Technician replicate {} - {}: {}".format(replicate, pmiDict[pmi], len(mainDict[element][replicate][pmi])))
                        if int(tech) == 1:
                            if element == 'map':
                                mainBioDict[(label, bio, pmi)] = [k for k, v in mainDict[element][replicate][pmi].items() if k in mainDict[element][replicate[:-1] + '2'][pmi] and v[0] == mainDict[element][replicate[:-1] + '2'][pmi][k][0]]
                            else:
                                mainBioDict[(label, bio, pmi)] = [key for key in mainDict[element][replicate][pmi] if key in mainDict[element][replicate[:-1] + '2'][pmi]]
                            supp.log_info("Biological replicate {} - {}: {}".format(label + '-' + bio, pmiDict[pmi], len(mainBioDict[(label, bio, pmi)])))
except:
    supp.log_error("Unhandled error!")

try:
    for replicate in ['1', '2']:
        label = '18-' + replicate
        for pmi in pmiDict:
            mainMergeDict[(label, pmi)] = [k for k in mainBioDict[('m18', replicate, pmi)] if k in mainBioDict[('n18', replicate, pmi)]]
            supp.log_info("Merge count {} - {}: {}".format(label, pmiDict[pmi], len(mainMergeDict[(label, pmi)])))
except:
    supp.log_error("Unhandled error in step {} {}".format(replicate, pmiDict[pmi]))

mainMergeDictFile = os.path.join(mapd, "mainMergeDict.csv")
with open(mainMergeDictFile, "wb") as handle:
    fieldnames = ['Replicate 1', 'Replicate 2']
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    for pmi in pmiDict:
        writer.writerow({"Replicate 1": pmiDict[pmi], "Replicate 2": ""})
        one = mainMergeDict[('18-1', pmi)]
        two = mainMergeDict[('18-2', pmi)]
        if len(one) < len(two):
            diff = [""] * (len(two) - len(one))
            one.sort()
            one.extend(diff)
            iterable = range(len(two))
        else:
            diff = [""] * (len(one) - len(two))
            two.sort()
            two.extend(diff)
            iterable = range(len(one))
        for item in iterable:
            writer.writerow({"Replicate 1": one[item], "Replicate 2": two[item]})

# A10-13 normalization index
# indexList = {"n18-1-1":"AGTCGCCG", "n18-1-2":"TAAACATC", "n18-2-1":"ACAATTCG", "n18-2-2":"TACTTGTC"}

# A20-23 expression index
# indexList = {"e18-1-1":"GGTATGTT", "e18-1-2":"GAGGGACC", "e18-2-1":"TAGCTCTA", "e18-2-2":"TAATTGCG"}

# a14, a15, a16 = get_sequence_count
gatc_path = '/home/anton/data/R-script/R-counts/GATCs_dm6.txt'
gatcDict = {}
with open(gatc_path, 'r') as handle:
    for line in handle:
        motif = line.split()[1:4]
        if gatcDict.get(motif[0]) is None:
            gatcDict[motif[0]] = motif[1:]
        else:
            gatcDict[motif[0]].extend(motif[1:])

supp.log_info("\nMatch GATC's\n")
unmatched_bc_gatcs = {}
for exp in mainDict['map']:
    for pmi in mainDict['map'][exp]:
        unmatch_count = 0
        for k, v in mainDict['map'][exp][pmi].items():
            chrom, coord, direction = v[0][0], v[0][1], v[0][2]
            if direction == "-":
                coord_list = [coord, str(int(coord) + 1), str(int(coord) - 1),  str(int(coord) + 4), str(int(coord) - 4)]
                if sum([1 for x in coord_list if x not in gatcDict[chrom]]) == 5:
                    # unmatch_count += 1
                    if unmatched_bc_gatcs.get(exp) is None:
                        unmatched_bc_gatcs[exp] = [k]
                    else:
                        unmatched_bc_gatcs[exp].append(k)
            else:
                if coord not in gatcDict[chrom]:
                    # unmatch_count += 1
                    if unmatched_bc_gatcs.get(exp) is None:
                        unmatched_bc_gatcs[exp] = [k]
                    else:
                        unmatched_bc_gatcs[exp].append(k)
        # supp.log_info("{}, {}. All gatcs: {}, unmatched gatcs: {}".format(exp, pmi, len(mainDict['map'][exp][pmi]), unmatch_count))

# supp.log_info("\nAlign in one genome position in mainBioDict\n")
# for pmi in pmiDict:
#     for bcOne in mainBioDict[('m18', '1', pmi)]:

gatc_mismatch_file = os.path.join(mapd, "gatc_mismatch.fastq")
counter = 0

with open(gatc_mismatch_file, "w") as handle:
    for title, seq, qual in FastqGeneralIterator(nopen(picks.input_file)):
        counter += 1
        if counter % 1000 == 0:
            supp.log_info("Prepared {} reads".format(counter))
        for repl in param.indexList:
            expr_repl = regex.compile("^{}.*".format(param.indexList[repl]))
            match_repl = expr_repl.match(seq)
            if match_repl is not None:
                for barcode in unmatched_bc_gatcs[repl]:
                    expr = regex.compile("({})({}){{s<={}}}{}.*".format(param.indexList[repl].upper(), param.const_1.upper(), str(param.const_1Error), reverseComplement(barcode)))
                    match = regex.search(expr, seq)
                    if match is not None:
                        seqStr = [title, seq, qual]
                        handle.write("@{}\n{}\n+\n{}\n".format(*seqStr))


# ###################################### Lib 29-36
# indexList = {"m1": "CAAGATAA", "m2": "GGACAACG"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_BC_Mut_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_BC_Mut_mapping"

indexList = {"29-36_m1": "CAAGATAA", "29-36_m2": "GGACAACG"}
lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_3__bcmutProb_80/Lib_29-36/sample_S1_L001_R1_001_Lib_29-36_mapping/Dump"
workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_3__bcmutProb_80/Lib_29-36/sample_S1_L001_R1_001_Lib_29-36_mapping"

# indexList = {"29-36_m1": "CAAGATAA", "29-36_m2": "GGACAACG"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_1__bcmutProb_95/sample_20171127_20180510_Lib_29-36_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_1__bcmutProb_95/sample_20171127_20180510_Lib_29-36_mapping"

# indexList = {"29-36_m1": "CAAGATAA", "29-36_m2": "GGACAACG"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_0__bcmutProb_95_bcError_0/Lib_29-36/sample_S1_L001_R1_001_Lib_29-36_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_0__bcmutProb_95_bcError_0/Lib_29-36/sample_S1_L001_R1_001_Lib_29-36_mapping"

# indexList = {"29-36_m1": "CAAGATAA", "29-36_m2": "GGACAACG"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_0__bcmutProb_80_bcError_0/Lib_29-36/sample_S1_L001_R1_001_Lib_29-36_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_0__bcmutProb_80_bcError_0/Lib_29-36/sample_S1_L001_R1_001_Lib_29-36_mapping"

# indexList = {"29-36_m1": "CAAGATAA", "29-36_m2": "GGACAACG"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_0__bcmutProb_50_bcError_0/Lib_29-36/sample_S1_L001_R1_001_Lib_29-36_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_0__bcmutProb_50_bcError_0/Lib_29-36/sample_S1_L001_R1_001_Lib_29-36_mapping"

# indexList = {"29-36_m1": "CAAGATAA", "29-36_m2": "GGACAACG"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_2__bcmutProb_95_bcError_2_bmc_genomics/Lib_29-36/sample_S1_L001_R1_001_Lib_29-36_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_2__bcmutProb_95_bcError_2_bmc_genomics/Lib_29-36/sample_S1_L001_R1_001_Lib_29-36_mapping"

indexList = {"29_36_m1": "CAAGATAA", "29_36_m2": "GGACAACG"}
lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_3__bcmutProb_80/Lib_29-36/sample_S1_L001_R1_001_Lib_29-36_mapping/Dump"
workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_3__bcmutProb_80/Lib_29-36/sample_S1_L001_R1_001_Lib_29-36_mapping"

# ################################################## Lib 33-40
# indexList = {"33_40_m1": "AGCGAGCT", "33_40_m2": "CTGCACGT"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_1__bcmutProb_80/Lib_33-40/1_S1_L001_R1_001_Lib_33-40_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_1__bcmutProb_80/Lib_33-40/1_S1_L001_R1_001_Lib_33-40_mapping"

# indexList = {"33_40_m1": "AGCGAGCT", "33_40_m2": "CTGCACGT"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_3__bcmutProb_95/Lib_33-40/1_S1_L001_R1_001_Lib_33-40_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_3__bcmutProb_95/Lib_33-40/1_S1_L001_R1_001_Lib_33-40_mapping"

# indexList = {"33_40_m1": "AGCGAGCT", "33_40_m2": "CTGCACGT"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_3__bcmutProb_95_bcError_0/Lib_33-40/1_S1_L001_R1_001_Lib_33-40_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_3__bcmutProb_95_bcError_0/Lib_33-40/1_S1_L001_R1_001_Lib_33-40_mapping"

# indexList = {"33_40_m1": "AGCGAGCT", "33_40_m2": "CTGCACGT"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_0__bcmutProb_50_bcError_0/Lib_33-40/1_S1_L001_R1_001_Lib_33-40_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_0__bcmutProb_50_bcError_0/Lib_33-40/1_S1_L001_R1_001_Lib_33-40_mapping"

# indexList = {"33_40_m1": "AGCGAGCT", "33_40_m2": "CTGCACGT"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_0__bcmutProb_80_bcError_0/Lib_33-40/1_S1_L001_R1_001_Lib_33-40_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_0__bcmutProb_80_bcError_0/Lib_33-40/1_S1_L001_R1_001_Lib_33-40_mapping"

# indexList = {"33_40_m1": "AGCGAGCT", "33_40_m2": "CTGCACGT"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_0__bcmutProb_95_bcError_0/Lib_33-40/1_S1_L001_R1_001_Lib_33-40_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_0__bcmutProb_95_bcError_0/Lib_33-40/1_S1_L001_R1_001_Lib_33-40_mapping"

# indexList = {"33_40_m1": "AGCGAGCT", "33_40_m2": "CTGCACGT"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_1__bcmutProb_80/Lib_33-40_/1_S1_L001_R1_001_Lib_33-40_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_1__bcmutProb_80/Lib_33-40_/1_S1_L001_R1_001_Lib_33-40_mapping"

# indexList = {"33_40_m1": "AGCGAGCT", "33_40_m2": "CTGCACGT"}
# lib_mapping_dump = "/home/anton/backup/input/trip/RUN_2018-05-10/results/1_S1_L001_R1_001_Lib_33-40_mapping/Dump"
# workdir = "/home/anton/backup/input/trip/RUN_2018-05-10/results/1_S1_L001_R1_001_Lib_33-40_mapping"

control_barcodes = {"wt1": ["TTCCAAGTGCAGGTTAGGCG", "TTACGCAT"],
                    "wt2": ["TGTGTACGGCTTGCTCTCAA", "TTACGCAT"],
                    "dc3": ["GAGCCCGGATCCACTCCAAG", "TTAGCATG"],
                    "dc4": ["TGTCACGTCAGCTAACCCAC", "TTAGCATG"]}
lib_bcDict, lib_seqDict, lib_resultDict = {}, {}, {}
similarity = 1


def count_hybrids():
    for item in indexList:
        lib_bcDict[item] = Pload(item+"_bcDict", lib_mapping_dump)
        lib_seqDict[item] = Pload(item+"_seqDict", lib_mapping_dump)
        lib_resultDict[item] = Pload(item+"_resultDict", lib_mapping_dump)
        lib_statDict, hybrids_summary_stat = {}, []
        for bc in lib_resultDict[item]:
            main_mutation_seq = lib_resultDict[item][bc][0]
            main_mutation_count = lib_resultDict[item][bc][1]
            temp_mutation, temp_statistic = [], {}
            for k in lib_seqDict[item][bc]:
                temp_mutation.extend(lib_seqDict[item][bc][k])
            # '''
            # x[0]: this is sequence of mutation
            # x[1]: this is count of this mutation
            # The element of temp_mutation is tuple: ('mutation', count)
            # '''
            temp_mutation = dict(Counter([x[0] for x in temp_mutation for _ in range(x[1])]))
            temp_statistic[main_mutation_seq] = temp_mutation[main_mutation_seq]
            for mutations in temp_mutation:
                if mutations != main_mutation_seq:
                    distance = Levenshtein.distance(main_mutation_seq, mutations)
                    if distance <= similarity:
                        temp_statistic[main_mutation_seq] += temp_mutation[mutations]
                    else:
                        temp_statistic[mutations] = temp_mutation[mutations]
            hybrid_portion = round(1 - (float(max(temp_statistic.values())) / float(sum(temp_statistic.values()))), 6)
            hybrids_summary_stat.append(hybrid_portion)
            lib_statDict[bc] = dict(hybrids=hybrid_portion, mutations=temp_statistic)
        hybrids_summary = str(round(mean(hybrids_summary_stat) * 100, 4)) + "%"
        # plt.figure(figsize=(10, 10))
        # plt.title("Histogram density plot for hybrids portion in {}".format(item))
        # plt.xlabel("Hybrids portion in pct")
        # plt.ylabel("Portion of hybrids portion from all values")
        # x = [x * 100 for x in hybrids_summary_stat]
        # sns.kdeplot(x)
        # plt.savefig(os.path.join(lib_mapping_dump, "Histogram density plot for hybrids portion in {}.pdf".format(item)), fmt='pdf')
        index = indexList[item]
        indexFile = os.path.join(workdir, item, "index_{}.fastq".format(index.upper()))
        csv_file = WriteResultsToFile(lib_resultDict[item], lib_bcDict[item], lib_seqDict[item], lib_statDict, os.path.join(workdir, item), indexFile, customTxt="with_hybrids_pct")
        supp.log_info("Report write to {}\nHybrids pct: {}".format(csv_file, hybrids_summary))

count_hybrids()
