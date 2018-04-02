#C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8

import os, pickle. param
import SupportFunc as supp
from TripMain_0_2 import Pload, Pdump

# import progressbar, csv, Levenshtein, regex, string, argparse, random, subprocess
# from collections import Counter
# from Bio import SeqIO
# from Bio.SeqIO.QualityIO import FastqGeneralIterator
# import itertools
# from itertools import islice
# from toolshed import nopen

mapd, normd, exprd = "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_Genome_Mapping_A1-4/Dump", "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_Genome_Norm_A10-13/Dump", "/home/anton/backup/input/trip/RUN_2017-11-27/results/sample_S1_L001_R1_001_Genome_Expr_A20-23/Dump"

namesList = ["18-1-1", "18-1-2", "18-2-1", "18-2-2"]

mlist, nlist, elist = ["m" + i for i in namesList], ["n" + i for i in namesList], ["e" + i for i in namesList]

bc, res = "bcDictPI", "resultDictPI"

mapBCDict, normBCDict, mapResDict = {}, {}, {}
for item in namesList:
    mapBCDict["m"+item] = Pload("m"+item+"_"+bc, mapd)
    normBCDict["n"+item] = Pload("n"+item+"_"+bc, normd)
    mapResDict["m"+item] = Pload("m"+item+"_"+res, mapd)

# A10-13 normalization index
# indexList = {"n18-1-1":"AGTCGCCG", "n18-1-2":"TAAACATC", "n18-2-1":"ACAATTCG", "n18-2-2":"TACTTGTC"}

#A20-23 expression index
# indexList = {"e18-1-1":"GGTATGTT", "e18-1-2":"GAGGGACC", "e18-2-1":"TAGCTCTA", "e18-2-2":"TAATTGCG"}

a14, a15, a16 = GetTotalSeqRecords