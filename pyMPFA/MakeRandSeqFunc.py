# C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8
import random
import Levenshtein
"""
**** MAKE RANDOM SEQUENCES FUNCTIONS
"""


def RandomSeqGenerator(count, length):
    y = 1
    seqList = []
    while y <= count:
        i = 1
        seq = ""
        while i <= length:
            seq += random.choice("ATGC")
            i += 1
        seqList.append(seq)
        y += 1
    return seqList


def CheckRandomSeq(seqList):
    chk = []
    for i in seqList:
        for y in seqList:
            if i != y:
                if Levenshtein.distance(i, y) <= int(len(i) * 0.35):
                    chk.append(False)
                else:
                    chk.append(True)
    return chk


def mainRandom(count, length):
    seq = RandomSeqGenerator(count, length)
    chk = CheckRandomSeq(seq)
    while False in chk:
        seq = RandomSeqGenerator(count, length)
        chk = CheckRandomSeq(seq)
        print("False in chk")
    print("Okay, I'm finding correct seq!")
    return seq
