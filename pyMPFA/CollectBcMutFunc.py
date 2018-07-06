#C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8
import regex, progressbar, os
import multiprocessing as mp
import subprocess
import SupportFunc as supp
from toolshed import nopen
from collections import Counter
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from PairedEndFunc import reverseComplement
from CheckBarcodesFunc import mainCheckBarcodeInDict, checkPMI
from ReliableCombBcMutFunc import SelectionReliableBarcode, SelectionReliableBarcodeMutation, SelectionReliableBarcodeGenome
"""
**** COLLECT BARCODE OR BARCODE/MUTATION FUNCTIONS
"""

def CollectBarcodeMutation(indexFile, barcodeLength, mutationLength, readsValue, barcodeError, const_2, const_3, const_2Error, const_3Error, regExpBcMut, reverseBC):
    # print("Start collect barcodes...\nI'm using next regular exxpression for search: {}".format(regExpBcMut))
    bcList, bcMutList = [], []
    non_matched_reads = os.path.join(os.path.dirname(indexFile), "undef_{}".format(indexFile))
    records = supp.GetTotalSeqRecords(indexFile)
    bar = progressbar.ProgressBar(maxval=records, widgets=[progressbar.Bar(left='<', marker='.', right='>')]).start()
    t=0.0
    expr = regex.compile(regExpBcMut)
    with nopen(indexFile) as handle:
        with open(non_matched_reads, "wb") as undef:
            for seq_record in SeqIO.parse(handle, "fastq"):
                bar.update(t)
                t+=1
                match = expr.match(str(seq_record.seq))
                if match is not None:
                    if int(barcodeLength*0.9) <= len(match.group("barcode")) <= int(barcodeLength*1.1) and len(match.group("mutation")) == mutationLength:
                        if "N" not in match.group("barcode").upper() and "N" not in match.group("mutation").upper():
                            if reverseBC:
                                bcList.append(reverseComplement(match.group("barcode")))
                                bcMutList.append((reverseComplement(match.group("barcode")), reverseComplement(match.group("mutation"))))
                            else:
                                bcList.append(match.group("barcode"))
                                bcMutList.append((match.group("barcode"), match.group("mutation")))
                else:
                    SeqIO.write(seq_record, undef, "fastq")
    bar.finish()
    bcCount = Counter(bcList)
    bcMutCount = dict(Counter(bcMutList))
    # print("        Total founded barcodes: {} items.\nAnd total combinations barcode-mutation: {} items".format(len(bcCount), len(bcMutCount)))
    bcDict = SelectionReliableBarcode(bcCount, readsValue, barcodeError)
    if len(bcDict) <= 10**5: 
        # print("        Checking barcodes ... Estimated time: ~ {}".format(supp.EstimateCalculationTime(bcDict)))
        mainCheckBarcodeInDict(bcDict, barcodeError)
    seqDict = SelectionReliableBarcodeMutation(bcMutCount, bcDict)
    return (bcDict, seqDict)

def CollectBarcodeMutationGenome(indexFile, barcodeLength, readsValue, barcodeError, const_2, const_2Error, regExpBcMut, merge_indexes, pmi, pmiLength, pmiSubst, reverseBC, bwIndex, rfplIndex):
    bcListDict, bcGenomeListDict = {p:[] for p in pmi}, {p:{} for p in pmi}
    bcDictPI, seqDictPI = {}, {}
    bcDictTmp, alignDict = {}, {}
    expr = regex.compile(regExpBcMut)
    for pmiItem in pmi:
        pmiWD = os.path.join(os.path.dirname(indexFile), pmiItem)
        pmiFile = os.path.join(pmiWD, pmiItem + ".fastq")
        countAllRds, countConst3, countLess20 = 0, 0, 0
        supp.LogInfo("  Collect from promotor index {}".format(pmiItem))
        if not os.path.exists(pmiWD): os.makedirs(pmiWD)    
        with nopen(indexFile) as handle, open(pmiFile, "w") as pmiHandle:      
            for title, seq, qual in FastqGeneralIterator(handle):
                match = expr.match(seq)
                if match is not None:
                    countAllRds += 1
                    if int(barcodeLength*0.9) <= len(match.group("barcode")) <= int(barcodeLength*1.1) and len(match.group("pIndex")) == pmiLength:
                        if "N" not in match.group("barcode").upper() and "N" not in match.group("pIndex").upper():
                            if match.group("genome1") is not None:
                                g = match.group("genome1")
                            else:
                                g = match.group("genome2")
                                countConst3 += 1
                            if reverseBC:
                                b = reverseComplement(match.group("barcode"))
                                p = reverseComplement(match.group("pIndex"))
                            else:
                                b = match.group("barcode")
                                p = match.group("pIndex")
                            matchedPMI = checkPMI(p, pmi, pmiSubst)
                            if matchedPMI in pmi and matchedPMI == pmiItem:
                                if len(g) >= 20:
                                    bcListDict[matchedPMI].append(b)
                                    bcGenomeListDict[matchedPMI][title.split(" ")[0]] = b
                                    FstGenomeCoord = len(b) + 16
                                    SndGenomeCoord = FstGenomeCoord + len(g)
                                    seqStr = [title, seq[FstGenomeCoord:SndGenomeCoord], qual[FstGenomeCoord:SndGenomeCoord]]
                                    pmiHandle.write("@{}\n{}\n+\n{}\n".format(*seqStr))
                                else:
                                    countLess20 += 1
        supp.LogInfo("      All reads: {} // Reads with GENOME length less 20 bp: {} // Reads with CONST_3 (20bp): {}".format(countAllRds, countLess20, countConst3))
        bcCount = Counter(bcListDict[pmiItem])
        bcDictTmp[pmiItem] = SelectionReliableBarcode(bcCount, readsValue, barcodeError)
        if len(bcDictTmp[pmiItem]) <= 10**5 and len(bcDictTmp[pmiItem]) > 0: mainCheckBarcodeInDict(bcDictTmp[pmiItem], barcodeError)
        alignDict[pmiItem], bwtAlignerDict = AlignReadsToGenome(pmiFile, pmiWD, pmiItem, bwIndex, rfplIndex)
        supp.makeStatFromBowtieAlign(bwtAlignerDict)
        mergeList = [(bcGenomeListDict[pmiItem][k], alignDict[pmiItem][k]) for k in bcGenomeListDict[pmiItem] if k in alignDict[pmiItem]]
        mergeListCount = dict(Counter(mergeList))
        bcDictPI[pmiItem], seqDictPI[pmiItem] = SelectionReliableBarcodeGenome(mergeListCount, bcDictTmp[pmiItem])
    return bcDictPI, seqDictPI

def AlignReadsToGenome(pmiFile, pmiWD, pmiItem, bwIndex, rfplIndex):
    cpus = mp.cpu_count()
    gDict, tmpbwt, bwtAlignerDict = {}, {}, {}
    filteredFQ, filterSam = os.path.join(pmiWD, "filt_" + os.path.basename(pmiFile)), os.path.join(pmiWD, os.path.basename(pmiFile).split('.')[0] + "_filter.sam")
    bamOut, bedOut = os.path.join(pmiWD, "dmel.bam"), os.path.join(pmiWD, "dmel_uniq.bed")
    bwtCore = "(bowtie2 -k 3 -p " + str(cpus) + " -t --phred33 --local -x "
    filterCmd = "--un " + " ".join([filteredFQ, rfplIndex, pmiFile]) + " -S " + filterSam + ")"
    alignCmd = " ".join([bwIndex, filteredFQ]) + " | samtools view -bS - -o " + bamOut + ")"
    samtoolsFilter = "(samtools sort -@ " + str(cpus) + " -o " + os.path.join(pmiWD, "dmel.bam") + " dmel-sort | samtools view -h -F 4 - | samtools view -S -q 25 -F 256 - | grep -Pv 'XS:i' | sam2bed > " + bedOut + ") 2>/dev/null"
    runFilter = subprocess.Popen(bwtCore + filterCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (tmpbwt[pmiItem  + '_outF'], tmpbwt[pmiItem  + '_errF']) = runFilter.communicate()
    runAlign = subprocess.Popen(bwtCore + alignCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (tmpbwt[pmiItem  + '_outAlg'], tmpbwt[pmiItem  + '_errAlg']) = runAlign.communicate()
    for k, v in tmpbwt.items():
        if v is not None: bwtAlignerDict[k] = v.splitlines()
    subprocess.call(samtoolsFilter, shell=True)
    with open(bedOut) as dmel:
        for line in dmel:
            chrName, start, end, title, strand = [line.split()[x] for x in [0, 1, 2, 3, 5]]
            if strand == "-":
                coord = end
            else:
                coord = start
            gDict[title] = (chrName, coord, strand)
    return gDict, bwtAlignerDict

def CollectBarcodeGenome(indexFile, barcodeLength, readsValue, barcodeError, const_2, const_2Error, regExpBc, mergeBC, reverseBC, pmi, pmiLength, pmiSubst):
    bcListDict = {p:[] for p in pmi}
    bcDictPI = {}
    records = supp.GetTotalSeqRecords(indexFile)
    bar = progressbar.ProgressBar(maxval=records, widgets=[progressbar.Bar(left='<', marker='.', right='>')]).start()
    t=0.0
    expr = regex.compile(regExpBc)
    with nopen(indexFile) as handle:
        for seq_record in SeqIO.parse(handle, "fastq"):
            bar.update(t)
            t+=1
            match = expr.match(str(seq_record.seq))
            if match is not None:
                if int(barcodeLength*0.9) <= len(match.group("barcode")) <= int(barcodeLength*1.1) and len(match.group("pIndex")) == pmiLength:
                    if "N" not in match.group("barcode").upper() and "N" not in match.group("pIndex").upper():
                        if reverseBC:
                            b = reverseComplement(match.group("barcode"))
                            p = reverseComplement(match.group("pIndex"))
                        else:
                            b = match.group("barcode")
                            p = match.group("pIndex")
                        matchedPMI = checkPMI(p, pmi, pmiSubst)
                        if matchedPMI in bcListDict:
                            bcListDict[matchedPMI].append(b)
    bar.finish()
    if mergeBC:
        return bcListDict
    for pI in bcListDict:
        bcCount = Counter(bcListDict[pI])
        bcDictPI[pI] = SelectionReliableBarcode(bcCount, readsValue, barcodeError)
        if len(bcDictPI[pI]) <= 10**5 and len(bcDictPI[pI]) > 0: 
            # print("        Checking barcodes for promotor index {} ... Estimated time: ~ {}".format(pI, supp.EstimateCalculationTime(bcDictPI[pI])))
            mainCheckBarcodeInDict(bcDictPI[pI], barcodeError)
    return bcDictPI

def CollectBarcode(indexFile, barcodeLength, readsValue, barcodeError, const_2, const_2Error, regExpBc, mergeBC, reverseBC):
    bcList = []
    records = supp.GetTotalSeqRecords(indexFile)
    bar = progressbar.ProgressBar(maxval=records, widgets=[progressbar.Bar(left='<', marker='.', right='>')]).start()
    t=0.0
    expr = regex.compile(regExpBc)
    with nopen(indexFile) as handle:
        for seq_record in SeqIO.parse(handle, "fastq"):
            bar.update(t)
            t+=1
            match = expr.match(str(seq_record.seq))
            if match is not None:
                if int(barcodeLength*0.9) <= len(match.group("barcode")) <= int(barcodeLength*1.1):
                    if "N" not in match.group("barcode").upper():
                        if reverseBC:
                            bcList.append(reverseComplement(match.group("barcode")))
                        else:
                            bcList.append(match.group("barcode"))
    bar.finish()
    if mergeBC:
        return bcList
    bcCount = Counter(bcList)
    bcDict = SelectionReliableBarcode(bcCount, readsValue, barcodeError)
    if len(bcDict) <= 10**5: 
        # print("        Checking barcodes ... Estimated time: ~ {}".format(supp.EstimateCalculationTime(bcDict)))
        mainCheckBarcodeInDict(bcDict, barcodeError)
    return bcDict