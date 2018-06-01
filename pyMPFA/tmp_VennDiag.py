#!/usr/bin/env python
import numpy as np
probability = np.arange(0.1, 1, 0.1)
attempt = 10
indexFile = "/home/anton/data/R-script/trip/python/randomRead2/index_NNNNNNNN.fastq"
indexFileRand = os.path.join(os.path.dirname(indexFile), "random_"+os.path.basename(indexFile))
reportOut = os.path.join(os.path.dirname(indexFile), "report.txt")
report={}

def randomread():
    for prob in probability:
        i = 0
        records = 0
        print("Start count {} probability".format(str(prob)))
        while i < attempt:
            print("Iterations number {}...".format(str(i)))
            RandomReadIndexes(indexFile, indexFileRand, prob)
            bcList = CollectBarcode(indexFileRand, param.barcodeLength, param.readsValue, param.barcodeError, param.const_2, param.const_2Error, param.regExpBc, args.merge_indexes, args.reverse_barcode)
            bcCount = Counter(bcList)
            bcCountSorted = sorted(bcCount.items(), key=lambda x:x[1], reverse=True)
            bcLengthMax = max(Counter([len(x[0]) for x in bcCountSorted]).items(), key=lambda x:x[1])
            bcLengthRatio = float(bcLengthMax[1]) / float(len(bcCountSorted))
            if bcLengthRatio >= 0.95: bcCountSorted = [x for x in bcCountSorted if len(x[0]) == bcLengthMax[0]]
            bcCountSorted = [x for x in bcCountSorted if x[1] > param.readsValue ]
            records += len(bcCountSorted)
            i+=1
        report[str(prob)] = str(float(records) / float(attempt))

def writeResults():
    with open(reportOut, "wb") as handle:
        fieldnames = ["Probability", "BarcodeCounts"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=';')
        writer.writeheader()
        for r in report:
            writer.writerow({
                "Probability": r,
                "BarcodeCounts": report[r]})

if __name__ == '__main__':
    randomread()
    writeResults()


from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles

vennDict = {}
indexNames = {k: v for (k, v) in zip(param.indexList, ['A'+str(x) for x in range(10, 17)])}
for index in param.indexList:
    indexFile = os.path.join(workdir, "index_{}.fastq".format(index.upper()))
    if not os.path.exists(indexFile) or os.stat(indexFile).st_size == 0:
        SplitFastqByIndexes(input_file, indexFile, index.upper(), param.indexError, param.const_1.upper(), param.const_1Error, param.regExpIndex, args.no_trim_index)
    bcList = (CollectBarcode(indexFile, param.barcodeLength, param.readsValue, param.barcodeError, param.const_2.upper(), param.const_2Error, param.regExpBc, args.merge_indexes, args.reverse_barcode))
    bcCount = Counter(bcList)
    bcCountSorted = sorted(bcCount.items(), key=lambda x:x[1], reverse=True)
    bcLengthMax = max(Counter([len(x[0]) for x in bcCountSorted]).items(), key=lambda x:x[1])
    bcLengthRatio = float(bcLengthMax[1]) / float(len(bcCountSorted))
    if bcLengthRatio >= 0.95: bcCountSorted = [x for x in bcCountSorted if len(x[0]) == bcLengthMax[0]]
    bcCountSorted = [x for x in bcCountSorted if x[1] > param.readsValue ]
    bcOutput = {x:i for x, i in bcCountSorted}
    vennDict[index] = bcOutput

indexCombinations = list(itertools.combinations_with_replacement(param.indexList, 3))
for ic in indexCombinations:
    sets = {}
    for i in ic:
        sets[i] = set([x for x in vennDict[i]])
    iNameTuple = (indexNames[x] for x in ic)
    iNameFig = "venn_" + "_".join([indexNames[x] for x in ic]) + ".png"
    plt.figure(figsize=(10,10))
    plt.title("Venn Diagram for " + " ".join([indexNames[x] for x in ic]))
    venn3([sets[x] for x in sets.keys()], iNameTuple)
    plt.savefig(os.path.join(workdir, iNameFig), fmt='png')

