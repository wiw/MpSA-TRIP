#C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8
import os, csv
"""
**** WRITE RESULTS TO FILE FUNCTIONS
"""

def WriteResultsToFile(resultDict, bcDict, seqDict, workdir, indexFile, customTxt=''):
    indexString = os.path.basename(indexFile).split(".")[0].split("_")[1]
    csvFile = os.path.join(workdir, "{}_barcode-mutation_count_{}.csv".format(indexString, customTxt))
    with open(csvFile, "wb") as handle:
        fieldnames = ['Barcode', 'Mutation', 'MutationCount', 'BCSequence', 'BCSequenceCount', 'MutatedBCassociatedWith', 'MutationVariants', 'Frequency', 'LostBarcodes']
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for barcodeID in resultDict:
            BCSequence = ""
            BCSequenceCount = ""
            BCSequence = '\n'.join(str(x[0]) for x in bcDict[barcodeID]) + "\n"
            BCSequenceCount = '\n'.join(str(x[1]) for x in bcDict[barcodeID]) + "\n"
            TotalBCSequence = sum([x[1] for x in bcDict[barcodeID]])
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
                    tmplist=[]
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
                'LostBarcodes': TotalBCSequence-resultDict[barcodeID][1]})
    return os.path.basename(csvFile)

def WriteBcDictToFile(bcDict, workdir, indexFile, customTxt=''):
    indexString = os.path.basename(indexFile).split(".")[0].split("_")[1]
    csvFile = os.path.join(workdir, "{}_barcodeDictionary_{}.csv".format(indexString, customTxt))
    with open(csvFile, "wb") as handle:
        fieldnames = ['Barcode', 'BCSequence', 'BCSequenceCount']
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for barcodeID in bcDict:
            BCSequence = ""
            BCSequenceCount = ""
            BCSequence = '\n'.join(str(x[0]) for x in bcDict[barcodeID]) + "\n"
            BCSequenceCount = '\n'.join(str(x[1]) for x in bcDict[barcodeID]) + "\n"
            writer.writerow({
                'Barcode': barcodeID,
                'BCSequence': BCSequence,
                'BCSequenceCount': BCSequenceCount})
    return os.path.basename(csvFile)

def SimpleCsvWriter(resultDict, bcDict, workdir, indexFile, customTxt=''):
    indexString = os.path.basename(indexFile).split(".")[0].split("_")[1]
    csvFile = os.path.join(workdir, "{}_for_R_statistics{}.csv".format(indexString, customTxt))
    if resultDict is None:
        resultDict = bcDict
    with open(csvFile, "wb") as handle:
        fieldnames = ['Barcode', 'BCSequenceCount']
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for barcodeID in resultDict:
            writer.writerow({
                'Barcode': barcodeID,
                'BCSequenceCount': sum([x[1] for x in bcDict[barcodeID]])})
    return csvFile