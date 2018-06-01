#C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8

import os
from collections import Counter
import SupportFunc as supp
import ReadIndexesFunc as rind
import CollectBcMutFunc as colb
import ReliableCombBcMutFunc as relc
import CheckBarcodesFunc as chkb
import WriteFunc as wrt
import param
import picks
from MainTrip_0_2 import Pload, Pdump

workdir = Pload("variableSet")["workdir"]
input_file = Pload("variableSet")["input_file"]
PdumpDir = Pload("variableSet")["PdumpDir"]

def main():
    bcList = []
    for index in param.indexList:
        indexFile = os.path.join(workdir, "index_{}.fastq".format(index.upper()))
        if not os.path.exists(indexFile) or os.stat(indexFile).st_size == 0:
            rind.SplitFastqByIndexes(input_file, indexFile, index.upper(), param.indexError, param.const_1.upper(), param.const_1Error, param.regExpIndex, picks.no_trim_index)
        if picks.random_read:
            indexFileRand = os.path.join(workdir, "random_index_{}.fastq".format(index.upper()))
            rind.RandomReadIndexes(indexFile, indexFileRand, param.probability)
            indexFile = indexFileRand
        print("\n\nEnd splitting.\n\n#####################################\n")
        print('''Processing on: '{}'.\n
    Total reads in file '{}': {} reads.\n
    Generate dictionary of barcodes.\n'''.format(os.path.basename(indexFile), os.path.basename(indexFile), supp.GetTotalSeqRecords(indexFile)))
        bcList.extend(colb.CollectBarcode(indexFile, param.barcodeLength, param.readsValue, param.barcodeError, param.const_2.upper(), param.const_2Error, param.regExpBc, picks.merge_indexes, picks.reverse_barcode))
    '''
    **** NOT FINISHED!
    '''

    # bcCount = Counter(bcList)
    # indexMerged = "_".join(param.indexList)
    # bcDict = relc.SelectionReliableBarcode(bcCount, param.readsValue, param.barcodeError)
    # if len(bcDict) <= 10**5: 
    #     print("        Checking barcodes ... This can take a long time!")
    #     chkb.mainCheckBarcodeInDict(bcDict, param.barcodeError)
    # seqDict = relc.SelectionReliableBarcodeMutation(bcMutCount, bcDict)
    print('''        I have get a {} unique barcodes.\n
    Generate  dictionary of <<Barcode vs mutation>>.\n'''.format(len(bcDict)))
    print('''        Unique combinations of barcode-mutation: {} items.\n
    Select the very probable mutations.\n'''.format(len(seqDict)))
    resultDict = relc.SelectTheMostProbableMutation(seqDict, param.mutationProbability)
    csvFile = wrt.WriteResultsToFile(resultDict, bcDict, seqDict, workdir, indexFile)
    print('''        I had select the {} unique mutations.\n
    Results writing to file '{}'
    in your working directory: '{}'\n'''.format(len(resultDict), csvFile, workdir))
    print("End processing with merged indexes: '{}'.\n\n".format(os.path.basename(indexMerged)))


if __name__ == "__main__":
    main()