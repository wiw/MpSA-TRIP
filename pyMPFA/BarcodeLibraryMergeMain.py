# C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8

import os
import subprocess
import SupportFunc as supp
import ReadIndexesFunc as rind
import CollectBcMutFunc as colb
import CheckBarcodesFunc as chkb
import WriteFunc as wrt
import param
import picks
from MainTrip_0_2 import Pload, Pdump

workdir = Pload("variableSet")["workdir"]
input_file = Pload("variableSet")["input_file"]
PdumpDir = Pload("variableSet")["PdumpDir"]


'''
NOT FINISHED!
'''

# def main():
#     for index in param.indexList:
#         indexFile = os.path.join(workdir, "index_{}.fastq".format(index.upper()))
#         if not os.path.exists(indexFile) or os.stat(indexFile).st_size == 0:
#             rind.SplitFastqByIndexes(input_file, indexFile, index.upper(), param.indexError, param.const_1.upper(), param.const_1Error, param.regExpIndex, picks.no_trim_index)
#         if picks.random_read:
#             indexFileRand = os.path.join(workdir, "random_index_{}.fastq".format(index.upper()))
#             rind.RandomReadIndexes(indexFile, indexFileRand, param.probability)
#             indexFile = indexFileRand
#         print("\n\nEnd splitting.\n\n#####################################\n")
#         print('''Processing on: '{}'.\n
#     Total reads in file '{}': {} reads.\n
#     Generate dictionary of barcodes.\n'''.format(os.path.basename(indexFile), os.path.basename(indexFile), supp.get_sequence_count(indexFile)))
#         bcDict = colb.CollectBarcode(indexFile, param.barcodeLength, param.readsValue, param.barcodeError, param.const_2.upper(), param.const_2Error, param.regExpBc, picks.merge_indexes, picks.reverse_barcode)
#         Pdump(bcDict, os.path.join(PdumpDir, "bcDict"))
#         if len(bcDict) <= 10**5:
#             print("        Checking barcodes ... This can take a long time!")
#             chkb.mainCheckBarcodeInDict(bcDict, param.barcodeError)
#         csvFile = wrt.WriteBcDictToFile(bcDict, workdir, indexFile)
#         csvFile_R = wrt.SimpleCsvWriter(None, bcDict, workdir, indexFile)
#         print('''        I had select the {} unique barcodes.\n
#     Results writing to file '{}'
#     in your working directory: '{}'\n'''.format(len(bcDict), csvFile, workdir))
#         if os.path.exists(param.rscript):
#             pathToScript = os.path.join(os.getcwd(), "trip_Rstat_{}.R".format(index))
#             option = [csvFile_R, os.path.dirname(csvFile_R), index]
#             cmd = [param.rscript, pathToScript] + option
#             subprocess.call(cmd)
#         else:
#             print("You do not have installed R-session, or you incorrectly specified the path to the Rscript.\nStatistics on barcodes will not be displayed.")
#         print("End processing with: '{}'.\n\n".format(os.path.basename(indexFile)))

if __name__ == "__main__":
    main()
