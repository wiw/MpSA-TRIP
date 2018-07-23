#C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8

import os
# import subprocess
import SupportFunc as supp
import ReadIndexesFunc as rind
import CollectBcMutFunc as colb
import WriteFunc as wrt
import param
import picks
from TripMain_0_2 import Pdump

def main():
    supp.setup_logging()
    for name in param.indexList:
        index = param.indexList[name]
        if not os.path.exists(os.path.join(picks.workdir, name)): os.makedirs(os.path.join(picks.workdir, name))
        indexFile = os.path.join(picks.workdir, name, "index_{}.fastq".format(index.upper()))
        # indexFiltFile = os.path.join(picks.workdir, name, "filt_index_{}.fastq".format(index.upper()))
        if not os.path.exists(indexFile) or os.stat(indexFile).st_size == 0:
            rind.SplitFastqByIndexes(picks.input_file, indexFile, index.upper(), param.indexError, param.const_1.upper(), param.const_1Error, param.regExpIndex, picks.no_trim_index)
        if picks.random_read:
            indexFileRand = os.path.join(picks.workdir, name, "random_index_{}.fastq".format(index.upper()))
            rind.RandomReadIndexes(indexFile, indexFileRand, param.probability)
            indexFile = indexFileRand
        supp.LogInfo("\n\nEnd splitting.\n\n#####################################\n")
        # readsStat[name] = rind.filterShadyReads(indexFile, param.reFilter, indexFiltFile)
        # indexFile = indexFiltFile
        # supp.LogInfo("Filter before: {}, after: {}\n indexFile - {}, indexFiltFile - {}".format(readsStat[name][0], readsStat[name][1], indexFile, indexFiltFile))
        supp.LogInfo('''Processing on: '{}'.\n
        Total reads in file '{}': {} reads.\n
        Generate dictionary of barcodes.\n'''.format(os.path.basename(indexFile), os.path.basename(indexFile), supp.GetTotalSeqRecords(indexFile)))
        bcDictPI = colb.CollectBarcodeGenome(indexFile, param.barcodeLength, param.readsValue, param.barcodeError, param.const_2.upper(), param.const_2Error, param.regExpBc, picks.merge_indexes, picks.reverse_barcode, param.pmi, param.pmiLength, param.pmiSubst)
        Pdump(bcDictPI, name + "_bcDictPI", picks.PdumpDir)
        # Pdump(readsStat, name + "_readsStat", picks.PdumpDir)
        for pI in bcDictPI:
            csvFile = wrt.WriteBcDictToFile(bcDictPI[pI], os.path.join(picks.workdir, name), indexFile, pI)
            # csvFile_R = wrt.SimpleCsvWriter(None, bcDictPI[pI], os.path.join(picks.workdir, name), indexFile, pI)
            supp.LogInfo('''        I had select the {} unique barcodes.\n
            Results writing to file '{}'
            in your working directory: '{}'\n'''.format(len(bcDictPI[pI]), csvFile, os.path.join(picks.workdir, name)))
            # if os.path.exists(param.rscript):
            #     pathToScript = os.path.join(os.getcwd(), "trip_Rstat.R")
            #     option = [csvFile_R, os.path.dirname(csvFile_R), index]
            #     cmd = [param.rscript, pathToScript] + option
            #     subprocess.call(cmd)
            # else:
            #     print("You do not have installed R-session, or you incorrectly specified the path to the Rscript.\nStatistics on barcodes will not be displayed.")
        supp.LogInfo("End processing with: '{}'.\n\n".format(os.path.basename(indexFile)))

if __name__ == "__main__":
    main()