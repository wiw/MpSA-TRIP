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
from TripMain_0_2 import Pdump


def main():
    supp.setup_logging()
    for name in param.indexList:
        index = param.indexList[name]
        if not os.path.exists(os.path.join(picks.workdir, name)):
            os.makedirs(os.path.join(picks.workdir, name))
        indexFile = os.path.join(
            picks.workdir, name, "index_{}.fastq".format(index.upper()))
        if not os.path.exists(indexFile) or os.stat(indexFile).st_size == 0:
            rind.SplitFastqByIndexes(picks.input_file, indexFile, index.upper(
            ), param.indexError, param.const_1.upper(), param.const_1Error, param.regExpIndex, picks.no_trim_index)
        if picks.random_read:
            indexFileRand = os.path.join(
                picks.workdir, name, "random_index_{}.fastq".format(index.upper()))
            rind.RandomReadIndexes(indexFile, indexFileRand, param.probability)
            indexFile = indexFileRand
        supp.log_info("End splitting.")
        supp.log_info('''Processing on: '{}'.\n
    Total reads in file '{}': {} reads.\n
    Generate dictionary of barcodes.\n'''.format(os.path.basename(indexFile), os.path.basename(indexFile), supp.get_sequence_count(indexFile)))
        bcDict = colb.CollectBarcode(indexFile, param.barcodeLength, param.readsValue, param.barcodeError, param.const_2.upper(
        ), param.const_2Error, param.regExpBc, picks.merge_indexes, picks.reverse_barcode)
        pickle_dump_file = Pdump(bcDict, name + "_bcDict", picks.PdumpDir)
        length_of_bcdict = supp.make_histogramm_plot(pickle_dump_file)
        if len(bcDict) <= 10**5:
            supp.log_info("        Checking barcodes ... Estimated time: ~ {}".format(
                supp.estimate_calculation_time(bcDict)))
            chkb.mainCheckBarcodeInDict(bcDict, param.barcodeError)
        csvFile = wrt.WriteBcDictToFile(
            bcDict, os.path.join(picks.workdir, name), indexFile)
        csvFile_R = wrt.SimpleCsvWriter(
            None, bcDict, os.path.join(picks.workdir, name), indexFile)
        supp.log_info('''        I had select the {} unique barcodes.\nResults writing to file '{}'in your working directory: '{}'\nDump files:\nbarcodes dictionary: {}\n'''.format(
            length_of_bcdict, csvFile, os.path.join(picks.workdir, name), pickle_dump_file))
        if os.path.exists(param.rscript):
            pathToScript = os.path.join(os.getcwd(), "trip_Rstat.R")
            option = [csvFile_R, os.path.dirname(csvFile_R), index]
            cmd = [param.rscript, pathToScript] + option
            subprocess.call(cmd)
        else:
            supp.log_info(
                "You do not have installed R-session, or you incorrectly specified the path to the Rscript.\nStatistics on barcodes will not be displayed.")
        supp.log_info("End processing with: '{}'.\n\n".format(
            os.path.basename(indexFile)))


if __name__ == "__main__":
    main()
