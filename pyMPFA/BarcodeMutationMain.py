# C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8

import os
import subprocess
import SupportFunc as supp
import ReadIndexesFunc as rind
import CollectBcMutFunc as colb
import ReliableCombBcMutFunc as relc
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
        supp.log_info(
            "\n\nEnd splitting.\n\n#####################################\n")
        supp.log_info('''Processing on: '{}'.\n
    Total reads in file '{}': {} reads.\n
    Generate dictionary of barcodes.\n'''.format(os.path.basename(indexFile), os.path.basename(indexFile), supp.get_sequence_count(indexFile)))
        (bcDict, seqDict) = colb.CollectBarcodeMutation(indexFile, param.barcodeLength, param.mutationLength, param.readsValue, param.barcodeError,
                                                        param.const_2.upper(), param.const_3.upper(), param.const_2Error, param.const_3Error, param.regExpBcMut, picks.reverse_barcode)
        pickle_dump_bcDict = Pdump(bcDict, name + "_bcDict", picks.PdumpDir)
        pickle_dump_seqDict = Pdump(seqDict, name + "_seqDict", picks.PdumpDir)
        length_of_bcdict = supp.make_histogramm_plot(pickle_dump_bcDict)
        supp.log_info('''        I have get a {} unique barcodes.\n
    Generate  dictionary of <<Barcode vs mutation>>.\n'''.format(length_of_bcdict))
        supp.log_info('''        Unique combinations of barcode-mutation: {} items.\n
    Select the very probable mutations.\n'''.format(len(seqDict)))
        resultDict = relc.SelectTheMostProbableMutation(
            seqDict, param.mutationProbability)
        pickle_dump_resultDict = Pdump(
            resultDict, name + "_resultDict", picks.PdumpDir)
        csvFile = wrt.WriteResultsToFile(
            resultDict, bcDict, seqDict, os.path.join(picks.workdir, name), indexFile)
        csvFile_R = wrt.SimpleCsvWriter(
            resultDict, bcDict, os.path.join(picks.workdir, name), indexFile)
        supp.log_info('''        I had select the {} unique mutations.\nResults writing to file '{}'in your working directory: '{}'\nDump files:\nbarcodes dictionary: {}\nsequences dictionary: {}\nresults dictionary: {}'''.format(
            len(resultDict), csvFile, os.path.join(picks.workdir, name), pickle_dump_bcDict, pickle_dump_seqDict, pickle_dump_resultDict))
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
    supp.log_info("... End script work ...")


if __name__ == "__main__":
    main()
