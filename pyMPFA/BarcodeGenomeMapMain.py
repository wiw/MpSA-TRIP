#!/usr/bin/env python
# encoding: utf-8
import os
import param
import picks
import SupportFunc as supp
import ReadIndexesFunc as rind
import CollectBcMutFunc as colb
import ReliableCombBcMutFunc as relc
import WriteFunc as wrt
from TripMain_0_2 import Pdump


def main():
    supp.setup_logging()
    for name in param.indexList:
        supp.log_info("Working in {}".format(name))
        index = param.indexList[name]
        resultDictPI = {}
        if not os.path.exists(os.path.join(picks.workdir, name)):
            os.makedirs(os.path.join(picks.workdir, name))
        indexFile = os.path.join(
            picks.workdir, name, "index_{}.fastq".format(index.upper()))
        if not os.path.exists(indexFile) or os.stat(indexFile).st_size == 0:
            rind.SplitFastqByIndexes(picks.input_file, indexFile, index.upper(
            ), param.indexError, param.const_1.upper(), param.const_1Error, param.regExpIndex, picks.no_trim_index)
            supp.log_info("  End splitting from {}".format(name))
        if picks.random_read:
            indexFileRand = os.path.join(
                picks.workdir, name, "random_index_{}.fastq".format(index.upper()))
            rind.RandomReadIndexes(indexFile, indexFileRand, param.probability)
            indexFile = indexFileRand
        supp.log_info('  Total reads in file {}: {} reads.\nGenerate dictionary of barcodes.\n'.format(
            os.path.basename(indexFile), supp.get_sequence_count(indexFile)))
        bcDictPI, seqDictPI = colb.CollectBarcodeMutationGenome(indexFile, param.barcodeLength, param.readsValue, param.barcodeError, param.const_2.upper(
        ), param.const_2Error, param.regExpBcMut, picks.merge_indexes, param.pmi, param.pmiLength, param.pmiSubst, picks.reverse_barcode, param.bwIndex, param.rfplIndex)
        Pdump(bcDictPI, name + "_bcDictPI", picks.PdumpDir)
        Pdump(seqDictPI, name + "_seqDictPI", picks.PdumpDir)
        for pI in param.pmi:
            resultDictPI[pI] = relc.SelectTheMostProbableMutation(
                seqDictPI[pI], param.mutationProbability)
            csvFile = wrt.WriteResultsToFile(resultDictPI[pI], bcDictPI[pI], seqDictPI[pI], os.path.join(
                picks.workdir, name), indexFile, customTxt=pI)
            supp.log_info('  Working in {}\n         Number of unique barcode: {}.\n         Number of unique combinations of barcode-mutation: {}'.format(
                pI, len(bcDictPI[pI]), len(resultDictPI[pI])))
        Pdump(resultDictPI, name + "_resultDictPI", picks.PdumpDir)
    supp.log_info("... End of program ...")


if __name__ == "__main__":
    main()
