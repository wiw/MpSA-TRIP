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
from pprint import pprint as view
from glob import glob


def main():
    # Configure environment
    supp.setup_logging()
    options = {attribute: value for attribute, value in param.__dict__.items() if attribute[:2] != "__" and not callable(value)}
    options.update({attribute: value for attribute, value in picks.__dict__.items() if attribute[:2] != "__" and not callable(value)})
    collection_of_output_data = {}
    for name, index in options["indexList"].items():
        collection_of_output_data.setdefault(name, {})
        supp.LogInfo("Working in {}".format(name))
        resultDict = {}
        experiment_dir = os.path.join(options["workdir"], name)
        if not os.path.exists(experiment_dir):
            os.makedirs(experiment_dir)
        indexFile = os.path.join(experiment_dir, "index_{}.fastq".format(index.upper()))
        if not os.path.exists(indexFile) or os.stat(indexFile).st_size == 0:
            paired_indexes = rind.SplitFastqByIndexesPaired(index, indexFile, options)
            supp.LogInfo("  End splitting from {}".format(name))
        else:
            files_list = glob(os.path.join(experiment_dir, "*fastq"))
            if len(files_list) > 2:
                supp.LogErr("Too much files. Exiting...")
                break
            else:
                paired_indexes = dict(zip(["fwd", "rev"], sorted(files_list)))
        supp.LogInfo('  {}  {}\
                Generate dictionary of barcodes.\n\
                '.format(*["Total reads in file {}: {} reads.\n\
                    ".format(os.path.basename(value), supp.GetTotalSeqRecords(value)) for value in paired_indexes.values()]))
        bcDict, seqDict, main_paired_stat = colb.main_paired(paired_indexes, options)
        Pdump(bcDict, name + "_bcDictPI", options["PdumpDir"])
        Pdump(seqDict, name + "_seqDictPI", options["PdumpDir"])
        for pmi_item in options["pmi"]:
            resultDict[pmi_item] = relc.SelectTheMostProbableMutation(seqDict[pmi_item], options["mutationProbability"])
            csvFile = wrt.WriteResultsToFile(resultDict[pmi_item], bcDict[pmi_item], seqDict[pmi_item], os.path.join(options["workdir"], name), indexFile, customTxt=pmi_item)
            supp.LogInfo('\                  Working in {}\n         Number of unique barcode: {}.\n\
                       Number of unique combinations of barcode-mutation: {}\
                       '.format(pmi_item, len(bcDict[pmi_item]), len(resultDict[pmi_item])))
            collection_of_output_data[name].setdefault("csvFile", {})
            collection_of_output_data[name]["csvFile"].setdefault(pmi_item, csvFile)
        Pdump(resultDict, name + "_resultDictPI", options["PdumpDir"])
        for big_data in ["main_paired_stat", "paired_indexes", "bcDict", "seqDict", "resultDict"]:
            collection_of_output_data[name].setdefault(big_data, eval(big_data))
    supp.main_report(collection_of_output_data, options)
    supp.LogInfo("... End of program ...")

if __name__ == "__main__":
    main()
