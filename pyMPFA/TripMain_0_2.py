#!/usr/bin/env python
# encoding: utf-8
"""
    trip_0.1.py
    programm software version 0.1
    Search for barcodes and associated mutations.
    Report on the work done and statistics output.\n
    author: Anton V. Ivankin
    e-mail: anton.ivankin.gmail.com
    source: github.com/.....
"""
import os
import pickle
import subprocess
import argparse
import json
import logging
import PairedEndFunc as pend
import SupportFunc as supp
import param

config = {}
Logger = logging.getLogger(__name__)
"""
**** SAVE/LOAD FUNCTIONS
"""


def load_main_config(config_path):
    if os.path.exists(config_path) and os.path.isfile(config_path):
        with open(config_path, "rt") as handle:
            config = json.load(handle)
            return config
    else:
        raise IOError
        Logger.exception("Don't load config file from '{}'".format(config_path))


def SaveDictToPy(dictVar, filename):
    with open(filename + ".py", "wb") as handle:
        for k, v in dictVar.items():
            if type(v) == str:
                handle.write(str(k) + " = '" + str(v) + "'\n")
            else:
                handle.write(str(k) + " = " + str(v) + "\n")


def Pdump(obj, name, folder):
    locationObj = os.path.join(folder, name)
    filename = locationObj + ".pickle"
    with open(filename, 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return filename


def Pload(name, folder):
    locationObj = os.path.join(folder, name)
    with open(locationObj + '.pickle', 'rb') as handle:
        unserialized_data = pickle.load(handle)
        return unserialized_data


MODULES = ["PairedEndFunc", "SupportFunc", "ReadIndexesFunc", "CollectBcMutFunc",
           "MakeRandSeqFunc", "CheckBarcodesFunc", "ReliableCombBcMutFunc", "WriteFunc",
           "BarcodeMutationMain", "BarcodeLibraryMain", "BarcodeGenomeLibMain", "BarcodeGenomeMapMain"]


"""
**** PARSE ARGUMENTS
"""


def ParseArguments():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", "-i", help="fastq file for TRIP experiment")
    p.add_argument("--forward", "-f", help="mapping forward fastq")
    p.add_argument("--reverse", "-r", help="mapping reverse fastq")
    p.add_argument("--mode", "-m", choices=["single", "paired", "fwd", "rev", "genome"],
                   default="single", help="prints <mode> to stdout where \
                        R1 is equivalent to reads without a pair in R2 [%(default)s]")
    p.add_argument("--output", "-o", help="Set location of output folder. Use absolute path!")
    p.add_argument("--no_trim_index", "-T", action="store_true", help="Don't trim indexes with first concatenated \
                                                                constant part. Trim by default. \
                                                                If you use this option, be sure to check \
                                                                the regular expressions in your parameter file - option '-pf'.")
    p.add_argument("--reverse_complement", "-rc", action="store_false",
                   default=True, help="reverse complement R2 when joining [%(default)s]. Default value is [%(default)s]")
    p.add_argument("--gap_size", "-G", default=0, help="If you have gap between forward and reverse reads,\
                                     then please specify gap size in b.p. Default value is [%(default)s]")
    p.add_argument("--experiment_label", "-l", help="A label that is added to the name of the working folder \
                                                        to separate your experiments.")
    p.add_argument("--random_read", "-R", action="store_true", help="Random read by sequence file.")
    p.add_argument("--make_barcode_library", "-B", action="store_true", help="Make only barcode library.")
    p.add_argument("--merge_indexes", "-M", action="store_true", help="Merge data into one table for all available indexes")
    p.add_argument("--reverse_barcode", "-rb", action="store_true", help="Reverse your barcode if you use reversed reads or barcodes")
    p.add_argument("--combine_paired", "-C", action="store_true", help="Merge paired end sequences")
    args = p.parse_args()
    return args

"""
**** CHECK ARGUMENTS
"""


def CheckModules(MODULES):
    report = {}
    for module in MODULES:
        modulePy = module + ".py"
        if os.path.exists(modulePy) and os.path.isfile(modulePy):
            report[module] = True
        else:
            report[module] = False
    if all(report):
        result = True
        wrongImport = ""
    else:
        result = False
        wrongImport = ", ".join([i + ".py" for i in report if report[i] == False])
    return (result, wrongImport)


def CheckArgs(args):
    variableSet = {}
    if args.mode == "single" or args.mode == "genome":
        try:
            if os.path.isabs(args.input):
                variableSet["input_file"] = args.input
                if args.output:
                    if os.path.isabs(args.output):
                        variableSet["inputLocation"] = args.output
                    else:
                        supp.log_error("\nPlease specify the absolute path to the OUTPUT folder\n")
                else:
                    variableSet["inputLocation"] = os.path.dirname(variableSet["input_file"])
                    supp.log_error("\nYou are not specify path to OUTPUT folder!\n\nThe result of the work will be recorded in the folder with the source file!\n")
            else:
                supp.log_error("\nPlease specify the absolute path to the fastq file!\n")
        except:
            supp.log_error("\nYou are in {} mode!\n\nPlease choose ONE fastq file with argument \"--input (-i) ...\"!\n\nExit programm...\n".format(args.mode))
            raise EOFError
    else:
        try:
            if os.path.isabs(args.forward):
                variableSet["r1"] = args.forward
            else:
                supp.log_error("\nPlease specify the absolute path to the forward fastq file!\n")
            if os.path.isabs(args.forward):
                variableSet["r2"] = args.reverse
            else:
                supp.log_error("\nPlease specify the absolute path to the reversed fastq file!\n")
            if args.output:
                if os.path.isabs(args.output):
                    variableSet["inputLocation"] = args.output
                else:
                    supp.log_error("\nPlease specify the absolute path to the OUTPUT folder\n")
            else:
                variableSet["inputLocation"] = os.path.dirname(args.forward)
                supp.log_error("\nYou are not specify path to OUTPUT folder!\n\nThe result of the work will be recorded in the folder with the source file!\n")
        except:
            if args.mode == "paired":
                modeKeyword = args.mode
            elif args.mode == "fwd":
                modeKeyword = "single forward"
            else:
                modeKeyword = "single reversed"
            supp.log_error("\nYou are in \"{}\" mode!\n\n\
                Please choose TWO fastq file with arguments \"--forward (-f) ... --reverse (-r) ...\"!\n\n\
                Exit programm...\n".format(modeKeyword))
            raise EOFError
        if args.combine_paired:
            supp.log_info('''\n            All checks are completed!
            Start working with paired reads...\n\n''')
            variableSet["input_file"] = pend.FastqJoinPaired(variableSet["r1"],
                                                             variableSet["r2"],
                                                             output_dir=variableSet["inputLocation"],
                                                             gap_size=args.gap_size,
                                                             separator=param.separator,
                                                             mode=args.mode,
                                                             reverse_complement=args.reverse_complement)
            supp.log_info('''\n End of combine reads...\n\n''')
    # make fake input_file for paired mode
    if args.mode == "paired":
        if variableSet.get("input_file") is None and variableSet.get("r1") is not None:
            variableSet["input_file"] = variableSet["r1"]
    supp.log_info('''\n            START MAIN PROGRAMM!\n
    author: Anton V. Ivankin\n
    e-mail: anton.ivankin.gmail.com\n
    source: github.com/.....\n\n
    #####################################\n
    Total reads count in your file: {} reads.\n
    Start splitting source file by index.\n\n'''.format(supp.get_sequence_count(variableSet["input_file"])))
    if args.experiment_label:
        variableSet["workdir"] = os.path.join(variableSet["inputLocation"], os.path.basename(variableSet["input_file"]).split(".")[0]+"_{}".format(args.experiment_label))
    else:
        variableSet["workdir"] = os.path.join(variableSet["inputLocation"], os.path.basename(variableSet["input_file"]).split(".")[0])
    supp.log_info("You are working directory in {}\n".format(variableSet["workdir"]))
    variableSet["PdumpDir"] = os.path.join(variableSet["workdir"], "Dump")
    if not os.path.exists(variableSet["PdumpDir"]):
        os.makedirs(variableSet["PdumpDir"])
    return variableSet

"""
**** RUN FUNCTION
"""


def main(args):
    supp.setup_logging()
    (result, wrongImport) = CheckModules(MODULES)
    if result:
        import param
        supp.log_info("        Parameters loaded successfully!\n\
            indexes: {}\n\
            index error: {}\n\
            barcode error: {}\n\
            barcode length: {}\n\
            mutation length: {}\n\
            reads value: {}\n\
            mutation probability: {}\n\
            separator: \"{}\"\n\
            const_1: {}\n\
            const_2: {}\n\
            const_3: {}\n\
            const_1Error: {}\n\
            const_2Error: {}\n\
            const_3Error: {}\n\
            regExpIndex: {}\n\
            ".format(", ".join([i for i in param.indexList.values()]),
                     param.indexError,
                     param.barcodeError,
                     param.barcodeLength,
                     param.mutationLength,
                     param.readsValue,
                     param.mutationProbability,
                     param.separator,
                     param.const_1.upper(),
                     param.const_2.upper(),
                     param.const_3.upper(),
                     param.const_1Error,
                     param.const_2Error,
                     param.const_3Error,
                     ', '.join([i for i in param.regExpIndex.values()])))
    else:
        supp.log_error("Attention! Not all the necessary files are in the folder with the program. \
            Please place '{}' into directory '{}' with the program.".format(wrongImport, os.getcwd()))
    variableSet = CheckArgs(args)
    variableSet.update(vars(args))
    SaveDictToPy(variableSet, "picks")
    if variableSet["mode"] == "single":
        if args.merge_indexes:
            if args.make_barcode_library:
                # scriptPy = "BarcodeLibraryMergeMain.py"
                supp.log_error("This options is not avalable now...")
                pass
            else:
                # scriptPy = "BarcodeMutationMergeMain.py"
                supp.log_error("This options is not avalable now...")
                pass
        else:
            if args.make_barcode_library:
                scriptPy = "BarcodeLibraryMain.py"
            else:
                scriptPy = "BarcodeMutationMain.py"
    elif variableSet["mode"] == "genome":
        if args.make_barcode_library:
            scriptPy = "BarcodeGenomeLibMain.py"
        else:
            scriptPy = "BarcodeGenomeMapMain.py"
    elif variableSet["mode"] == "paired":
        if args.make_barcode_library:
            # scriptPy = "BarcodePairedLibMain.py"
            supp.log_error("This options is not avalable now...")
            pass
        else:
            scriptPy = "BarcodePairedMapMain.py"
    else:
        # scriptPy = "AnotherCrazyFile.py"
        supp.log_error("This options is not avalable now...")
        pass
    cmd = ["python", scriptPy]
    subprocess.call(cmd)

if __name__ == "__main__":
    args = ParseArguments()
    if args.input or args.forward or args.reverse:
        main(args)
    else:
        print("You are in incorrect mode. Use '-h' then to get help.")

# if args.merge_indexes: bcList = []
# for index in param.indexList:
#     indexFile = os.path.join(workdir, "index_{}.fastq".format(index.upper()))
#     if not os.path.exists(indexFile) or os.stat(indexFile).st_size == 0:
#         rind.SplitFastqByIndexes(input_file, indexFile, index.upper(), param.indexError, param.const_1.upper(), param.const_1Error, param.regExpIndex, args.no_trim_index)
#     if args.random_read:
#         indexFileRand = os.path.join(workdir, "random_index_{}.fastq".format(index.upper()))
#         rind.RandomReadIndexes(indexFile, indexFileRand, param.probability)
#         indexFile = indexFileRand
#     print("\n\nEnd splitting.\n\n#####################################\n")
#     print('''Processing on: '{}'.\n
# Total reads in file '{}': {} reads.\n
# Generate dictionary of barcodes.\n'''.format(os.path.basename(indexFile), os.path.basename(indexFile), get_sequence_count(indexFile)))
#     if args.make_barcode_library:
#         if args.merge_indexes:
#             bcList.extend(CollectBarcode(indexFile, param.barcodeLength, param.readsValue, param.barcodeError, param.const_2.upper(), param.const_2Error, param.regExpBc, args.merge_indexes, args.reverse_barcode))
#         else:
#             bcDict = CollectBarcode(indexFile, param.barcodeLength, param.readsValue, param.barcodeError, param.const_2.upper(), param.const_2Error, param.regExpBc, args.merge_indexes, args.reverse_barcode)
#             # mainCheckBarcodeInDict(bcDict, param.barcodeError)
#             csvFile = WriteBcDictToFile(bcDict, workdir, indexFile)
#             csvFile_R = SimpleCsvWriter(None, bcDict, workdir, indexFile)
#             print('''        I had select the {} unique barcodes.\n
#         Results writing to file '{}'
#         in your working directory: '{}'\n'''.format(len(bcDict), csvFile, workdir))
#     else:
#         (bcDict, seqDict) = CollectBarcodeMutation(indexFile, param.barcodeLength, param.mutationLength, param.readsValue, param.barcodeError, param.const_2.upper(), param.const_3.upper(), param.const_2Error, param.const_3Error, param.regExpBcMut, args.reverse_barcode)
#         print('''        I have get a {} unique barcodes.\n
#     Generate  dictionary of <<Barcode vs mutation>>.\n'''.format(len(bcDict)))
#         print('''        Unique combinations of barcode-mutation: {} items.\n
#     Select the very probable mutations.\n'''.format(len(seqDict)))
#         resultDict = SelectTheMostProbableMutation(seqDict, param.mutationProbability)
#         csvFile = WriteResultsToFile(resultDict, bcDict, seqDict, workdir, indexFile)
#         csvFile_R = SimpleCsvWriter(resultDict, bcDict, workdir, indexFile)
#         print('''        I had select the {} unique mutations.\n
#     Results writing to file '{}'
#     in your working directory: '{}'\n'''.format(len(resultDict), csvFile, workdir))
#     if not args.merge_indexes:
#         if os.path.exists(param.rscript):
#             pathToScript = os.path.join(os.getcwd(), "trip_Rstat_{}.R".format(index))
#             option = [csvFile_R, os.path.dirname(csvFile_R), index]
#             cmd = [param.rscript, pathToScript] + option
#             subprocess.call(cmd)
#         else:
#             print("You do not have installed R-session, or you incorrectly specified the path to the Rscript.\nStatistics on barcodes will not be displayed.")
#     print("End processing with: '{}'.\n\n".format(os.path.basename(indexFile)))
# if args.merge_indexes:
#     bcCount = Counter(bcList)
#     if args.temp_option:
#         pass
#     else:
#         indexMerged = "_".join(param.indexList)
#         bcDict = SelectionReliableBarcode(bcCount, param.readsValue, param.barcodeError)
#         # mainCheckBarcodeInDict(bcDict, param.barcodeError)
#         csvFile = WriteBcDictToFile(bcDict, workdir, indexFile)
#         csvFile_R = SimpleCsvWriter(None, bcDict, workdir, indexFile)
#         print('''        I had select the {} unique barcodes.\n
#         Results writing to file '{}'
#         in your working directory: '{}'\n'''.format(len(bcDict), csvFile, workdir))
#         if os.path.exists(param.rscript):
#             pathToScript = os.path.join(os.getcwd(), "trip_Rstat_{}.R".format(index))
#             option = [csvFile_R, os.path.dirname(csvFile_R), indexMerged]
#             cmd = [param.rscript, pathToScript] + option
#             subprocess.call(cmd)
#         else:
#             print("You do not have installed R-session, or you incorrectly specified the path to the Rscript.\nStatistics on barcodes will not be displayed.")
