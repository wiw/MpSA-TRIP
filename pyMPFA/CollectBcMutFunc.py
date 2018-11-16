#!/usr/bin/env python
# encoding: utf-8
import regex
import progressbar
import os
import multiprocessing as mp
import subprocess
import SupportFunc as supp
from toolshed import nopen
from collections import Counter
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from PairedEndFunc import reverseComplement
from CheckBarcodesFunc import mainCheckBarcodeInDict, checkPMI
from ReliableCombBcMutFunc import SelectionReliableBarcode, SelectionReliableBarcodeMutation, SelectionReliableBarcodeGenome
"""
**** MPSA section
"""


def CollectBarcodeMutation(indexFile, barcodeLength, mutationLength, readsValue, barcodeError, const_2, const_3, const_2Error, const_3Error, regExpBcMut, reverseBC):
    # print("Start collect barcodes...\nI'm using next regular exxpression for search: {}".format(regExpBcMut))
    bcList, bcMutList = [], []
    non_matched_reads = os.path.join(os.path.dirname(
        indexFile), "undef_{}".format(os.path.basename(indexFile)))
    records = supp.get_sequence_count(indexFile)
    bar = progressbar.ProgressBar(maxval=records, widgets=[
                                  progressbar.Bar(left='<', marker='.', right='>')]).start()
    t = 0
    expr = regex.compile(regExpBcMut)
    with nopen(indexFile) as handle:
        with open(non_matched_reads, "wb") as undef:
            for seq_record in SeqIO.parse(handle, "fastq"):
                bar.update(t)
                t += 1
                match = expr.match(str(seq_record.seq))
                if match is not None:
                    if int(barcodeLength * 0.9) <= len(match.group("barcode")) <= int(barcodeLength * 1.1) and len(match.group("mutation")) == mutationLength:
                        if "N" not in match.group("barcode").upper() and "N" not in match.group("mutation").upper():
                            if reverseBC:
                                bcList.append(reverseComplement(
                                    match.group("barcode")))
                                bcMutList.append((reverseComplement(match.group(
                                    "barcode")), reverseComplement(match.group("mutation"))))
                            else:
                                bcList.append(match.group("barcode"))
                                bcMutList.append(
                                    (match.group("barcode"), match.group("mutation")))
                else:
                    SeqIO.write(seq_record, undef, "fastq")
    bar.finish()
    bcCount = Counter(bcList)
    bcMutCount = dict(Counter(bcMutList))
    # print("        Total founded barcodes: {} items.\nAnd total combinations barcode-mutation: {} items".format(len(bcCount), len(bcMutCount)))
    bcDict = SelectionReliableBarcode(bcCount, readsValue, barcodeError)
    if len(bcDict) <= 10**5:
        # print("        Checking barcodes ... Estimated time: ~ {}".format(supp.estimate_calculation_time(bcDict)))
        mainCheckBarcodeInDict(bcDict, barcodeError)
    seqDict = SelectionReliableBarcodeMutation(bcMutCount, bcDict)
    return (bcDict, seqDict)


def CollectBarcode(indexFile, barcodeLength, readsValue, barcodeError, const_2, const_2Error, regExpBc, mergeBC, reverseBC):
    bcList = []
    records = supp.get_sequence_count(indexFile)
    bar = progressbar.ProgressBar(maxval=records, widgets=[
                                  progressbar.Bar(left='<', marker='.', right='>')]).start()
    t = 0
    expr = regex.compile(regExpBc)
    with nopen(indexFile) as handle:
        for seq_record in SeqIO.parse(handle, "fastq"):
            bar.update(t)
            t += 1
            match = expr.match(str(seq_record.seq))
            if match is not None:
                if int(barcodeLength * 0.9) <= len(match.group("barcode")) <= int(barcodeLength * 1.1):
                    if "N" not in match.group("barcode").upper():
                        if reverseBC:
                            bcList.append(reverseComplement(
                                match.group("barcode")))
                        else:
                            bcList.append(match.group("barcode"))
    bar.finish()
    if mergeBC:
        return bcList
    bcCount = Counter(bcList)
    bcDict = SelectionReliableBarcode(bcCount, readsValue, barcodeError)
    if len(bcDict) <= 10**5:
        # print("        Checking barcodes ... Estimated time: ~ {}".format(supp.estimate_calculation_time(bcDict)))
        mainCheckBarcodeInDict(bcDict, barcodeError)
    return bcDict


"""
**** TRIP section
"""


def CollectBarcodeMutationGenome(indexFile,
                                 barcodeLength,
                                 readsValue,
                                 barcodeError,
                                 const_2,
                                 const_2Error,
                                 regExpBcMut,
                                 merge_indexes,
                                 pmi,
                                 pmiLength,
                                 pmiSubst,
                                 reverseBC,
                                 bwIndex,
                                 rfplIndex):
    bcListDict, bcGenomeListDict = {p: [] for p in pmi}, {p: {} for p in pmi}
    bcDictPI, seqDictPI = {}, {}
    bcDictTmp, alignDict = {}, {}
    expr = regex.compile(regExpBcMut)
    for pmiItem in pmi:
        pmiWD = os.path.join(os.path.dirname(indexFile), pmiItem)
        pmiFile = os.path.join(pmiWD, pmiItem + ".fastq")
        countAllRds, countConst3, countLess20 = 0, 0, 0
        supp.log_info("  Collect from promotor index {}".format(pmiItem))
        if not os.path.exists(pmiWD):
            os.makedirs(pmiWD)
        with nopen(indexFile) as handle, open(pmiFile, "w") as pmiHandle:
            for title, seq, qual in FastqGeneralIterator(handle):
                match = expr.match(seq)
                if match is not None:
                    countAllRds += 1
                    if int(barcodeLength * 0.9) <= len(match.group("barcode")) <= int(barcodeLength * 1.1) and len(match.group("pIndex")) == pmiLength:
                        if "N" not in match.group("barcode").upper() and "N" not in match.group("pIndex").upper():
                            if match.group("genome1") is not None:
                                g = match.group("genome1")
                            else:
                                g = match.group("genome2")
                                countConst3 += 1
                            if reverseBC:
                                b = reverseComplement(match.group("barcode"))
                                p = reverseComplement(match.group("pIndex"))
                            else:
                                b = match.group("barcode")
                                p = match.group("pIndex")
                            matchedPMI = checkPMI(p, pmi, pmiSubst)
                            if matchedPMI in pmi and matchedPMI == pmiItem:
                                if len(g) >= 20:
                                    bcListDict[matchedPMI].append(b)
                                    bcGenomeListDict[matchedPMI][title.split(" ")[
                                        0]] = b
                                    FstGenomeCoord = len(b) + 16
                                    SndGenomeCoord = FstGenomeCoord + len(g)
                                    seqStr = [title, seq[FstGenomeCoord:SndGenomeCoord], qual[
                                        FstGenomeCoord:SndGenomeCoord]]
                                    pmiHandle.write(
                                        "@{}\n{}\n+\n{}\n".format(*seqStr))
                                else:
                                    countLess20 += 1
        supp.log_info("      All reads: {} // Reads with GENOME length less 20 bp: {} // Reads with CONST_3 (20bp): {}".format(
            countAllRds, countLess20, countConst3))
        bcCount = Counter(bcListDict[pmiItem])
        bcDictTmp[pmiItem] = SelectionReliableBarcode(
            bcCount, readsValue, barcodeError)
        if len(bcDictTmp[pmiItem]) <= 10**5 and len(bcDictTmp[pmiItem]) > 0:
            mainCheckBarcodeInDict(bcDictTmp[pmiItem], barcodeError)
        alignDict[pmiItem], bwtAlignerDict = AlignReadsToGenome(
            pmiFile, pmiWD, pmiItem, bwIndex, rfplIndex)
        supp.make_stat_from_bowtie(bwtAlignerDict)
        mergeList = [(bcGenomeListDict[pmiItem][k], alignDict[pmiItem][k])
                     for k in bcGenomeListDict[pmiItem] if k in alignDict[pmiItem]]
        mergeListCount = dict(Counter(mergeList))
        bcDictPI[pmiItem], seqDictPI[pmiItem] = SelectionReliableBarcodeGenome(
            mergeListCount, bcDictTmp[pmiItem])
    return bcDictPI, seqDictPI


def AlignReadsToGenome(pmiFile, pmiWD, pmiItem, bwIndex, rfplIndex):
    cpus = mp.cpu_count()
    gDict, tmpbwt, bwtAlignerDict = {}, {}, {}
    filteredFQ, filterSam = os.path.join(pmiWD, "filt_" + os.path.basename(
        pmiFile)), os.path.join(pmiWD, os.path.basename(pmiFile).split('.')[0] + "_filter.sam")
    bamOut, bedOut = os.path.join(
        pmiWD, "dmel.bam"), os.path.join(pmiWD, "dmel_uniq.bed")
    bwtCore = "(bowtie2 -k 3 -p " + str(cpus) + " -t --phred33 --local -x "
    filterCmd = "--un " + \
        " ".join([filteredFQ, rfplIndex, pmiFile]) + " -S " + filterSam + ")"
    alignCmd = " ".join([bwIndex, filteredFQ]) + \
        " | samtools view -bS - -o " + bamOut + ")"
    samtoolsFilter = "(samtools sort -@ " + str(cpus) + " -o " + os.path.join(pmiWD, "dmel.bam") + \
        " dmel-sort | samtools view -h -F 4 - | samtools view -S -q 25 -F 256 - | grep -Pv 'XS:i' | sam2bed > " + \
        bedOut + ") 2>/dev/null"
    runFilter = subprocess.Popen(bwtCore + filterCmd, shell=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (tmpbwt[pmiItem + '_outF'], tmpbwt[pmiItem + '_errF']
     ) = runFilter.communicate()
    runAlign = subprocess.Popen(bwtCore + alignCmd, shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (tmpbwt[pmiItem + '_outAlg'],
     tmpbwt[pmiItem + '_errAlg']) = runAlign.communicate()
    for k, v in tmpbwt.items():
        if v is not None:
            bwtAlignerDict[k] = v.splitlines()
    subprocess.call(samtoolsFilter, shell=True)
    with open(bedOut) as dmel:
        for line in dmel:
            chrName, start, end, title, strand = [
                line.split()[x] for x in [0, 1, 2, 3, 5]]
            if strand == "-":
                coord = end
            else:
                coord = start
            gDict[title] = (chrName, coord, strand)
    return gDict, bwtAlignerDict


def CollectBarcodeGenome(indexFile, barcodeLength, readsValue, barcodeError, const_2, const_2Error, regExpBc, mergeBC, reverseBC, pmi, pmiLength, pmiSubst):
    bcListDict = {p: [] for p in pmi}
    bcDictPI = {}
    records = supp.get_sequence_count(indexFile)
    bar = progressbar.ProgressBar(maxval=records, widgets=[
                                  progressbar.Bar(left='<', marker='.', right='>')]).start()
    t = 0
    expr = regex.compile(regExpBc)
    with nopen(indexFile) as handle:
        for seq_record in SeqIO.parse(handle, "fastq"):
            bar.update(t)
            t += 1
            match = expr.match(str(seq_record.seq))
            if match is not None:
                if int(barcodeLength * 0.9) <= len(match.group("barcode")) <= int(barcodeLength * 1.1) and len(match.group("pIndex")) == pmiLength:
                    if "N" not in match.group("barcode").upper() and "N" not in match.group("pIndex").upper():
                        if reverseBC:
                            b = reverseComplement(match.group("barcode"))
                            p = reverseComplement(match.group("pIndex"))
                        else:
                            b = match.group("barcode")
                            p = match.group("pIndex")
                        matchedPMI = checkPMI(p, pmi, pmiSubst)
                        if matchedPMI in bcListDict:
                            bcListDict[matchedPMI].append(b)
    bar.finish()
    if mergeBC:
        return bcListDict
    for pI in bcListDict:
        bcCount = Counter(bcListDict[pI])
        bcDictPI[pI] = SelectionReliableBarcode(
            bcCount, readsValue, barcodeError)
        if len(bcDictPI[pI]) <= 10**5 and len(bcDictPI[pI]) > 0:
            # print("        Checking barcodes for promotor index {} ... Estimated time: ~ {}".format(pI, supp.estimate_calculation_time(bcDictPI[pI])))
            mainCheckBarcodeInDict(bcDictPI[pI], barcodeError)
    return bcDictPI


"""
**** paired TRIP section
"""


def CollectBarcodeMutationGenomePaired(paired_indexes, options):
    collected_data, bowtie_args = {}, {}
    where_is_title, four_letters_seq_collection = {}, {}
    for side in ["fwd", "rev"]:
        supp.log_info("  Use {} reads".format(side))
        side_path = paired_indexes[side]
        collected_data.setdefault(side, {})
        bowtie_args.setdefault(side, {})
        expr = regex.compile(options["regExpBcMut_" + side])
        for pmi_item in options["pmi"]:
            supp.log_info("  Collect from promotor index {}".format(pmi_item))
            pmi_wd = os.path.join(os.path.dirname(side_path), pmi_item)
            pmi_file = os.path.join(pmi_wd, side + "_" + pmi_item + ".fastq")
            bowtie_args[side].setdefault(pmi_item, pmi_file)
            if not os.path.exists(pmi_wd):
                os.makedirs(pmi_wd)
            with nopen(side_path) as handle, open(pmi_file, "w") as pmi_handle:
                for title, seq, qual in FastqGeneralIterator(handle):
                    match = expr.match(seq)
                    _title = title.split()[0]
                    if match is not None:
                        match_result = checkup_regex(
                            match.groupdict(), options)
                        if match_result is not False:
                            if side == "fwd":
                                four_letters_seq_collection.setdefault(
                                    side, {})
                                _pmi = match_result.get("pIndex", False)
                                matched_pmi = checkPMI(_pmi, pmi_item, options)
                                if matched_pmi is not False:
                                    # from this stage get 4-letters seq, and collect to
                                    # dict
                                    four_letters_seq = [
                                        x for x in match.groups() if len(x) == 4]
                                    four_letters_seq = four_letters_seq[
                                        0] if len(four_letters_seq) == 1 else ""
                                    if matched_pmi == pmi_item:
                                        four_letters_seq_collection[side].setdefault(
                                            pmi_item, []).append(four_letters_seq)
                                    # end
                                    collected_data[side].setdefault(
                                        matched_pmi, {})
                                    collected_data[side][
                                        matched_pmi].setdefault(_title, {})
                                    where_is_title.setdefault(
                                        _title, matched_pmi)
                                else:
                                    continue
                            else:
                                if _title in where_is_title:
                                    matched_pmi = where_is_title[_title]
                                    collected_data[side].setdefault(
                                        matched_pmi, {})
                                    collected_data[side][
                                        matched_pmi].setdefault(_title, {})
                                else:
                                    continue
                            available_match = [
                                k for k in match_result if k != "pIndex"]
                            if matched_pmi == pmi_item:
                                for _mk in available_match:
                                    if match_result.get(_mk) is not None:
                                        collected_data[side][matched_pmi][
                                            _title].setdefault(_mk, match_result[_mk])
                                    if _mk == 'genome':
                                        _1st, _2nd = match.span(_mk)
                                        seq_str = [title, match_result[
                                            _mk], qual[_1st:_2nd]]
                                        pmi_handle.write(
                                            "@{}\n{}\n+\n{}\n".format(*seq_str))
                        else:
                            continue
    return collected_data, bowtie_args, four_letters_seq_collection


def checkup_regex(match, options):
    match_result = {}
    bc_len, pmi_len, mut_len = options["barcodeLength"], options[
        "pmiLength"], options["mutationLength"]
    for item, value in match.items():
        checkup = []
        if item == "barcode":
            range_item = range(int(int(bc_len) * 0.9),
                               1 + int(int(bc_len) * 1.1))
            checkup.append(True if len(value) in range_item else False)
            checkup.append(True if "N" not in value.upper() else False)
        elif item == "pIndex":
            checkup.append(True if len(value) == pmi_len else False)
            checkup.append(True if "N" not in value.upper() else False)
        elif item == "genome":
            checkup.append(True if type(value) == str else False)
            checkup.append(True if len(value) >= 20 else False)
        elif item == "mutation":
            checkup.append(True if type(value) == str else False)
            checkup.append(True if len(value) == mut_len else False)
        if all(checkup):
            if item in ["barcode", "pIndex"]:
                match_result[item] = reverseComplement(value) if options[
                    "reverse_barcode"] else value
            elif item in ["genome", "mutation"]:
                match_result[item] = value
    if len(match_result) != len(match):
        return False
    return match_result


def comm_fastq(bowtie_args, options):
    comm_stat = {}
    for pmi_item in options["pmi"]:
        _sides = [bowtie_args[side][pmi_item] for side in ["fwd", "rev"]]
        titles = {}
        comm_stat.setdefault(pmi_item, {})
        for _item in _sides:
            with open(_item, "r") as handle:
                for t, s, q in FastqGeneralIterator(handle):
                    item_label = os.path.basename(_item)
                    comm_stat[pmi_item][item_label] = comm_stat[
                        pmi_item].setdefault(item_label, 0) + 1
                    titles.setdefault(_item, []).append(t.split()[0])
        titles_sorted = sorted(titles.iteritems())
        fst, snd = [files[1] for files in titles_sorted]
        common_list = set(fst) & set(snd)
        comm_stat[pmi_item].setdefault("comm", len(common_list))
        for _item in _sides:
            _tmp_file = os.path.join(os.path.dirname(_item), "tmp.fastq")
            with open(_item, "r") as handle, open(_tmp_file, "w") as t_handle:
                for t, s, q in FastqGeneralIterator(handle):
                    if t.split()[0] in common_list:
                        t_handle.write("@{}\n{}\n+\n{}\n".format(t, s, q))
            os.rename(_tmp_file, _item)
        supp.log_info("Complete {} 'comm' files".format(pmi_item))
    return comm_stat


def AlignReadsToGenomePaired(bowtie_args, options):
    cpus = mp.cpu_count()
    max_dist = 5000
    genome_data = {}
    tmp_stat = {}
    bwt_aligner_stat = {}
    for pmi_item in options["pmi"]:
        supp.log_info("Start align paired-end {} file".format(pmi_item))
        fwd, rev = [bowtie_args[side][pmi_item] for side in ["fwd", "rev"]]
        wd = os.path.dirname(fwd)
        genome_data.setdefault(pmi_item, {})
        bam_out, bed_out = os.path.join(
            wd, "dmel.bam"), os.path.join(wd, "dmel_uniq.bed")
        cmd_0 = "bowtie2 -k 3 -p {} -t --phred33 --local -X {} -x {} -1 {} -2 {} | \
            samtools view -bS - -o {}".format(str(cpus),
                                              str(max_dist),
                                              options["bwIndex"],
                                              fwd,
                                              rev,
                                              bam_out)
        run_align = subprocess.Popen(
            cmd_0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (tmp_stat[pmi_item + '_outF'],
         tmp_stat[pmi_item + '_errF']) = run_align.communicate()
        cmd_1 = "(samtools sort -@ {} -o {} dmel-sort | \
            samtools view -h - | \
            samtools view -S -q 25 - | \
            grep -Pv 'XS:i' | \
            grep -P 'AS:i.*YS:i.*YT:Z:CP' | \
            sam2bed > {}) 2>/dev/null".format(str(cpus), bam_out, bed_out)
        subprocess.call(cmd_1, shell=True)
        for k, v in tmp_stat.items():
            if v is not None:
                bwt_aligner_stat[k] = v.splitlines()
        with open(bed_out, "r") as bed_handle:
            header = ['chr', 'start', 'end', 'title',
                      'mapq', 'strand', 'flag', 'length']
            for line in bed_handle:
                bar = dict(zip(header, line.split()[:8]))
                bar["length"] = bar["length"][:-1]
                if get_length(bar["length"]) >= 50:
                    coord = bar["start"] if bar["strand"] == "+" else bar["end"]
                    genome_data[pmi_item].setdefault(
                        bar["title"], (bar["chr"], coord, bar["strand"]))
    return genome_data, bwt_aligner_stat

# The function converts a string with letters and numbers into a list with numbers, which are then added together.
# The function is resistant to entering spaces and control characters.
# Example:
# >>> s_0 = "5S6M30D"
# >>> print(get_length(s_0))
# 41
# >>> s_1 = "5Sqwdcc';6M 30D**/-"
# >>> print(get_length(s_1))
# 41


def get_length(s):
    length = 0
    s = str(s)
    splitted_string = regex.split("(\\d+)", s)
    for item in splitted_string:
        if item.isdigit():
            length += int(item)
    return length


def intersect_collected_and_genome(collected_data, genome_data, options):
    bcDict, seqDict = {}, {}
    for pmi_item in options["pmi"]:
        collected_titles = [k for k in collected_data["fwd"][pmi_item]]
        genome_titles = [k for k in genome_data[pmi_item]]
        intersected_titles = set(collected_titles) & set(genome_titles)
        bc_list = [collected_data["fwd"][pmi_item][tt]["barcode"]
                   for tt in intersected_titles]
        bc_genome_list = [(collected_data["fwd"][pmi_item][tt]["barcode"], genome_data[
                           pmi_item][tt]) for tt in intersected_titles]
        bc_count = Counter(bc_list)
        bcDict.setdefault(pmi_item, SelectionReliableBarcode(
            bc_count, options["readsValue"], options["barcodeError"]))
        bc_genome_count = dict(Counter(bc_genome_list))
        seqDict.setdefault(pmi_item, SelectionReliableBarcodeMutation(
            bc_genome_count, bcDict[pmi_item]))
    return bcDict, seqDict


def main_paired(paired_indexes, options):
    main_paired_stat = {}
    collected_data, bowtie_args, main_paired_stat[
        "four_letters_seq_collection"] = CollectBarcodeMutationGenomePaired(paired_indexes, options)
    main_paired_stat.setdefault("comm_stat", comm_fastq(bowtie_args, options))
    genome_data, main_paired_stat[
        "bwt_aligner_stat"] = AlignReadsToGenomePaired(bowtie_args, options)
    bcDict, seqDict = intersect_collected_and_genome(
        collected_data, genome_data, options)
    return bcDict, seqDict, main_paired_stat
