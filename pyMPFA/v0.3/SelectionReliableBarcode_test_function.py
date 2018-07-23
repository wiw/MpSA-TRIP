#C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8

import math
from collections import OrderedDict, Counter
import progressbar

barcodeError = 2
barcodeLength = 18

def build_hash_word(word):
    nc_hash = {
        "A": 1,
        "T": 2,
        "G": 3,
        "C": 4
    }
    weight_letter_list = [(1 / float(len(word))) * (x + 1) for x in range(len(word))]
    _hash = [nc_hash[word[w]] * weight_letter_list[w] for w in range(len(word))]
    return _hash

# def rough_score(barcodeError, barcodeLength):
#     _tpl_min_score = "A" * int(barcodeLength)
#     _tpl_max_score = str(_tpl_min_score) + "C" * int(barcodeError)
#     min_hash, max_hash = build_hash_word(_tpl_min_score), build_hash_word(_tpl_max_score)
#     _score = float(max_hash - min_hash)
#     return _score

def select_reliable_barcode(bcCount, readsValue, barcodeError):
    _bcDict = {}
    # Sorting
    bcCountSorted = sorted(bcCount.items(), key=lambda x:x[1], reverse=True)
    # Get the most common length of barcodes
    bcLengthMax = max(Counter([len(x[0]) for x in bcCountSorted]).items(), key=lambda x:x[1])
    # Estimate often part  of barcode length
    bcLengthRatio = float(bcLengthMax[1]) / float(len(bcCountSorted))
    # If part more then 99% that barcodes to leave with the same length
    if bcLengthRatio >= 0.99:
        bcCountSorted = [x for x in bcCountSorted if len(x[0]) == bcLengthMax[0]]
    # Select barcodes with reads count more then readsValue
    bcCountSorted = [x for x in bcCountSorted if x[1] > readsValue]
    if barcodeError == 0:
        _bcDict = {bc: [(bc, seq)] for bc, seq in bcCountSorted}
    else:
        # Convert 'bcCountSorted' in dict
        _bc_count = {bc:seq for bc, seq in bcCountSorted}
        # Sorting by value the dictionary '_bc_count'
        _ordered_bc_count = OrderedDict(sorted(_bc_count.items(), key = lambda k: k[1], reverse = True))
        # Get hash for each bc's in dict
        _bc_scoring = {bc: build_hash_word(bc) for bc in _ordered_bc_count}
        # Visual progress bar
        records = len(_ordered_bc_count)
        bar = progressbar.ProgressBar(maxval=records, widgets=[progressbar.SimpleProgress()]).start()
        t=0
        for bc in _ordered_bc_count:
            bar.update(t)
            t += 1
            # Testing existance of bc in hash dictionary
            if bc in _bc_scoring:
                # Calculate hash vector for barcode
                root_bc_hash = _bc_scoring[bc]
                # Obtaining a list of similar barcodes, including the barcode itself
                mutated_bc_list = [mutated_bc for mutated_bc, _hash in _bc_scoring.iteritems() if len(set(root_bc_hash).difference(_hash)) <= barcodeError]
                # Assign count for barcodes and sorting
                mutated_bc_count_list = [(m_bc, _ordered_bc_count[m_bc]) for m_bc in mutated_bc_list]
                mutated_bc_count_list = sorted(mutated_bc_count_list, key = lambda k: k[1], reverse = True)
                # Remove founded barcodes from hash dictionary
                for used_bc in mutated_bc_list:
                    _bc_scoring.pop(used_bc, None)
                # Write result to output dictionary
                _bcDict.setdefault(bc, []).extend(mutated_bc_count_list)
        bar.finish()
    return _bcDict

def main():
    pass

if __name__ == '__main__':
    main()

bbmut, bbumut = {}, {}
for bc, bcmut in aamut.items():
    bbmut.setdefault(bc, []).extend([build_hash_word(x) for x, y in bcmut])
    bc, count = bcmut[0]
    bbumut.setdefault(bc, build_hash_word(bc))

ddmut = {}
for bc, bcmutcount in bbmut.items():
    mainbc = bcmutcount[0]
    minbc = min(bcmutcount)
    maxbc = max(bcmutcount)
    rslt = max([bcmutcount[0]-minbc, maxbc - bcmutcount[0]])
    ddmut.setdefault(bc, rslt)

_bins = math.log(len(ddmut), 2) + 1
n, bins, patches = plt.hist([v for k,v in ddmut.items()], _bins)
plt.show()