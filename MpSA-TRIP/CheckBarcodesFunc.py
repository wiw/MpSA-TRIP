#C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8
import Levenshtein, progressbar, itertools
import SupportFunc as supp
"""
**** CHECKING BARCODES FOR UNIQUENESS FUNCTIONS
"""

def workerCheckBarcode(bcDict, barcodeError, bcCombinationList):
    stat = 0
    for comb in bcCombinationList:
        if Levenshtein.distance(comb[0], comb[1]) <= barcodeError:
            for key in comb:
                bcDict.pop(key, None)
            stat += 1
    return stat

def mainCheckBarcodeInDict(bcDict, barcodeError):
    # print("        Check unique barcodes for error...\n")
    bcRanges = [(x, x + 1000) for x in range(0, len(bcDict), 1000)]
    partCombinations = list(itertools.combinations_with_replacement(bcRanges, 2))
    bar = progressbar.ProgressBar(maxval=len(partCombinations), widgets=[progressbar.Bar(left='<', marker='.', right='>')]).start()
    t = 0
    stat = 0
    for combItem in partCombinations:
        if combItem[0] == combItem[1]:
            a, b = combItem[0][0], combItem[0][1]
            _comb_ = list(itertools.combinations(bcDict.keys()[a:b], 2))
        else:
            a, b, c, d = combItem[0][0], combItem[0][1], combItem[1][0], combItem[1][1]
            _comb_ = list(itertools.product(bcDict.keys()[a:b], bcDict.keys()[c:d]))
        stat += workerCheckBarcode(bcDict, barcodeError, _comb_)
        bar.update(t)
        t += 1
    bar.finish()
    if stat == 0:
        supp.LogInfo("        I did not find errors in the analysis. All barcodes are correct.\n")
    else:
        supp.LogInfo("        There were {} errors. Do not worry, they are already deleted.\n".format(stat))

def checkPMI(exp, pmi, pmiSubst):
    chkResult = []
    for i in pmi:
        if Levenshtein.distance(str(exp), str(i)) <= int(pmiSubst):
            chkResult.append(str(i))
    if len(chkResult) == 1:
        return chkResult[0]
    return False
