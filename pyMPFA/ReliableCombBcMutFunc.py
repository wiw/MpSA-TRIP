#C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8
from collections import Counter
import progressbar, Levenshtein
# import SupportFunc as supp
"""
**** CHOICE OF RELIABLE COMBINATIONS BARCODE/MUTATION FUNCTIONS
"""

def SelectionReliableBarcode(bcCount, readsValue, barcodeError):
    bcDict = {}
    bcHashDict = {}
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
        bcDict = {bc:seq for bc, seq in bcCountSorted}
    else:
        bcSortedList = [x[0] for x in bcCountSorted]
        records = len(bcCountSorted)
        bar = progressbar.ProgressBar(maxval=records, widgets=[progressbar.SimpleProgress()]).start()
        t=0.0
        for bc in bcCountSorted:
            bar.update(t)
            t+=1
            chunks, chunk_size = len(bc[0]), 4
            bcKeyListOne = [bc[0][i:i+chunk_size] for i in range(0, chunks, chunk_size-1)]
            bcKeyListTwo = [i for i in bcKeyListOne if i in bcHashDict.keys()]
            if any(bcKeyListTwo):
                # True
                tmpListBc = []
                for hashKey in bcKeyListTwo: tmpListBc.extend([bcMut for bcMut in bcHashDict[hashKey] if Levenshtein.distance(bc[0], bcMut) <= barcodeError])
                tmpListBc = list(set(tmpListBc))
                if any(tmpListBc):
                    # True
                    if len(tmpListBc) > 1:
                        # We get the frequency of occurrence of the barcode from the dictionary, for those barcodes that are in tmpListBc
                        # Get the maximum frequency for bcFreq
                        # Find the bar code with the maximum frequency of occurrence
                        bcUniqMut = [bcc[0] for bcc in [bcu for bcu in bcCountSorted if bcu[0] in tmpListBc] if bcc[1] == max([m[1] for m in [bcu for bcu in bcCountSorted if bcu[0] in tmpListBc]])]
                        # If the maximum frequency is inherent in two or more barcodes, then we take the bar code with the minimum index
                        if len(bcUniqMut) > 1:
                            bcUniqMut = bcCountSorted[min([bcSortedList.index(x) for x in tmpListBc])][0]
                        else:
                            bcUniqMut = "".join(bcUniqMut)
                        # Write a mutant barcode to a string with a unique bar code
                        bcDict[bcUniqMut].append(bc)
                    else:
                        bcDict[tmpListBc[0]].append(bc)
                else:
                    # False
                    difHash = [y for y in bcKeyListOne if y not in bcKeyListTwo]
                    sameHash = [y for y in bcKeyListOne if y in bcKeyListTwo]
                    if any(difHash): bcHashDict.update(dict.fromkeys(difHash, [bc[0]]))
                    for sameHashKey in sameHash: bcHashDict[sameHashKey].extend([bc[0]])
                    bcDict[bc[0]] = [bc]
            else:
                # False
                bcHashDict.update(dict.fromkeys(bcKeyListOne, [bc[0]]))
                bcDict[bc[0]] = [bc]
        bar.finish()
    return bcDict

def SelectionReliableBarcodeMutation(bcMutCount, bcDict):
    mDict = {}
    for x, y in bcMutCount.items():
        if mDict.get(x[0]) == None:
            mDict[x[0]] = [(x[1], y)]
        else:
            mDict[x[0]].append((x[1], y))
    seqDict = {x: {y[0]: mDict[y[0]] for y in bcDict[x]} for x in bcDict.keys()}
    return seqDict

def SelectionReliableBarcodeGenome(bcGenomeCount, bcDict):
    gDict = {}
    for x, y in bcGenomeCount.items():
        if gDict.get(x[0]) == None:
            gDict[x[0]] = [(x[1], y)]
        else:
            gDict[x[0]].append((x[1], y))
    bcDictAlign = {k: [(bc[0], bc[1]) for bc in v if bc[0] in gDict] for k, v in bcDict.items() if k in gDict}
    seqDict = {x: {y[0]: gDict[y[0]] for y in bcDictAlign[x]} for x in bcDictAlign.keys()}
    return bcDictAlign, seqDict

def SelectTheMostProbableMutation(seqDict, mutationProbability=0.99):
    resultDict = {}
    for barcodeID in seqDict:
        numberMut = sum([len(x) for x in seqDict[barcodeID].values()])
        if numberMut > 1:
            # Count number of mutations for unique barcode
            mutLi = {}
            for key in seqDict[barcodeID].values():
                for m in key:
                    if mutLi.get(m[0]) != None:
                        mutLi[m[0]] += m[1]
                    else:
                        mutLi[m[0]] = m[1]
            # Get frequently encountered mutation and enter this in 99% percentile
            mutCoverage = round(float(max(mutLi.values())) / float(sum(mutLi.values())), 2)
            if mutCoverage >= mutationProbability:
                mutZero, value = max(mutLi.iteritems(), key=lambda x:x[1])
                MutatedBarcodes = {key: value for key, value in seqDict[barcodeID].iteritems() if key != barcodeID}
                resultDict[barcodeID] = [mutZero, value, MutatedBarcodes]
        elif numberMut == 1:
            resultDict[barcodeID] = [seqDict[barcodeID].values()[0][0][0], seqDict[barcodeID].values()[0][0][1]]
    return resultDict