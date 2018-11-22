#!/usr/bin/env python
# encoding: utf-8
import os
import string
from itertools import islice
from toolshed import nopen
"""
**** PAIRED-END FUNCTIONS
"""


def read_fastq(fq):
    r"""fastq parser that returns name, seq, and qual."""
    while True:
        values = list(islice(fq, 4))
        if len(values) == 4:
            id1, seq, id2, qual = values
        elif len(values) == 0:
            raise StopIteration
        else:
            raise EOFError("unexpected end of file")
        assert id1.startswith('@'),\
            ">> Fastq out of sync at read:\n%s\n" % id1
        assert id2.startswith('+'),\
            ">> Fastq out of sync at read:\n%s\n" % id1
        assert len(seq) == len(qual),\
            ">> Sequence and Quality are not the same length \
            for read:\n%s\n" % id1
        yield id1[1:-1], seq[:-1], qual[:-1]


def fastqtodict(fastq, separator):
    """returns dict of read name to sequence"""
    fdict = {}
    with nopen(fastq) as fq:
        for name, seq, qual in read_fastq(fq):
            # explicitly state space to facilitate future changes
            fdict[name.split(separator)[0]] = [seq, qual]
    return fdict


def reverseComplement(x):
    complement = string.maketrans('ACGTNSRYMKWHBVD', 'TGCANSRYMKWHBVD')
    x = x.translate(complement)[::-1]
    return x


def FastqJoinPaired(r1, r2, output_dir, gap_size, separator, mode="paired", reverse_complement=False):
    # set dictionary based on mode
    if mode == "R2":
        fastqdict = fastqtodict(r1, separator)
        fastq = r2
    else:
        fastqdict = fastqtodict(r2, separator)
        fastq = r1
    gap_bind, gap_qual = "N" * int(gap_size), "*" * int(gap_size)
    p_out, unq_out = os.path.join(output_dir, "output_paired.fastq"), os.path.join(
        output_dir, "output_unique.fastq")
    with nopen(fastq) as fq:
        with open(p_out, "w") as handle_p:
            with open(unq_out, "w") as handle_unq:
                for name, seq, qual in read_fastq(fq):
                    try:
                        # explicitly state space to facilitate future changes
                        name = name.split(" ")[0]
                        cseq = fastqdict.get(name)[0]
                        cqual = fastqdict.get(name)[1]
                        if reverse_complement:
                            cseq = reverseComplement(cseq)
                        handle_p.write(
                            "@{}\n{}{}{}\n+\n{}{}{}\n".format(name, seq, gap_bind, cseq, qual, gap_qual, cqual))
                    except KeyError:
                        # without pairs
                        if not mode == "paired":
                            handle_unq.write(
                                "@{}\n{}\n+\n{}\n".format(name, seq, qual))
    return p_out
