#!/usr/bin/env python
import sys
import argparse

import ktoolu_io as KTIO

def readClassification(fn):
    with open(fn) as fi:
        return set(line.strip().split()[:2] for line in fi if line.strip().startswith('C'))

def computeSets(resultsA, resultsB):
    setA = readClassification(resultsA)
    setB = readClassification(resultsB)
    setAB = setA.intersection(setB)
    setA.difference_update(setB)
    setB.difference_update(setA)
    return setA, setB, setAB

def assignSequences(sets, fileInfo):
    assert fileInfo.input_format in ('fq', 'fa')
    if fileInfo.input_format == 'fq':
        getID, getSeqs, nlines = KTIO.getFastqIdentifier, KTIO.readFastq, 4
    else:
        getID, getSeqs, nlines = KTIO.getFastaIdentifier, KTIO.readFasta, 2

    R1gen, R2gen = getSeqs(fileInfo.inR1), None
    R1A, R1B, R1AB, R1U = map(lambda x:open(x, 'wb'), [fileInfo.outAR1, fileInfo.outBR1, fileInfo.outABR1, fileInfo.outUR1])
    R2A, R2B, R2AB, R2U = None, None, None, None

    if inR2 is not None:
        R2gen = getSeqs(fileInfo.inR2)
        R2A, R2B, R2AB, R2U = map(lambda x:open(x, 'wb'), [fileInfo.outAR2, fileInfo.outBR2, fileInfo.outABR2, fileInfo.outUR2])

    fxid1, fxid2 = None, None
    while 1:
        try:
            R1rec = R1gen.next()
        except:
            break
        fxid1, fxid2 = getID(R1rec[0]), None

        if R2gen is not None:
            try:
                R2rec = R2gen.next()
            except:
                break
            fxid2 = getID(R2rec[0])

        assert fxid1 == fxid2 or fxid2 is None

        # set order is A, B, AB
        if fxid1 in sets[0]:
            dest = R1A, R2A
        elif fxid1 in sets[1]:
            dest = R1B, R2B
        elif fxid1 in sets[2]:
            dest = R1AB, R2AB
        else:
            dest = R1U, R2U

        dest[0].write(('%s\n' * nlines) % R1rec
        if dest[1] is not None:
            dest[1].write(('%s\n' * nlines) % R2rec

    map(lambda x:x.close(), [R1A, R1B, R1AB, R1U])
    if R2A is not None:
        map(lambda x:x.close(), [R2A, R2B, R2AB, R2U)

    pass


def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input-format', help='Input sequences stored in Fasta (fa) or Fastq (fq) file(s).', default='fq')
    parser.add_argument('--inR1', help='The R1-file (single-end reads or forward paired-end reads).')
    parser.add_argument('--inR2', help='The R2-file (reverse paired-end reads)')
    parser.add_argument('--kraken-resultsA', type=str, help='A file containing the results of classification A for the input sequences.')
    parser.add_argument('--kraken-resultsB', type=str, help='A file containing the results of classification B for the input sequences.')
    parser.add_argument('--outAR1', type=str, help='')
    parser.add_argument('--outAR2', type=str, help='')
    parser.add_argument('--outBR1', type=str, help='')
    parser.add_argument('--outBR2', type=str, help='')
    parser.add_argument('--outABR1', type=str, help='')
    parser.add_argument('--outABR2', type=str, help='')
    parser.add_argument('--outUnclassifiedR1', type=str, help='')
    parser.add_argument('--outUnclassifiedR2', type=str, help='')
    args = parser.parse_args()

    assert 'kraken_resultsA' in args and os.path.exists(args.kraken_resultsA)
    assert 'kraken_resultsB' in args and os.path.exists(args.kraken_resultsB)
    assert args.inR1 and os.path.exists(args.inR1) and verifyFileFormat(args.inR1, args.input_format)
    assert (not args.inR2) or (os.path.exists(args.inR2) and verifyFileFormat(args.inR2, args.input_format))
    assert args.outAR1 is not None and args.outBR1 is not None and args.outABR1 is not None and args.outUR1 is not None
    assert (not args.inR2) or (args.outAR2 is not None and args.outBR2 is not None and args.outABR2 is not None and args.outUR2 is not None)

    sets = computeSets(args.kraken_resultsA, args.kraken_resultsB)
    assignSequences(sets, args)








    pass

if __name__ == '__main__': main()
