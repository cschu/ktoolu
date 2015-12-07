#!/usr/bin/env python
import sys
import argparse

import ktoolu_io as KTIO

def readClassification(fn):
    with KTIO.openFile(fn) as fi:
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
        getID, getSeqs, nlines, outfmt = KTIO.getFastqIdentifier, KTIO.readFastq, 4, '%s\n%s\n+\n%s\n'
    else:
        getID, getSeqs, nlines, outfmt = KTIO.getFastaIdentifier, KTIO.readFasta, 2, '%s\n%s\n'

    if fileInfo.gz_output:
        ffmt = 'gz'
    elif fileInfo.bz2_output:
        ffmt = 'bz2'
    else:
        ffmt = None

    R1gen, R2gen = getSeqs(fileInfo.inR1), None
    R1A, R1B, R1AB, R1U = map(lambda x:KTIO.openFile(x, fmt=ffmt, mode='wb'), [fileInfo.outAR1, fileInfo.outBR1, fileInfo.outABR1, fileInfo.outUR1])
    R2A, R2B, R2AB, R2U = None, None, None, None

    if fileInfo.inR2 is not None:
        R2gen = getSeqs(fileInfo.inR2)
        R2A, R2B, R2AB, R2U = map(lambda x:KTIO.openFile(x, fmt=ffmt, mode='wb'), [fileInfo.outAR2, fileInfo.outBR2, fileInfo.outABR2, fileInfo.outUR2])

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
            dest, destid = (R1A, R2A), 'A'
        elif fxid1 in sets[1]:
            dest, destid = (R1B, R2B), 'B'
        elif fxid1 in sets[2]:
            dest, destid = (R1AB, R2AB), '+'
        else:
            dest, destid = (R1U, R2U), 'U'

        sys.stdout.write('\t'.join([destid, fxid1]) + '\n')

        dest[0].write(outfmt % R1rec)
        if dest[1] is not None:
            dest[1].write(outfmt % R2rec)

    map(lambda x:x.close(), [R1A, R1B, R1AB, R1U])
    if R2A is not None:
        map(lambda x:x.close(), [R2A, R2B, R2AB, R2U])

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
    parser.add_argument('--gz-output', action='store_true')
    parser.add_argument('--bz2-output', action='store_true')
    args = parser.parse_args()

    assert args.kraken_resultsA is not None and os.path.exists(args.kraken_resultsA)
    assert args.kraken_resultsB is not None and os.path.exists(args.kraken_resultsB)
    assert args.inR1 is not None and os.path.exists(args.inR1) and verifyFileFormat(args.inR1, args.input_format)
    assert (args.inR2 is None) or (os.path.exists(args.inR2) and verifyFileFormat(args.inR2, args.input_format))
    assert args.outAR1 is not None and args.outBR1 is not None and args.outABR1 is not None and args.outUR1 is not None
    assert (args.inR2 is None) or (args.outAR2 is not None and args.outBR2 is not None and args.outABR2 is not None and args.outUR2 is not None)
    def xor(a,b):
        return (a and not b) or (not a and b)
    assert xor(args.gz_output, args.bz2_output) or not(args.gz_output or args.bz2_output)

    sets = computeSets(args.kraken_resultsA, args.kraken_resultsB)
    assignSequences(sets, args)
    pass

if __name__ == '__main__': main()
