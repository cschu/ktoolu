#/usr/bin/env python
import os
import sys
import argparse

import ktoolu_io as KTIO
from ktoolu_taxonomy import ktTaxonomyTree as ktree

def compileValidTaxIDs(db, wantedTaxIDs=[], unwantedTaxIDs=[], vipTaxIDs=[], logfile=None):
    keepTaxIDs = set()

    if logfile is not None:
        [logfile.write('Reading taxonomy...\n'), logfile.flush()]
    taxTree = ktree(db)
    if logfile is not None:
        [logfile.write('Traversing requested taxonomy branch(es)...\n'), logfile.flush()]

    for taxID in wantedTaxIDs:
        keepTaxIDs.update(taxTree.getDescendents(abs(taxID)))
    for taxID in unwantedTaxIDs:
        keepTaxIDs.difference_update(taxTree.getDescendents(abs(taxID)))

    return keepTaxIDs.union(set(vipTaxIDs))

def filterSequences(db, f_inputClassification, keepTaxIDs, allowUnclassified=False, logfile=None):
    if logfile is not None:
        [logfile.write('Filtering sequences...\n'), logfile.flush()]
    assert keepTaxIDs or allowUnclassified
    nseqs = 0
    keepSequences = set()

    with KTIO.openFile(f_inputClassification) as fi:
        for line in fi:
            nseqs += 1
            line = line.strip().split()

            """
            if taxID > 0:
                # extract reads from branch
                # don't allow unclassified reads
                takeUnclassified = allowUnclassified and line[0] == 'U'
                # only allow reads that have been assigned a taxonomy id belonging to the branch
                takeClassified = line[0] == 'C' and int(line[2]) in keepTaxIDs
            elif taxID < 0:
                # ignore reads from branch
                # do not allow unclassified reads
                takeUnclassified = allowUnclassified and line[0] == 'U'
                # allow only reads that have been assigned a taxonomy id outside of the branch
                takeClassified = line[0] == 'C' and int(line[2]) not in keepTaxIDs
            else:
                # extract unclassified
                # allow unclassified reads
                takeUnclassified = allowUnclassified and line[0] == 'U'
                # don't allow classified reads
                takeClassified = False
            """
            takeClassified = line[0] == 'C' and int(line[2]) in keepTaxIDs
            takeUnclassified = allowUnclassified and line[0] == 'U'

            if takeUnclassified or takeClassified:
                keepSequences.add(line[1].strip())

        if logfile is not None:
            [logfile.write('Keeping %i of %i sequences (%.1f).\n' % (len(keepSequences), nseqs, float(len(keepSequences))/nseqs)), logfile.flush()]

    return keepSequences


def extractSequences(keepSequences, fileInfo):
    assert fileInfo.input_format in ('fq', 'fa')
    if fileInfo.input_format == 'fq':
        getID, getSeqs, nlines = KTIO.getFastqIdentifier, KTIO.readFastq, 4
    else:
        getID, getSeqs, nlines = KTIO.getFastaIdentifier, KTIO.readFasta, 2

    if fileInfo.gz_output:
        ffmt = 'gz'
    elif fileInfo.bz2_output:
        ffmt = 'bz2'
    else:
        ffmt = None

    fwdOut, fwdGen = KTIO.openFile(fileInfo.outR1, mode='wb', fmt=ffmt), getSeqs(fileInfo.inR1)
    revOut, revGen = None, None

    if args.outR2 is not None and args.inR2 is not None:
        revOut, revGen = KTIO.openFile(args.outR2, mode='wb', fmt=ffmt), getSeqs(args.inR2)

    fxid1, fxid2 = None, None
    while 1:
        try:
            fwdRecord = fwdGen.next()
        except:
            break
        fxid1 = getID(fwdRecord[0])
        if revGen is not None:
            try:
                revRecord = revGen.next()
            except:
                break
            fxid2 = getID(revRecord[0])

        assert fxid1 == fxid2 or fxid2 is None

        if fxid1 in keepSequences:
            fwdOut.write(('%s\n' * nlines) % fwdRecord)
            if revOut is not None:
                revOut.write(('%s\n' * nlines) % revRecord)
        else:
            pass
    fwdOut.close()
    if revOut is not None:
        revOut.close()

    pass



def main(argv):

    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--db', help='Path to the root directory of a kraken database')
    parser.add_argument('--inR1', help='The r1-file (single-end reads or forward paired-end reads).')
    parser.add_argument('--inR2', help='The r2-file (reverse paired-end reads)')
    parser.add_argument('--keep-taxids', type=str, default='', help='A comma-separated list of taxonomy ids. These ids and all descending ids will be kept unless they are specified in --drop-taxids.')
    parser.add_argument('--drop-taxids', type=str, default='', help='A comma-separated list of taxonomy ids. These ids and all descending ids will be dropped.')
    parser.add_argument('--vip-taxids', type=str, default='', help='A comma-separated list of taxonomy ids. These ids and all descending ids will always be kept.')

    parser.add_argument('--kraken-results', type=str, help='A file containing kraken classification results for the input sequences.')
    parser.add_argument('--input-format', help='Input sequences stored in Fasta (fa) or Fastq (fq) file(s).', default='fq')
    parser.add_argument('--include-unclassified', action='store_true', help='Extract unclassified sequences.')

    parser.add_argument('--outR1', type=str, help='The r1-output file.')
    parser.add_argument('--outR2', type=str, help='The r2-output file.')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--logfile', type=str, help='A logfile.', default='kt_extract.log')

    parser.add_argument('--gz-output', action='store_true')
    parser.add_argument('--bz2-output', action='store_true')
    args = parser.parse_args()

    # for k in sorted(args):
    #    print '%s\t%s' % (k, args[k])
    # return None
    # print args
    # return None

    assert 'db' in args and os.path.exists(args.db)
    assert 'kraken_results' in args and os.path.exists(args.kraken_results)
    assert args.inR1 and os.path.exists(args.inR1) and KTIO.verifyFileFormat(args.inR1, args.input_format) and args.outR1
    assert (not args.inR2) or (os.path.exists(args.inR2) and KTIO.verifyFileFormat(args.inR2, args.input_format) and args.outR2)

    def xor(a,b):
        return (a and not b) or (not a and b)
    assert xor(args.gz_output, args.bz2_output) or not(args.gz_output or args.bz2_output)

    try:
        wantedTaxIDs = map(int, args.keep_taxids.replace(' ', '').split(','))
    except:
        # by default just take the whole tree
        if args.includeUnclassified:
            wantedTaxIDs = []
        else:
            wantedTaxIDs = [1]
    try:
        unwantedTaxIDs = map(int, args.drop_taxids.replace(' ', '').split(','))
    except:
        unwantedTaxIDs = []
    try:
        vipTaxIDs = map(int, args.vip_taxids.replace(' ', '').split(','))
    except:
        vipTaxIDs = []

    keepTaxIDs = compileValidTaxIDs(args.db, wantedTaxIDs=wantedTaxIDs, unwantedTaxIDs=unwantedTaxIDs, vipTaxIDs=vipTaxIDs)
    keepSequences = filterSequences(args.db, args.kraken_results, keepTaxIDs, allowUnclassified=args.include_unclassified)
    KTIO.extractSequences(keepSequences, args)
    pass

if __name__ == '__main__': main(sys.argv[1:])
