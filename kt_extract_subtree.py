#/usr/bin/env python
import sys

import ktoolu_io as KTIO
from ktoolu_taxonomy import ktTaxonomyTree as ktree

def filterSequences(db, taxID, f_inputClassification, logfile=None):
    keepSequences = set()
    if taxID != 0:
        if logfile is not None:
            [logfile.write('Reading taxonomy...\n'), logfile.flush()]
        taxTree = ktree(db)
        if logfile is not None:
            [logfile.write('Traversing requested taxonomy branch...\n'), logfile.flush()]
        subtreeTaxIDs = taxTree.getDescendents(abs(taxID))
        if logfile is not None:
            [logfile.write('Family tree of %s has %i descendents.\n' % (str(abs(taxID)), len(subtreeTaxIDs))), logfile.flush()]
    else:
        subtreeTaxIDs = set()

    if logfile is not None:
        [logfile.write('Filtering sequences...\n'), logfile.flush()]

    nseqs = 0
    with open(f_inputClassification) as fi:
        for line in fi:
            nseqs += 1
            line = line.strip().split()

            if taxID > 0:
                # extract reads from branch
                # don't allow unclassified reads
                takeUnclassified = False
                # only allow reads that have been assigned a taxonomy id belonging to the branch
                takeClassified = line[0] == 'C' and int(line[2]) in subtreeTaxIDs
            elif taxID < 0:
                # ignore reads from branch
                # allow unclassified reads
                takeUnclassified = line[0] == 'U'
                # allow only reads that have been assigned a taxonomy id outside of the branch
                takeClassified = line[0] == 'C' and int(line[2]) not in subtreeTaxIDs
            else:
                # extract unclassified
                # allow unclassified reads
                takeUnclassified = taxID == 0 and line[0] == 'U'
                # don't allow classified reads
                takeClassified = False

            if takeUnclassified or takeClassified:
                keepSequences.add(line[1].strip())

        if logfile is not None:
            [logfile.write('Keeping %i of %i sequences (%.1f).\n' % (len(keepSequences), nseqs, float(len(keepSequences))/nseqs)), logfile.flush()]

    return keepSequences


def extractSequences(keepSequences, inputR1, outputR1, inputR2=None, outputR2=None, fmt='fq'):
    assert fmt in ('fq', 'fa')
    if fmt == 'fq':
        getID, getSeqs, nlines = KTIO.getFastqIdentifier, KTIO.readFastq, 4
    else:
        getID, getSeqs, nlines = KTIO.getFastaIdentifier, KTIO.readFasta, 2

    fwdOut, fwdGen = open(outputR1, 'wb'), getSeqs(inputR1)
    revOut, revGen = None, None
    revSid, revSeq = None, None

    if outputR2 is not None and inputR2 is not None:
        revOut, revGen = open(outputR2, 'wb'), getSeqs(inputR2)

    fxid1, fxid2 = None, None
    while True:
        try:
            rwdRecord = fwdGen.next()
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
            fwdOut.write(('%s\n' * nlines) % fwdRecord
            if revOut is not None:
                revOut.write(('%s\n' * nlines) % revRecord
        else:
            pass
    fwdOut.close()
    if revOut is not None:
        revOut.close()

    pass
