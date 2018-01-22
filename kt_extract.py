#!/usr/bin/env python
#
# Copyright (c) 2014-2016 Christian Schudoma, The Sainsbury Laboratory
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import os
import sys
import csv
import argparse

# import ktoolu_io as KTIO
from ktoolu_taxonomy import ktTaxonomyTree as ktree

import ktio.ktio as KTIO
# from ktoolu.ktoolu_taxonomy import ktTaxonomyTree as ktree



def compileValidTaxIDs(db, wantedTaxIDs=[], unwantedTaxIDs=[], vipTaxIDs=[], logfile=None):
    keepTaxIDs = set()

    if logfile is not None:
        [logfile.write('Reading taxonomy...\n'), logfile.flush()]
    taxTree = ktree(db)
    if logfile is not None:
        [logfile.write('Traversing requested taxonomy branch(es)...\n'), logfile.flush()]

    # print(wantedTaxIDs)
    for taxID in wantedTaxIDs:
        keepTaxIDs.update(taxTree.getDescendents(abs(taxID)))
    # print(unwantedTaxIDs)
    for taxID in unwantedTaxIDs:
        keepTaxIDs.difference_update(taxTree.getDescendents(abs(taxID)))

    return keepTaxIDs.union(set(vipTaxIDs))

def filterSequences(db, f_inputClassification, keepTaxIDs, allowUnclassified=False, logfile=None, isKraken=True):
    if logfile is not None:
        [logfile.write('Filtering sequences...\n'), logfile.flush()]
    assert keepTaxIDs or allowUnclassified
    nseqs = 0
    keepSequences = set()
    # need to keep track of all dropped sequences in case we want unclassified but used the --only-classified-output switch 
    # when running kraken (unclassified sequences are unmarked and need to be distinguished from unwanted sequences.)
    dropSequences = set()
    n_unclassified_keep, n_classified_keep, n_unclassified_discard, n_classified_discard = 0, 0, 0, 0
    
    with KTIO.openFile(f_inputClassification) as fi:
        if isKraken:
            for line in fi:
                nseqs += 1
                line = line.strip().split()
                takeClassified = line[0] == 'C' and int(line[2]) in keepTaxIDs
                takeUnclassified = allowUnclassified and line[0] == 'U'
                # print(line, takeClassified, takeUnclassified)
    
                if takeUnclassified or takeClassified:
                    keepSequences.add(line[1].strip())
                elif allowUnclassified:
                    # we want unclassified, but current line was classified and rejected
                    # myself 2.something years later: there is probably a reason for doing it this way...
                    dropSequences.add(line[1].strip())
        else:
            next(fi)
            for nseqs, row in enumerate(csv.reader(fi, delimiter='\t'), start=1):
                keep = False
                isUnclassified = row[1] == 'unclassified'
                if isUnclassified:
                    if allowUnclassified:
                        keep = True
                        n_unclassified_keep += 1
                    else:
                        n_unclassified_discard += 1
                else:
                    if int(row[2]) in keepTaxIDs:
                        keep = True
                        n_classified_keep += 1
                    else:
                        n_classified_discard += 1
                if keep:
                    keepSequences.add(row[0])
                elif allowUnclassified:
                    dropSequences.add(row[0])

                """
                takeClassified = row[1] != 'unclassified' and int(row[2]) in keepTaxIDs
                takeUnclassified = allowUnclassified and row[1] == 'unclassified'

                if takeUnclassified or takeClassified:
                    keepSequences.add(row[0])
                elif allowUnclassified:
                    dropSequences.add(row[0])
                """

        if logfile is not None:
            logfile.write('Keeping %i of %i sequences (%.1f).\n' % (len(keepSequences), nseqs, float(len(keepSequences))/nseqs))
            logfile.write('Dropping %i of %i sequences (%.1f).\n' % (len(dropSequences), nseqs, float(len(dropSequences))/nseqs))
   
            logfile.write('C:keep={} ({}), C:discard={} ({}), U:keep={} ({}), U:discard={} ({})\n'.format(n_classified_keep, n_classified_keep/nseqs, n_classified_discard, n_classified_discard/nseqs, n_unclassified_keep, n_unclassified_keep/nseqs, n_unclassified_discard, n_unclassified_discard/nseqs))
            logfile.flush()

    return keepSequences, dropSequences


def extractSequences_obsolete(keepSequences, fileInfo):
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
    parser.add_argument('--centrifuge-mode', action='store_true')
    args = parser.parse_args()

    assert args.db and os.path.exists(args.db)
    assert args.kraken_results and os.path.exists(args.kraken_results)
    input_exists = args.inR1 and os.path.exists(args.inR1)
    fformat_matches = KTIO.verifyFileFormat(args.inR1, args.input_format)
    assert input_exists and fformat_matches and args.outR1
    input_exists = args.inR2 and os.path.exists(args.inR2)
    fformat_matches = not input_exists or KTIO.verifyFileFormat(args.inR2, args.input_format) 
    assert (not input_exists) or (input_exists and fformat_matches and args.outR2)

    def xor(a,b):
        return (a and not b) or (not a and b)
    assert xor(args.gz_output, args.bz2_output) or not(args.gz_output or args.bz2_output)

    try:
        # Let's see if we have one or more roots specified to extract taxonomic subtrees.
        wantedTaxIDs = list(map(int, args.keep_taxids.replace(' ', '').split(',')))
    except:
        # If not,
        if args.include_unclassified:
            # and include-unclassified is True, then we assume we only want unclassified sequences.
            wantedTaxIDs = []
        else:
            # Otherwise just take the whole tree.
            wantedTaxIDs = [1]
    try:
        unwantedTaxIDs = list(map(int, args.drop_taxids.replace(' ', '').split(',')))
    except:
        unwantedTaxIDs = []
    try:
        vipTaxIDs = list(map(int, args.vip_taxids.replace(' ', '').split(',')))
    except:
        vipTaxIDs = []

    logfile = sys.stdout
    keepTaxIDs = compileValidTaxIDs(args.db, wantedTaxIDs=wantedTaxIDs, unwantedTaxIDs=unwantedTaxIDs, vipTaxIDs=vipTaxIDs, logfile=logfile)
    keepSequences, dropSequences = filterSequences(args.db, args.kraken_results, keepTaxIDs, allowUnclassified=args.include_unclassified, isKraken=(not args.centrifuge_mode), logfile=logfile)
    KTIO.extractSequences(keepSequences, args, rejected=dropSequences)
    print("DONE", file=logfile)
    pass

if __name__ == '__main__': main(sys.argv[1:])


__author__ = "Christian Schudoma"
__copyright__ = "Copyright 2014-2016, Christian Schudoma, The Sainsbury Laboratory"
__credits__ = ["Pirasteh Pahlavan", "Agathe Jouet", "Yogesh Gupta"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Christian Schudoma"
__email__ = "cschu1981@gmail.com"
__status__ = "Development"
