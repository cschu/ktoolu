import os
import sys
import csv
import argparse

# import ktoolu_io as KTIO
from ktoolu_taxonomy import ktTaxonomyTree as ktree

def filterSequences(db, _in_class_stream, keepTaxIDs, allowUnclassified=False, logfile=None, class_fmt="kraken"):
    if logfile is not None:
        [logfile.write('Filtering sequences...\n'), logfile.flush()]
    assert keepTaxIDs or allowUnclassified
    nseqs = 0
    keepSequences = set()
    # need to keep track of all dropped sequences in case we want unclassified but used the --only-classified-output switch 
    # when running kraken (unclassified sequences are unmarked and need to be distinguished from unwanted sequences.)
    dropSequences = set()
    n_unclassified_keep, n_classified_keep, n_unclassified_discard, n_classified_discard = 0, 0, 0, 0
    
    with _in_class_stream as fi: 
        if class_fmt == "tabular":
            for nseqs, line in enumerate(fi, start=1):
                line = line.strip().split()
                takeClassified = line[1].isdigit() and int(line[1]) in keepTaxIDs 
                takeUnclassified = allowUnclassified and line[1] == "N/A"

                if takeUnclassified or takeClassified:
                    keepSequences.add(line[0].strip())
                elif allowUnclassified:
                    dropSequences.add(line[1].strip())
        elif class_fmt == "kraken":
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
        elif class_fmt == "centrifuge":
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




if __name__ == "__main__":
    from ktoolu_taxonomy import compileValidTaxIDs, ktTaxonomyTree
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("taxid", type=int, default=1)
    ap.add_argument("in_format", choices=("kraken", "centrifuge", "tabular"), default="kraken")
    ap.add_argument("taxdb")
    ap.add_argument("--allow-unclassified", action="store_true")
    args = ap.parse_args()

    keepTaxIDs = compileValidTaxIDs(args.taxdb, wantedTaxIDs=[args.taxid])
    with open(args.input) as _in_class_stream:
        keep, drop = filterSequences(args.taxdb, _in_class_stream, keepTaxIDs, class_fmt=args.in_format, allowUnclassified=args.allow_unclassified)
    for seq in keep:
        print(seq)
    
