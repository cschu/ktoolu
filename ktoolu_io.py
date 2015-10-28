#!/usr/bin/env python
import sys
import unittest

def isPreCassava18(string):
    return string.endswith('/1') or string.endswith('/2')
def getFastqIdentifier(string):
    return string[1:-2] if isPreCassava18(string) else string.split()[0][1:]
def getFastaIdentifier(string):
    string = string.split()[0][1:]
    return string[:-2] if isPreCassava18(string) else string    
def verifyFileFormat(fn, fileFormat):
    firstChar = open(fn).read(1)
    verifiedFastq = firstChar == '@' and fileFormat == 'fq'
    verifiedFasta = firstChar == '>' and fileFormat == 'fa'
    return verifiedFastq or verifiedFasta


def readFasta(fn):
    """
    Provides a generator to the contents of a Fasta file.
    The function can handle single seq/multi seq, formatted/unformatted Fasta.

    Input: a filename
    Output: a generator object providing identifier and sequence information
    for each record in the Fastq file
    """
    with open(fn) as fi:
        seq, identifier = [], ''
        for line in fi:
            line = line.strip()
            if not line:
                # ignore empty lines
                continue
            if line.startswith('>'):
                # new record starts with >identifier
                if seq:
                    # if there is data in seq, return id, seq and empty seq
                    yield (identifier, ''.join(seq))
                    seq = []
                identifier = line
            else:
                # current line contains sequence information
                seq.append(line)
        if seq:
            # return last sequence record
            yield (identifier, ''.join(seq))

def readFastq(fn):
    """
    Provides a generator to the contents of a Fastq file.
    The function can handle single seq/multi seq, formatted/unformatted Fastq.

    Input: a filename
    Output: a generator object providing identifier, sequence, and quality information
    for each record in the Fastq file
    """
    with open(fn) as fi:
        block, identifier = [], ''
        for line in fi:
            line = line.strip()
            if not line:
                # ignore empty lines
                continue
            if line.startswith('@'):
                # new record starts with @identifier
                identifier = line
            elif line.startswith('+'):
                # if seq/qual separator occurs,
                # read equal number of lines of quality information
                seq = ''.join(block)
                qual = ''.join([fi.next().strip() for row in block])
                yield (identifier, seq, qual)
                block, head = [], ''
            else:
                # current line contains sequence information
                block.append(line)





class ktFastxIterator(object):
    def __init__(self, **kwargs):
        fmt = kwargs.get('fmt', 'fq')
        assert fmt in ('fq', 'fa')
        if fmt == 'fa':
            self.getID = getFastaIdentifier
            self.getSeqs = readFasta
        else:
            self.getID = getFastqIdentifier
            self.getSeqs = readFastq

        assert 'in1' in kwargs and 'out1' in kwargs
        self.in1, self.out1 = self.getSeqs(kwargs['in1']), open(kwargs['out1'], 'wb')
        self.in2, self.out2 = kwargs.get('in2', None), kwargs.get('out2', None)
        if self.in2 is not None and self.out2 is not None:
            self.in2, self.out2 = self.getSeqs(self.in2), open(self.out2, 'wb')
        pass
    def next(self):
        assert self.in1
        y1, y2 = self.in1.next(), (self.in2.next() if self.in2 is not None else None)
        print 'y1', y1, type(y1)
        if y1:
            yield y1, y2
        #while 1:
        #    try:
        #         yield self.in1.next(), (self.in2.next() if self.in2 is not None else None)
        #    except:
        #        break

    pass


"""

    fwdOut = open(outputR1, 'wb')
    fwdGen = getSeqs(inputR1) # anabl_getSeqsFromFastX(inputR1, X=fastx)

    revOut, revGen = None, None
    revSid, revSeq = None, None
    if outputR2 is not None and inputR2 is not None:
        revOut = open(outputR2, 'wb')
        revGen = getSeqs(inputR2) # anabl_getSeqsFromFastX(inputR2, X=fastx)


    fxid1, fxid2 = None, None
    while True:
        try:
            fwdSid, fwdSeq = fwdGen.next()
        except:
            break
        logfile.write('*%s*\n' % fwdSid)
        fxid1 = getID(fwdSid)
        if revGen is not None:
            try:
                revSid, revSeq = revGen.next()
            except:
                break
            fxid2 = getID(revSid)

        if fxid1 != fxid2 and fxid2 is not None:
            sys.stderr.write('Error: fxid-mismatch %s %s.\n' % (fxid1, fxid2))
            sys.exit(1)
        if fxid1 in keepSequences:
            fwdOut.write(('%s\n' * fastx) % ((fwdSid,) + fwdSeq))
            if revOut is not None:
                revOut.write(('%s\n' * fastx) % ((revSid,) + revSeq))
        else:
            # sys.stdout.write('%s is not in keepSequences\n' % fqid1)
            pass
    fwdOut.close()
    if revOut is not None:
        revOut.close()
    pass













"""





"""
Unit tests
"""

FASTQ_TESTFILE_R1 = 'testdata/ktoolu_test.R1.fq'
FASTQ_TESTRECORDS = [22, 27, 32]
FASTA_TESTFILE = 'testdata/ktoolu_test.fa'
FASTA_TESTRECORDS = [22, 27, 32]

class readFastqTest(unittest.TestCase):
    def test(self):
        data = list(readFastq(FASTQ_TESTFILE_R1))
        self.assertEqual(len(data), len(FASTQ_TESTRECORDS))
        self.assertEqual(len(data[0][1]), FASTQ_TESTRECORDS[0])
        self.assertEqual(len(data[1][1]), FASTQ_TESTRECORDS[1])
        self.assertEqual(len(data[2][1]), FASTQ_TESTRECORDS[2])
class readFastaTest(unittest.TestCase):
    def test(self):
        data = list(readFasta(FASTA_TESTFILE))
        self.assertEqual(len(data), len(FASTA_TESTRECORDS))
        self.assertEqual(len(data[0][1]), FASTQ_TESTRECORDS[0])
        self.assertEqual(len(data[1][1]), FASTQ_TESTRECORDS[1])
        self.assertEqual(len(data[2][1]), FASTQ_TESTRECORDS[2])
class getFastqIdentifierTest(unittest.TestCase):
    def test(self):
        self.assertEqual(getFastqIdentifier('@HWUSI-EAS100R:6:73:941:1973#0/1'), 'HWUSI-EAS100R:6:73:941:1973#0')
        self.assertEqual(getFastqIdentifier('@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG'), 'EAS139:136:FC706VJ:2:2104:15343:197393')

if __name__ == '__main__':
    unittest.main()
    pass
