#!/usr/bin/env python
import os
import sys
import bz2
import gzip
import unittest


def isGZ(fn):
    """
    Tests whether a file is gz-compressed.

    :param fn: a filename
    :type fn: str

    :returns: True if fn is gz-compressed otherwise False
    """
    assert os.path.exists(fn)
    with open(fn, 'rb') as fi:
        b1, b2 = fi.read(1), fi.read(1)
        return b1 == '\x1f' and b2 == '\x8b'

def isBZ2(fn):
    """
    Tests whether a file is bz2-compressed.

    :param fn: a filename
    :type fn: str

    :returns: True if fn is bz2-compressed otherwise False
    """
    assert os.path.exists(fn)
    with open(fn, 'rb') as fi:
        return fi.read(10) == 'BZh91AY&SY'

def openFile(fn, fmt=None, mode='rb'):
    """
    Opens an uncompressed, gz-compressed or bz2-compressed file.

    :param fn: a filename
    :type fn: str
    :param fmt: a file compression format
    :type fmt: str {None, 'gz', or 'bz2'}
    :param mode: a file mode (append, read, write)
    :type mode: str {'a', 'ab', 'r', 'rb', 'w', 'wb'}

    :returns: a handle to fn
    """
    assert mode in ('a', 'ab', 'r', 'rb', 'w', 'wb')
    assert os.path.exists(fn) or mode in ('w', 'wb', 'a', 'ab')
    assert fmt in (None, 'gz', 'bz2')
    if fmt == 'gz' or (fmt is None and isGZ(fn)):
        return gzip.open(fn, mode)
    elif fmt == 'bz2' or (fmt is None and isBZ2(fn)):
        return bz2.BZ2File(fn, mode)
    else:
        return open(fn, mode)

def nRecords(fn, fmt='fq'):
    assert fmt in ('fa', 'fq')
    assert os.path.exists(fn)
    nlines = sum(1 for line in openFile)
    return nlines / (4 if fmt == 'fq' else 2)


def isPreCassava18(string):
    return string.endswith('/1') or string.endswith('/2')
def getFastqIdentifier(string):
    return string[1:-2] if isPreCassava18(string) else string.split()[0][1:]
def getFastaIdentifier(string):
    string = string.split()[0][1:]
    return string[:-2] if isPreCassava18(string) else string
def verifyFileFormat(fn, fileFormat):
    firstChar = openFile(fn).read(1)
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
    with openFile(fn) as fi:
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
    with openFile(fn) as fi:
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




"""
# Experimental code
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






"""
Unit tests
"""

GZ_TESTFILE = 'testdata/ktoolu_test.R1.fq.gz'
BZ2_TESTFILE = 'testdata/ktoolu_test.R1.fq.bz2'

FASTQ_TESTFILE_R1 = 'testdata/ktoolu_test.R1.fq'
FASTQ_TESTRECORDS = [22, 27, 32]
FASTA_TESTFILE = 'testdata/ktoolu_test.fa'
FASTA_TESTRECORDS = [22, 27, 32]
class compressedFileTest(unittest.TestCase):
    def test(self):
        self.assertTrue(isGZ(GZ_TESTFILE))
        self.assertFalse(isBZ2(GZ_TESTFILE))
        self.assertTrue(isBZ2(BZ2_TESTFILE))
        self.assertFalse(isGZ(BZ2_TESTFILE))
        self.assertFalse(isGZ(FASTA_TESTFILE))
        self.assertFalse(isGZ(FASTQ_TESTFILE_R1))
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
