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
    file_exists = os.path.exists(fn)
    assert mode in ('a', 'ab', 'r', 'rb', 'w', 'wb')
    assert file_exists or mode in ('w', 'wb', 'a', 'ab')
    assert fmt in (None, 'gz', 'bz2')
    if fmt == 'gz' or (fmt is None and (file_exists and isGZ(fn))):
        return gzip.open(fn, mode)
    elif fmt == 'bz2' or (fmt is None and (file_exists and isBZ2(fn))):
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


def readFasta(_in):
    if type(_in) is str:
        with openFile(_in) as fi:
            return processFasta(fi)
    else:
        return processFasta(_in)

def processFasta(fi):
    """
    Provides a generator to the contents of a Fasta file.
    The function can handle single seq/multi seq, formatted/unformatted Fasta.

    Input: a file/stream handle
    Output: a generator object providing identifier and sequence information
    for each record in the Fastq file
    """
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

def readFastq(_in):
    if type(_in) is str:
        with openFile(_in) as fi:
            return processFastq(fi)
    else:
        return processFastq(_in)

def processFastq(fi):
    """
    Provides a generator to the contents of a Fastq file.
    The function can handle single seq/multi seq, formatted/unformatted Fastq.

    Input: a file/stream handle
    Output: a generator object providing identifier, sequence, and quality information
    for each record in the Fastq file
    """
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


def extractSequences(keepSequences, fileInfo, rejected=set()):
    assert fileInfo.input_format in ('fq', 'fa')
    if fileInfo.input_format == 'fq':
        getID, getSeqs, nlines, outfmt = getFastqIdentifier, readFastq, 4, '%s\n%s\n+\n%s\n'
    else:
        getID, getSeqs, nlines, outfmt = getFastaIdentifier, readFasta, 2, '%s\n%s\n'

    if 'gz_output' in fileInfo and fileInfo.gz_output:
        ffmt = 'gz'
    elif 'bz2_output' in fileInfo and fileInfo.bz2_output:
        ffmt = 'bz2'
    else:
        ffmt = None
    
    fwdOut, fwdGen = openFile(fileInfo.outR1, mode='wb', fmt=ffmt), getSeqs(fileInfo.inR1)
    revOut, revGen = None, None

    if fileInfo.outR2 and fileInfo.inR2:
        revOut, revGen = openFile(fileInfo.outR2, mode='wb', fmt=ffmt), getSeqs(fileInfo.inR2)

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

        if fxid1 in keepSequences or (rejected and fxid1 not in rejected):  # the rejected-set allows extracting unclassified sequences if they are not in the kraken-output
            fwdOut.write(outfmt % fwdRecord)
            if revOut is not None:
                revOut.write(outfmt % revRecord)
        else:
            pass
    fwdOut.close()
    if revOut is not None:
        revOut.close()
    pass



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


__author__ = "Christian Schudoma"
__copyright__ = "Copyright 2014-2016, Christian Schudoma, The Sainsbury Laboratory"
__credits__ = ["Pirasteh Pahlavan", "Agathe Jouet", "Yogesh Gupta"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Christian Schudoma"
__email__ = "cschu1981@gmail.com"
__status__ = "Development"
