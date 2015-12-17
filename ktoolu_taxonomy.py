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


import sys
import os

class ktTaxonomyTree(object):
    def __init__(self, db):
        self._readNames(os.path.join(db, 'taxonomy', 'names.dmp'))
        self._readNodes(os.path.join(db, 'taxonomy', 'nodes.dmp'))
        pass
    def _readNames(self, fn):
        assert os.path.exists(fn)
        self.nodeNames = dict()
        with open(fn) as fi:
            for line in fi:
                line = line.strip().strip('|').strip()
                if not line:
                    break
                line = line.split('\t|\t')
                if line[3].strip() == 'scientific name':
                    self.nodeNames[int(line[0].strip())] = line[1].strip()
        pass
    def _readNodes(self, fn):
        assert os.path.exists(fn)
        self.nodeRanks, self.nodeChildren = dict(), dict()
        with open(fn) as fi:
            for line in fi:
                line = line.strip().strip('|').strip()
                if not line:
                    break
                line = map(lambda x:x.strip(), line.split('\t|\t'))
                line[:2] = map(int, line[:2])
                if line[0] == 1:
                    line[1] = 0
                try:
                    self.nodeChildren[line[1]].add(line[0])
                except:
                    self.nodeChildren[line[1]] = set([line[0]])
                self.nodeRanks[line[0]] = line[2]
        pass
    def getDescendents(self, taxID):
        assert taxID >= 0
        assert self.nodeChildren
        descendents, queue = set([taxID]), [taxID]
        while queue:
            children = self.nodeChildren.get(queue.pop(), set())
            if children:
                descendents.update(children)
                queue.extend(children)
            pass
        return descendents
    pass


__author__ = "Christian Schudoma"
__copyright__ = "Copyright 2014-2016, Christian Schudoma, The Sainsbury Laboratory"
__credits__ = ["Pirasteh Pahlavan", "Agathe Jouet", "Yogesh Gupta"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Christian Schudoma"
__email__ = "cschu1981@gmail.com"
__status__ = "Development"
