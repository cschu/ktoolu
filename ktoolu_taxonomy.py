#!/usr/bin/env python
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
