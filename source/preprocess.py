#! /usr/bin/env python
#coding: utf-8

"""
save np file of vectorized alignment (simply char->int, 0~5) and output(to std) the id of newick topology
./preprocess.py [alignmentFilepath] [newickFilepath]

TODO: directorypath as input (or similar), output format
"""

import itertools
import numpy as np
import sys

from Bio import SeqIO
from ete3 import Tree
from sklearn.ensemble import RandomForestClassifier

def alignment_preprcess(alignmentFilepath):
    def sequence_to_vector(sequence):
        sequenceVector = np.array(sequence)
        sequenceVector[np.where(sequenceVector == "_")] = 0
        sequenceVector[np.where(sequenceVector == "N")] = 1
        sequenceVector[np.where(sequenceVector == "A")] = 2
        sequenceVector[np.where(sequenceVector == "C")] = 3
        sequenceVector[np.where(sequenceVector == "T")] = 4
        sequenceVector[np.where(sequenceVector == "G")] = 5
        sequenceVector = sequenceVector.astype(int)
        return sequenceVector

    id_lst = []
    sequenceVector_lst = []
    for record in SeqIO.parse(alignmentFilepath, "fasta"):
        if "enriched" in record.id:
            break
        id_lst.append(record.id)
        sequenceVector_lst.append(sequence_to_vector(record.seq))
    return id_lst, sequenceVector_lst

def tree_preprocess(newickFilepath, id_lst):
    t = Tree(newickFilepath)
    for i, id in enumerate(id_lst):
        t.search_nodes(name=id)[0].name = "{0}".format(i)

    for i, nodeOrder in enumerate(list(itertools.permutations(("0", "1", "2", "3", "4", "5")))):
        tree_nwk = "({0[0]},({0[1]},({0[2]},({0[3]},({0[4]},{0[5]})))));".format(nodeOrder)
        tree = Tree(tree_nwk)
        if tree.get_topology_id() == t.get_topology_id():
            topologyId = int(i/2)
            break
    return topologyId

if __name__=="__main__":
    argvs = sys.argv
    argc = len(argvs)
    if (argc != 3):
        print('Usage: # python {} alignmentFilepath newickFilepath'.format(argvs[0]))
        quit()
        
    alignmentFilepath = argvs[1]
    newickFilepath = argvs[2]
    
    ids, sequenceVector_lst = alignment_preprcess(alignmentFilepath)
    np.save("alignment_vector.npy", np.array(sequenceVector_lst))
    print(tree_preprocess(newickFilepath, ids))