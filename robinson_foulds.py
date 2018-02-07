#! /usr/bin/env python
#coding: utf-8

"""
compute RF distance between pruned unfiltered tree and reference tree (for all files in the argument directory) 
"""

import glob
import re
import sys
from ete3 import Tree

def prune_and_compute_rf_distance(treeFilepath, referenceFilepath):
    tree = Tree(treeFilepath)
    ref = Tree(referenceFilepath)
    
    tree.prune(ref.get_leaf_names())
    return tree.robinson_foulds(ref)[0]

if __name__=="__main__":
    argvs = sys.argv
    argc = len(argvs)
    if (argc != 2):
        print('Usage: # python {} directoryPath'.format(argvs[0]))
        quit()
    
    refFilepath_lst = glob.glob(argvs[1] + "problem*.reference.nwk")
    for referenceFilepath in refFilepath_lst:
        treeFilepath = referenceFilepath.replace("reference", "prank.unfilteredPhyml")
        match = re.search(r'problem[0-9]*', treeFilepath)
        print(match.group(0), ":", prune_and_compute_rf_distance(treeFilepath, referenceFilepath))
    