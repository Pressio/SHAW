#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys, os, time
import numpy as np
from numpy import linspace, meshgrid
from matplotlib import cm
import collections
from argparse import ArgumentParser
import random
import scipy.sparse as sp
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee


# d is a dictionary with the full mesh graph
# returns the graph in sparse matrix format
def convertGraphDicToSparseMatrix(d):
  # vectorize all entries of the graph so that each entry in the new
  # arrays contains (row_index, col_index, 1) refereing to one entry in the matrix

  row_ind = [k for k, v in d.items() for _ in range(len(v))]
  col_ind = []
  for k, v in d.items():
    for i in v: col_ind.append(i)

  #col_ind = [int(i) for ids in d.values() for i in ids]
  val = np.ones(len(row_ind)) # can just put ones, since this is a graph

  # from the graph, create a sparse matrix
  spM = sp.csr_matrix(( val, (row_ind, col_ind)) )
  return spM
