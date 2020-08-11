#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys, os, time
import numpy as np
from numpy import linspace, meshgrid
from matplotlib import cm

# rCmb  ----------------
#       |
#       |
#  r,j  |
#       |
#       |
# rSurf |_______________
#       thL             thR
#             theta, i
#
#
#   o--*--o--*--o--*--o cmb
#   |     |     |     |
#   x     x     x     x
#   |     |     |     |
#   o--*--o--*--o--*--o
#   |     |     |     |
#   x     x     x     x
#   |     |     |     |
#   o--*--o--*--o--*--o surface
#
# r varies from the earth surface to cmb
#
# we use:
# vp  = v_r,phi          symbol = o
# srp = sigma_r,phi      symbol = x
# stp = sigma_theta,phi  symbol = *
