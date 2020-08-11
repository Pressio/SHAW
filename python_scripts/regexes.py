#!/usr/bin/env python

import numpy as np
import sys, re
from argparse import ArgumentParser

# ------------------------------------------------------------
# regex to use for detecting the output of the cumulative_energy script
# ------------------------------------------------------------
# \d{1,} match one or more (any) digits before the .
# \d{9,} match nine or more (any) digits
cumEnRegExp = re.compile(r'numBasis: \d{1,}')

# ------------------------------------------------------------
# regex to use for detecting the output of the cumulative_energy script
# ------------------------------------------------------------
cumEnRegExpVp = re.compile(r'Nvp: \d{1,}')
cumEnRegExpSp = re.compile(r'Nsp: \d{1,}')
