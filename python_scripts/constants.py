#!/usr/bin/env python

import sys, os, time
import shutil, re

thisDir = os.getcwd()

# this the directory of the script loading this one
# so for example, when using wf_reproductive_rom,
# thisDir is the directory where that script is located
shwaveRepoDir = thisDir
print("shwaveRepoDir = {}".format(shwaveRepoDir))

fomExeName = 'shwave_fom'
romExeName = 'shwave_rom'
svdExeName = 'computeThinSVD'
extractStateExeName = 'extractStateFromSnaps'
reconstructFomStateExeName = 'reconstructFomState'
reconstructSeismoExeName = 'reconstructSeismogram'

# where the mesh generator is located
meshExeDir = thisDir + '/../meshing'
meshExe  = 'create_single_mesh.py'