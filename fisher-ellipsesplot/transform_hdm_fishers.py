import glob

import numpy as np

from parameters import *

jac = np.loadtxt('jaco-w-ns-hdm-0s.txt')

hdm_files = glob.glob('./FMs/GC/FisherMatrix-GC-Euclid-fiducial$HDM0s$HDM*.txt') + glob.glob('./FMs/WL/FisherMatrix-Euclid-fiducial$HDM0s$HDM*.txt') \
	  + glob.glob('./FMs/GC/FisherMatrix-GC-SKA2-fiducial$HDM0s$HDM*.txt') + glob.glob('./FMs/WL/FisherMatrix-SKA2-fiducial$HDM0s$HDM*.txt')

W0_INDEX = 4

for f in hdm_files:
  if 'transformed' in f:
    continue
  with open(f, 'r') as fi:
    print f
    titelline = fi.readline().split()[1:]
    fiducialline = fi.readline().split()[1:]

    titelline[W0_INDEX] = 'ns'
    fiducialline[W0_INDEX] = HDM_FIDUCIAL['n_s']

    fisher = np.loadtxt(f)
    fisher_new = np.dot(jac.T, np.dot(fisher, jac))

    newname = f.replace('.txt', '-transformed.txt')

    np.savetxt(newname, fisher_new, header='{}\n{}'.format('\t'.join(titelline), '\t'.join([str(f) for f in fiducialline])))