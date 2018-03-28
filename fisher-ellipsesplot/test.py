import glob
import copy
import itertools

import numpy as np
import matplotlib.pyplot as plt

from fishertools import FileLoader
from fisherplots import HDMPlot, GridPlot, DiagonalPlot, ContourPlot

class TestException(Exception):
  pass


filename = './FMs/GC/FisherMatrix-GC-Euclid-fiducial${model}0s$HDM-{cluster}earPk-*final.txt'

outname = 'plots/ellipses_consistency_rel.pdf'

linfl = FileLoader(
	glob.glob(filename.format(model='HDM', cluster='lin'))
)

linmatrices = linfl.get_matrices()

if not np.isclose(linmatrices[0].matrix[0,0], 5.11297887199981e6):
  raise TestException('first element does not match')

linmatrices[0].marginalise_all(['omegac','omegab','10^9As','h','w0','M_nu'], update=True)

if not np.isclose(linmatrices[0].matrix[0,0], 1.79239619e+06):
  print linmatrices[0].matrix
  raise TestException('first element does not match after marg')

linmatrices[0].marginalise_all(['omegab','10^9As','h','w0','M_nu'], update=True)

if not np.isclose(linmatrices[0].matrix[0,0], 8.97101411e+06):
  print linmatrices[0].matrix
  raise TestException('first element does not match after marg')

for var in linmatrices[0].variables:
  print var.name, np.sqrt(var.uncertainty)
  if var.name == 'omegab':
    if not np.isclose(var.uncertainty**2, 1.20140974e-05):
      raise TestException("uncertainty is wrong")

ax = plt.gca()
p = DiagonalPlot(ax, linmatrices, '10^9As')
p.plot()
plt.title('Should be gauss around 2.15 with stddev 0.59')
plt.show()

plt.clf()

ax = plt.gca()
p = ContourPlot(ax, linmatrices, '10^9As', 'h')
p.plot()
plt.xlim(1.0, 3.0)
plt.ylim(0.5, 1.0)
plt.show()
