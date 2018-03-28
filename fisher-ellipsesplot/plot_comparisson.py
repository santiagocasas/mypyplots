import glob
import copy
import itertools

import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt

from fishertools import FileLoader, get_common_vars
from fisherplots import HDMPlot

filename = './FMs/GC/FisherMatrix-GC-{exp}-fiducial${model}0s${model}-{lin}earPk-*vfin2-final{transformed}.txt'
wlfilename = './FMs/WL/FisherMatrix-{exp}-fiducial${model}0s${model}-nonlinearPk-*vfin2{transformed}.txt'

outname = 'plots/ellipses_comparisson.pdf'

blue = plt.get_cmap('Blues')
red = plt.get_cmap('Reds')

colors = itertools.cycle([red, blue])

lcdmfl = FileLoader(
	glob.glob(filename.format(exp='DESI-ELG', model='LCDM', transformed='', lin='lin'))
	+ glob.glob(filename.format(exp='Euclid', model='LCDM', transformed='', lin='lin'))
	+ glob.glob(wlfilename.format(exp='Euclid', model='LCDM', transformed=''))
	+ glob.glob(filename.format(exp='SKA2', model='LCDM', transformed='', lin='nonlin'))
	+ glob.glob(wlfilename.format(exp='SKA2', model='LCDM', transformed=''))
)
hdmfl = FileLoader(
	glob.glob(filename.format(exp='DESI-ELG', model='HDM', transformed='-transformed', lin='lin'))
	+ glob.glob(filename.format(exp='Euclid', model='HDM', transformed='-transformed', lin='lin'))
	+ glob.glob(wlfilename.format(exp='Euclid', model='HDM', transformed='-transformed'))
	+ glob.glob(filename.format(exp='SKA2', model='HDM', transformed='-transformed', lin='nonlin'))
	+ glob.glob(wlfilename.format(exp='SKA2', model='HDM', transformed='-transformed'))
)


blue = plt.get_cmap('Blues')
red = plt.get_cmap('Reds')

lcdmmatrices = lcdmfl.get_matrices(colors=itertools.cycle([red, ]))
hdmmatrices = hdmfl.get_matrices(colors=itertools.cycle([blue, ]))

#for m in hdmmatrices:
#  m.fix(['N*'])

lcdmmcmc = lcdmfl.load_covariance_to_fisher('CMs/LCDM.covmat', 'LCDM MCMC', lcdmmatrices[0].variables, color=red)
hdmmcmc = hdmfl.load_covariance_to_fisher('CMs/HDM.covmat', 'HDM MCMC', hdmmatrices[0].variables, color=blue)

lcdmdesi = lcdmmatrices[0]
hdmdesi = hdmmatrices[0]

lcdmeuclidlin = lcdmmatrices[1].add(lcdmmatrices[2])
hdmeuclidlin = hdmmatrices[1].add(hdmmatrices[2])

lcdmska = lcdmmatrices[3].add(lcdmmatrices[4])
hdmska = hdmmatrices[3].add(hdmmatrices[4])

lcdmdesi = lcdmdesi.add(lcdmmcmc)
hdmdesi = hdmdesi.add(hdmmcmc)

lcdmeuclidlinplus = lcdmeuclidlin.add(lcdmmcmc)
hdmeuclidlinplus = hdmeuclidlin.add(hdmmcmc)

lcdmska = lcdmska.add(lcdmmcmc)
hdmska = hdmska.add(hdmmcmc)

common_vars = get_common_vars(lcdmmatrices[0], hdmmatrices[0])

matrices_to_plot = [[lcdmmcmc, hdmmcmc], [lcdmdesi, hdmdesi], [lcdmeuclidlinplus, hdmeuclidlinplus], [lcdmska, hdmska]]

for ml in matrices_to_plot:
  ml[0].line_kwargs['linestyle']='dashed'
  ml[0].patch_kwargs['linestyle']='dashed'
  for m in ml:
    m.marginalise_all(common_vars, update=True, get_covariance=True)

with open('errors_lcdm_hdm.dat', 'w') as f:
  f.write('\t'.join(common_vars))
  f.write('\n')
  for ml in matrices_to_plot:
    for m in ml:
      f.write(m.label + '\n')
      f.write('\t'.join(["{0:.2f}%".format(100.0 * m.get_variable(vname).uncertainty / m.get_variable(vname).fiducial) for vname in common_vars]))
      f.write('\n')

gp = HDMPlot(matrices_to_plot, ['MCMC', '+DESI-like', '+EUCLID-like',  '+SKA2-like'], varlist=['ns', 'omegac', '10^9As'], refvar = 'ns', varstoplot = ['omegac', '10^9As'], zoom={'ns': 1.25, 'omegac': 2.0})
gp.plot()
gp.save(bbox_inches='tight', name=outname, hspace=0.1)
