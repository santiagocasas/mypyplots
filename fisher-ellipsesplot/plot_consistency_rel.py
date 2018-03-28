import glob
import copy
import itertools

import numpy.linalg

import matplotlib.pyplot as plt

from fishertools import FileLoader, get_common_vars
from fisherplots import HDMPlot

gcfilename = './FMs/GC/FisherMatrix-GC-{exp}-fiducial${model}0s$HDM-{lin}earPk-*-vfin2-final.txt'
wlfilename = './FMs/WL/FisherMatrix-{exp}-fiducial${model}0s$HDM-nonlinearPk-*vfin2.txt'

outname = 'plots/ellipses_consistency_rel.pdf'

hdmfl = FileLoader(
	glob.glob(gcfilename.format(exp='DESI-ELG', model='HDM', lin='lin'))
	+ glob.glob(gcfilename.format(exp='Euclid', model='HDM', lin='lin'))
	+ glob.glob(wlfilename.format(exp='Euclid', model='HDM'))
	+ glob.glob(gcfilename.format(exp='SKA2', model='HDM', lin='nonlin'))
	+ glob.glob(wlfilename.format(exp='SKA2', model='HDM'))
)
wcdmfl = FileLoader(
	glob.glob(gcfilename.format(exp='DESI-ELG', model='WCDM', lin='lin'))
	+ glob.glob(gcfilename.format(exp='Euclid', model='WCDM', lin='lin'))
	+ glob.glob(wlfilename.format(exp='Euclid', model='WCDM'))
	+ glob.glob(gcfilename.format(exp='SKA2', model='WCDM', lin='nonlin'))
	+ glob.glob(wlfilename.format(exp='SKA2', model='WCDM'))
)

blue = plt.get_cmap('Blues')
green = plt.get_cmap('Greens')

wcdmmatrices = wcdmfl.get_matrices(colors=itertools.cycle([green, ]))
hdmmatrices = hdmfl.get_matrices(colors=itertools.cycle([blue, ]))

#for m in hdmmatrices:
#  m.fix(['N*'])

wcdmmcmc = wcdmfl.load_covariance_to_fisher('CMs/WCDM.covmat', 'WCDM MCMC', wcdmmatrices[0].variables, color=green)
hdmmcmc = hdmfl.load_covariance_to_fisher('CMs/HDM.covmat', 'HDM MCMC', hdmmatrices[0].variables, color=blue)

wcdmdesi = wcdmmatrices[0]
hdmdesi = hdmmatrices[0]

wcdmdesi = wcdmdesi.add(wcdmmcmc)
hdmdesi = hdmdesi.add(hdmmcmc)

wcdmeuclidlin = wcdmmatrices[1].add(wcdmmatrices[2])
hdmeuclidlin = hdmmatrices[1].add(hdmmatrices[2])

wcdmeuclidlinplus = wcdmeuclidlin.add(wcdmmcmc)
hdmeuclidlinplus = hdmeuclidlin.add(hdmmcmc)

wcdmska = wcdmmatrices[3].add(wcdmmatrices[4])
hdmska = hdmmatrices[3].add(hdmmatrices[4])

wcdmska = wcdmska.add(wcdmmcmc)
hdmska = hdmska.add(hdmmcmc)

common_vars = get_common_vars(wcdmmatrices[0], hdmmatrices[0])

matrices_to_plot = [[wcdmmcmc, hdmmcmc], [wcdmdesi, hdmdesi], [wcdmeuclidlinplus, hdmeuclidlinplus], [wcdmska, hdmska]]

for ml in matrices_to_plot:
  ml[0].line_kwargs['linestyle']='dotted'
  ml[0].patch_kwargs['linestyle']='dotted'
  for m in ml:
    m.marginalise_all(common_vars, update=True, get_covariance=True)

with open('errors_wcdm_hdm.dat', 'w') as f:
  f.write('\t'.join(common_vars))
  f.write('\n')
  for ml in matrices_to_plot:
    for m in ml:
      f.write(m.label + '\n')
      f.write('\t'.join(["{0:.2f}%".format(abs(100.0 * m.get_variable(vname).uncertainty / m.get_variable(vname).fiducial)) for vname in common_vars]))
      f.write('\n')

gp = HDMPlot(matrices_to_plot, ['MCMC', '+DESI-like', '+EUCLID-like', '+SKA2-like'], varlist=['w0', 'omegac', '10^9As'], refvar = 'w0', varstoplot = ['omegac', '10^9As'], zoom={'omegac': 2.0})
for i in range(gp.axarr.shape[0]):
  gp.axarr[i][0].axvline(-1,linestyle='-', linewidth=0.5, color=plt.get_cmap('Reds')(0.7))
  for j in range(1, gp.axarr.shape[1]):
    gp.axarr[i][j].axhline(-1,linestyle='-', linewidth=0.5, color=plt.get_cmap('Reds')(0.7))
gp.plot()
gp.save(bbox_inches='tight', name=outname, hspace=0.1)
