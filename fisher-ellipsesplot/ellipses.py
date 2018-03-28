import glob

from fishertools import FileLoader
from fisherplots import GridPlot


plots = (
    ('HDM', '0s$HDM', 'LCDM', '0s$LCDM'),
    ('HDM', '0s$HDM', 'WCDM', '0s$HDM'),
)

gcfilename = './FMs/GC/FisherMatrix-GC-Euclid-fiducial${model}{fiducial}-{cluster}earPk-*final.txt'
wlfilename = './FMs/WL/FisherMatrix-Euclid-fiducial${model}{fiducial}-{cluster}earPk-*.txt'
outname = 'plots/{run}/{cluster}-{marg}-{model1}-vs-{model2}-fidu1_{fiducial1}-fidu2_{fiducial2}.pdf'


for model1, fiducial1, model2, fiducial2 in plots:
  for cluster in ('nonlin', 'lin'):
    for treat_matrix in ('marg', 'fix'):
      filename1 = gcfilename.format(model=model1, fiducial=fiducial1, cluster=cluster)
      filename2 = gcfilename.format(model=model2, fiducial=fiducial2, cluster=cluster)

      gcfl = FileLoader(
	      glob.glob(filename1) + glob.glob(filename2),
	      ['{} - {}'.format(model1, fiducial1),'{} - {}'.format(model2, fiducial2)]
      )

      # there is only nonlin for WL
      filename1 = wlfilename.format(model=model1, fiducial=fiducial1, cluster='nonlin')
      filename2 = wlfilename.format(model=model2, fiducial=fiducial2, cluster='nonlin')

      wlfl = FileLoader(
	      glob.glob(filename1) + glob.glob(filename2),
	      ['{} - {}'.format(model1, fiducial1),'{} - {}'.format(model2, fiducial2)]
      )

      if treat_matrix == 'marg':
	gcmatrices = gcfl.get_matrices(marginalise_over_distinct=True)
	wlmatrices = wlfl.get_matrices(marginalise_over_distinct=True)
	mcmcmatrices = [gcfl.load_covariance_to_fisher('CMs/{}.covmat'.format(m), 'MCMC', gcmatrices[i].variables) for i,m in enumerate((model1, model2))]
      else:
	gcmatrices = gcfl.get_matrices(maximize_over_distinct=True)
	wlmatrices = wlfl.get_matrices(maximize_over_distinct=True)
	mcmcmatrices = [gcfl.load_covariance_to_fisher('CMs/{}.covmat'.format(m), 'MCMC', gcmatrices[i].variables) for i,m in enumerate((model1, model2))]

      for run in ('GC', 'total', 'total_mcmc'):
	matrices = []
	for gcm, wlm, mcm in zip(gcmatrices, wlmatrices, mcmcmatrices):
	  if run == 'GC':
	    matrices.append(gcm)
	  elif run == 'total':
	    matrices.append(gcm.add(wlm, label=gcm.label))
	  elif run == 'total_mcmc':
	    matrices.append(gcm.add(wlm).add(mcm, label=gcm.label))


	gp = GridPlot(matrices)
	gp.plot()
	gp.save(bbox_inches='tight', name=outname.format(run=run, cluster=cluster, marg=treat_matrix, model1=model1, model2=model2, fiducial1=fiducial1, fiducial2=fiducial2))