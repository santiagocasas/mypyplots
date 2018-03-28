from math import degrees
import copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patch
from matplotlib.ticker import StrMethodFormatter
from matplotlib.collections import PatchCollection

class DiagonalPlot(object):
  def __init__(self, ax, fishers, variable):
    self.fishers = fishers
    self.ax = ax
    self.variable = variable

  def gaussian(self, mu, sigma, x):
    return np.exp( - (x - mu)**2 / (2 * sigma**2) )

  def plot(self):
    for f in self.fishers:
      v = f.get_variable(self.variable)
      x = np.linspace(v.range[0], v.range[1], 500)
      y = self.gaussian(v.fiducial, v.uncertainty, x)
      self.ax.plot(x,y, label=f.label, **f.line_kwargs)
      x = np.linspace(v.fiducial - 2.0 * v.uncertainty, v.fiducial + 2.0 * v.uncertainty, 500)
      y = self.gaussian(v.fiducial, v.uncertainty, x)
      self.ax.fill_between(x,y, 0, facecolor=f.line_kwargs['color'], edgecolor=None, alpha=0.08)
      x = np.linspace(v.fiducial - 1.0 * v.uncertainty, v.fiducial + 1.0 * v.uncertainty, 500)
      y = self.gaussian(v.fiducial, v.uncertainty, x)
      self.ax.fill_between(x,y, 0, facecolor=f.line_kwargs['color'], edgecolor=None, alpha=0.04)

class ContourPlot(object):
  ellipses_kwargs = {
      'fill': True,
  }
  sigma_values = [1.52, 2.48]

  def get_patch_kwargs(self, fisher):
    kwargs = copy.copy(self.ellipses_kwargs)
    kwargs.update(fisher.patch_kwargs)
    return kwargs

  def __init__(self, ax, fishers, var1, var2, patch_kwargs = {}):
    self.ax = ax
    self.fishers = fishers
    self.variable1 = var1
    self.variable2 = var2
    self.ellipses_kwargs.update(patch_kwargs)

  def sigma1_ellipse(self, FM_matrix, v1, v2, scal=1.52):
    sub_matrix = FM_matrix.marginalise_all([v1.name,v2.name], get_covariance=True)
    if v1.index > v2.index:
        sub_matrix = np.rot90(sub_matrix, k=2)

    eigvals, eigvects = np.linalg.eig(sub_matrix)

    for j in range(eigvects.shape[1]): # check columns
        if (eigvects[0,j]<0 and eigvects[1,j]<0) or (eigvects[0,j]>0 and eigvects[1,j]<0):
            eigvects[:,j] = -eigvects[:,j]
    u, v = eigvects[:,0], np.array([1,0])
    c = np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v) # -> cosine of the angle between
    angle = np.arccos(np.clip(c, -1, 1))

    a, b = [scal*np.abs(np.sqrt(eigval)) for eigval in eigvals]

    if v1.index > v2.index:
        return 2*a, 2*b, np.pi/2. - angle
    else:
        return 2*b, 2*a, -(np.pi/2. - angle)

  def plot(self):
    ellipses = []
    for f in self.fishers:
      var1 = f.get_variable(self.variable1)
      var2 = f.get_variable(self.variable2)
      for s in self.sigma_values:
	results = self.sigma1_ellipse(f,var1, var2, scal=s)
	ellipses.append(patch.Ellipse(xy=(var1.fiducial, var2.fiducial), width=results[0], height=results[1], angle=degrees(results[2]), **self.get_patch_kwargs(f)))
    coll = PatchCollection(ellipses, match_original=True)
    self.ax.add_collection(coll)


class GridPlot(object):
  tickformat = "{x:.3f}"
  nr_ticks = 4
  tick_kwargs = {
    'labelsize': 5,
  }
  xtick_label_kwargs = {
    'rotation': 45,
  }
  ytick_label_kwargs = {

  }
  legend_kwargs = {
      'fontsize': 12,
  }

  def __init__(self, fishers, varlist=None, tick_kwargs={}, tick_label_kwargs={}, xtick_label_kwargs={}, ytick_label_kwargs={}, zoom={}, dim=None):
    self.fishers = fishers
    self.dim = fishers[0].matrix.shape
    self.tick_kwargs.update(tick_kwargs)
    xtick_label_kwargs.update(tick_label_kwargs)
    ytick_label_kwargs.update(tick_label_kwargs)
    self.xtick_label_kwargs.update(xtick_label_kwargs)
    self.ytick_label_kwargs.update(ytick_label_kwargs)
    self.zoom = zoom

    # TODO: Need to check if fishers have same vars
    for obj in fishers:
      f = obj.matrix
      if f.ndim != 2 or f.shape[0] != f.shape[1]:
	raise Exception("fishers needs to be a list of quadratic matrices")
      if not f.shape == self.dim:
	raise Exception("all fishers need to have the same dimension")

    if varlist is None:
      varlist = [var.name for var in fishers[0].variables]
    self.nvars = len(varlist)
    self.varlist = varlist

    dim = (self.nvars, self.nvars) if dim is None else dim

    self.figure, self.axarr = plt.subplots(*dim)

  def get_labels(self):
    labels = {}
    for name in self.varlist:
      labels[name] = self.fishers[0].get_variable(name).label
    return labels

  def get_ranges(self):
    uncertainties = {}
    ranges = {}
    # compute ranges
    for name in self.varlist:
      for f in self.fishers:
	var = f.get_variable(name)
	if var:
	  # just set the values to something, will be overwriten 5 lines below
	  if name not in uncertainties:
	    uncertainties[name] = var.uncertainty
	  if name not in ranges:
	    ranges[name] = [var.fiducial, var.fiducial]
	  # find the actual ranges
	  if var.uncertainty > uncertainties[name]:
	    uncertainties[name] = var.uncertainty
	  low_lim = var.range[0]
	  up_lim = var.range[1]
	  if low_lim < ranges[name][0]:
	    ranges[name][0] = low_lim
	  if up_lim > ranges[name][1]:
	    ranges[name][1] = up_lim
    # apply the zoom
    for name in self.varlist:
      fid = 0.5 * (ranges[name][0] + ranges[name][1])
      dif = fid - ranges[name][0]
      z = self.zoom[name] if name in self.zoom else 1.0
      ranges[name] = [fid - dif / z, fid + dif / z]
    # compute ticks
    ticks = {}
    for name in self.varlist:
      width = (ranges[name][1] - ranges[name][0]) / (self.nr_ticks+1)
      ticks[name] = np.linspace(ranges[name][0] + 0.5 * width, ranges[name][1] - 0.5 * width, 5)
    return ranges, ticks

  def plot(self):
    ranges, ticks = self.get_ranges()
    labels = self.get_labels()

    # diagonal plots
    for i, name in enumerate(self.varlist):
      ax = self.axarr[i][i]
      p = DiagonalPlot(ax, self.fishers, name)
      p.plot()
      ax.set_xlim(*ranges[name])
      ax.set_xticks(ticks[name])
      ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
      ax.set_xticklabels([])
      ax.set_yticklabels([])
      ax.tick_params(**self.tick_kwargs)

    # bottom left plots
    for i, i_name in enumerate(self.varlist):
      for j, j_name in enumerate(self.varlist):
	if j >= i:
	  continue
	# make a plot in bottom right
	ax = self.axarr[i][j]
	p = ContourPlot(ax, self.fishers, j_name, i_name)
	p.plot()
	ax.set_xlim(*ranges[j_name])
	ax.set_ylim(*ranges[i_name])
	ax.set_xticks(ticks[j_name])
	ax.set_yticks(ticks[i_name])
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.tick_params(**self.tick_kwargs)
	# turn of the corresponding plot in the top right corner
	self.axarr[j][i].axis('off')

    # restore the labels
    n = self.nvars - 1
    formater = StrMethodFormatter(self.tickformat)
    for i, name in enumerate(self.varlist):
      self.axarr[i][0].set_yticklabels([formater.format_data_short(k) for k in self.axarr[i][0].get_yticks()], **self.ytick_label_kwargs)
      self.axarr[i][0].set_ylabel('${}$'.format(labels[name]))
    for i, name in enumerate(self.varlist):
      self.axarr[n][i].set_xticklabels([formater.format_data_short(k) for k in self.axarr[n][i].get_xticks()], **self.xtick_label_kwargs)
      self.axarr[n][i].set_xlabel('${}$'.format(labels[name]))

    handles, labels = self.axarr[0][0].get_legend_handles_labels()
    self.axarr[0][n].legend(handles, labels, **self.legend_kwargs)

  def save(self, name='plot.pdf', wspace=0.0, hspace=0.0, **kwargs):
    plt.subplots_adjust(wspace=wspace, hspace=hspace)
    plt.savefig(name, **kwargs)


class HDMPlot(GridPlot):
  left_legend_kwargs = {
      'fontsize': 10,
  }
  def __init__(self, fishers, caselabels, varlist, refvar, varstoplot, **kwargs):
    self.cases = fishers
    self.caselabels = caselabels
    self.refvar = refvar
    self.varstoplot = varstoplot
    matrices = [m for matrix_set in fishers for m in matrix_set]
    dim = (len(fishers), len(varlist))
    super(HDMPlot, self).__init__(matrices, varlist=varlist, dim=dim, **kwargs)

  def get_ranges(self, *args):
    ranges, ticks = super(HDMPlot, self).get_ranges(*args)
    if 'w0' in ranges:
      ranges['w0'] = [-1.001, -0.968]
      ticks['w0'] = [-1, -0.99, -0.98, -0.97]
    return ranges, ticks

  def plot(self):
    ranges, ticks = self.get_ranges()
    labels = self.get_labels()

    # diagonal plots
    for i, case in enumerate(self.cases):
      ax = self.axarr[i][0]
      ax.set_xlim(*ranges[self.refvar])
      ax.set_xticks(ticks[self.refvar])
      ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
      ax.set_xticklabels([])
      ax.set_yticklabels([])
      ax.tick_params(**self.tick_kwargs)
      p = DiagonalPlot(ax, case, self.refvar)
      p.plot()

    for i, case in enumerate(self.cases):
      for j, j_name in enumerate(self.varstoplot):
	ax = self.axarr[i][j+1]
	ax.set_xlim(*ranges[j_name])
	ax.set_ylim(*ranges[self.refvar])
	ax.set_xticks(ticks[j_name])
	ax.set_yticks(ticks[self.refvar])
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.tick_params(**self.tick_kwargs)
	p = ContourPlot(ax, case, j_name, self.refvar)
	p.plot()

    # restore the labels
    n = self.nvars - 1
    formater = StrMethodFormatter(self.tickformat)
    for i, case, caselabel in zip(range(len(self.cases)), self.cases, self.caselabels):
      self.axarr[i][0].set_yticklabels([formater.format_data_short(k) for k in self.axarr[i][0].get_yticks()], **self.ytick_label_kwargs)
      self.axarr[i][0].set_ylabel(caselabel, **self.left_legend_kwargs)
      self.axarr[i][n].set_yticklabels([formater.format_data_short(k) for k in self.axarr[i][n].get_yticks()])
      self.axarr[i][n].tick_params(axis='y', labelleft=False, labelright=True)
      self.axarr[i][n].set_ylabel('${}$'.format(labels[self.refvar]), **self.legend_kwargs)
      self.axarr[i][n].get_yaxis().set_label_position("right")
    for j, name in enumerate([self.refvar, ] + self.varstoplot):
      self.axarr[0][j].tick_params(axis='x', labeltop=True, labelbottom=False)
      self.axarr[0][j].set_xticklabels([formater.format_data_short(k) for k in self.axarr[i][j].get_xticks()], **self.xtick_label_kwargs)
      self.axarr[len(self.cases) - 1][j].set_xticklabels([formater.format_data_short(k) for k in self.axarr[i][j].get_xticks()], **self.xtick_label_kwargs)
      self.axarr[len(self.cases) - 1][j].set_xlabel('${}$'.format(labels[name]))