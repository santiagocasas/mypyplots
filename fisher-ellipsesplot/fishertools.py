import copy
import itertools

import matplotlib.pyplot as plt
import numpy as np

class Variable:
    latex_labels = {
      '10e9As': '10^9A_s',
      '10^9As': '10^9A_s',
      'omegab': '\Omega_{b} h^2',
      'Omegab': '\Omega_{b}',
      'omegac': '\Omega_{cdm} h^2',
      'Omegac': '\Omega_{cdm}',
      'Omegad': '\Omega_{d}',
      'w0': 'w_0',
      'wa': 'w_a',
      'ns': 'n_s',
      'h': 'h',
      'omegam': '\Omega_{m} h^2',
      'Omegam': '\Omega_{m}',
      'sigma8': '\sigma_{8}',
      'gamma': '\gamma',
      'M_nu' : r'\Sigma m_{\nu}'
    }
    scales = {
      'omegab': 100.0,
      'Omegab': 100.0,

    }
    def __init__(self, name, index, fiducial, scaling=None):
        self.name = name
        self.index = index
        self.fiducial = fiducial
        self.uncertainty = 0.0
        self.range = [0,0]
        self.label = self.latex_labels.get(name, name)
        if scaling:
	  self.scaling = scaling
	else:
	  self.scaling = self.scales.get(name, 1.0)

def get_color_kwargs(color):
  alpha = list(color(0.8))
  no_alpha = list(color(0.5))
  alpha[3] = 0.8
  no_alpha[3] = 0.2
  return {
    'facecolor': tuple(no_alpha),
    'edgecolor': tuple(alpha), #tuple(no_alpha),
  }

def get_common_vars(f1, f2):
  vars1 = f1.varlist()
  vars2 = f2.varlist()
  common_vars = set.intersection(set(vars1), set(vars2))
  return list(common_vars)

class FileLoader:
  colors = itertools.cycle([plt.get_cmap('Blues'), plt.get_cmap('Blues'), plt.get_cmap('Blues')])

  def __init__(self, filenames, legendnames=None):
    self.filenames = filenames
    self.legendnames = legendnames if legendnames else filenames

  def load_vars_from_file(self, filename, filt=None):
    with open(filename, 'r') as f:
      varnames = f.readline().split()[1:]
      fiducials = f.readline().split()[1:]

      f.close()

      varlist = []

      for index, name, fid in zip(range(len(varnames)), varnames, fiducials):
	varlist.append(Variable(name, index, float(fid)))

      if filt:
	sort_varlist = []
	for name in filt:
	  v = None
	  for v1 in varlist:
	    if v1.name == name:
	      v = v1
	  if v:
	    sort_varlist.append(v)
	  else:
	    raise Exception('Unknown variable {}'.format(name))
	return sort_varlist
      else:
	return varlist

  def get_matrices(self, marginalise_over_distinct=False, maximize_over_distinct=False, colors=None):
    matrices = []
    varlists = []

    for fname, lname in zip(self.filenames, self.legendnames):
      d = np.loadtxt(fname)
      print fname
      varlist = self.load_vars_from_file(fname)
      varlists.append(set([v.name for v in varlist]))
      f = FisherMatrix(d, varlist, lname)
      f.compute_ranges_and_uncertainties()
      matrices.append(f)

      colors = self.colors if not colors else colors

      c = next(colors)
      f.set_patch_kwargs(get_color_kwargs(c))
      f.set_line_kwargs({'color': c(0.8)})

    common_vars = set.intersection(*varlists)
    common_vars = list(common_vars)
    self.common_vars = common_vars

    if marginalise_over_distinct or maximize_over_distinct:
      # establish a list of common vars
      for m in matrices:
          if marginalise_over_distinct:
              m.marginalise_all(common_vars, update=True)
          else:
              m.maximize_all(common_vars, update=True)
    return matrices

  def get_indices_for_varnames(self, f, names):
    with open(f, 'r') as fi:
        line = fi.readline().split()[1:]
        indices = []
        variables = []
        for var in names:
            for ind, s in enumerate(line):
                if var == s[:-1]:
                    indices.append(ind)
        if len(indices) != len(names):
	    print names
	    print f
            raise Exception('could not find all vars')
        return indices

  def translate_to_montepython(self, var):
    if var in ['ns', ]:
      return 'n_s'
    elif var in ['10^9As', ]:
      return 'A_s'
    elif var in ['omegac', ]:
      return 'omega_cdm'
    elif var in ['omegab', ]:
      return 'omega_b'
    elif var in ['M_nu', ]:
      return 'M_tot'
    elif var in ['w0', ]:
      return 'w0_fld'
    elif var in ['N*', ]:
      return 'N_star'
    else:
      return var

  def load_covariance_to_fisher(self, filename, label, variable_list, is_montepython=True, color=None):
    varnames = []
    for var in variable_list:
      varnames.append(self.translate_to_montepython(var.name))
    indices = self.get_indices_for_varnames(filename, varnames)
    m = np.loadtxt(filename)
    m = m[:, indices][indices, :]

    if is_montepython:
      as_index = varnames.index('A_s')
      scale = np.identity(len(indices))
      scale[as_index, as_index] = 1e9
      m = np.dot(np.dot(scale, m), scale)

    m = np.linalg.inv(m)

    f = FisherMatrix(m, copy.deepcopy(variable_list), label)

    f.compute_ranges_and_uncertainties()

    c = color if color else next(self.colors)

    f.set_patch_kwargs(get_color_kwargs(c))
    f.set_line_kwargs({'color': c(0.8)})
    return f



class FisherMatrix:
  line_kwargs = {}
  patch_kwargs = {}
  sigma_limit = 3.0

  def __init__(self, matrix, variables, label):
    self.matrix = matrix
    self.variables = variables
    self.label = label

  def __str__(self):
    return str(self.matrix)

  def varlist(self):
    return [var.name for var in self.variables]

  def compute_ranges_and_uncertainties(self):
    for v in self.variables:
      v.uncertainty = self.uncertainty(v.name)
      v.range = (v.fiducial - self.sigma_limit * v.uncertainty, v.fiducial + self.sigma_limit * v.uncertainty)

  def set_patch_kwargs(self, patch_kwargs):
    self.patch_kwargs = patch_kwargs

  def set_line_kwargs(self, line_kwargs):
    self.line_kwargs = line_kwargs

  def get_variable(self, varname):
    for v in self.variables:
      if v.name == varname:
	return v

  def get_index(self, varname):
    for v in self.variables:
      if v.name == varname:
	return v.index

  def add(self, matrix, label=None):
    for var1, var2 in zip(self.variables, matrix.variables):
      if var1.name != var2.name or var1.fiducial != var2.fiducial:
	print var1.name, var1.fiducial
	print var2.name, var2.fiducial
	raise Exception('Matrices need to have the same variables and fiducials')
    data = copy.copy(self.matrix) + copy.copy(matrix.matrix)
    if label is None:
      label = '{} + {}'.format(self.label, matrix.label)
    f = FisherMatrix(data, copy.deepcopy(self.variables), label)
    f.compute_ranges_and_uncertainties()
    f.set_patch_kwargs(self.patch_kwargs)
    f.set_line_kwargs(self.line_kwargs)
    return f

  def fix(self, varnames):
    num_vars = [self.get_index(n) for n in varnames]
    if len(num_vars) != len(varnames) or any(var is None for var in num_vars):
      print varnames, num_vars
      raise 'Unable to find all variables'
    m = np.delete(self.matrix,num_vars,axis=1)
    self.matrix = np.delete(m,num_vars,axis=0)
    for n in varnames:
      ind = None
      for i, var in enumerate(self.variables):
	if n == var.name:
	  ind = i
      del self.variables[ind]

  def marginalise(self, varnames):
    raise Exception('Not implemented')
    num_vars = []
    return self.marginalise_all()

  # takes the variables TO KEEP as an argument
  def marginalise_all(self, varnames, update=False, get_covariance=False):
    num_vars = [self.get_index(n) for n in varnames]
    m = np.linalg.inv(self.matrix)
    m = np.take(m,indices=num_vars,axis=1)
    m = np.take(m,indices=num_vars,axis=0)
    if update:
      new_vars = []
      for index, name in enumerate(varnames):
	var = self.get_variable(name)
	var.index = index
	new_vars.append(var)
      self.variables = new_vars
      self.matrix = np.linalg.inv(m)
    if get_covariance:
      return m
    else:
      return np.linalg.inv(m)

  def maximize_all(self, varnames, update=False):
    num_vars = [self.get_index(n) for n in varnames]
    m = self.matrix
    m = np.take(m,indices=num_vars,axis=1)
    m = np.take(m,indices=num_vars,axis=0)
    if update:
      new_vars = []
      for index, name in enumerate(varnames):
	var = self.get_variable(name)
	var.index = index
	new_vars.append(var)
      self.variables = new_vars
      self.matrix = m
    return m

  def uncertainty(self, varname):
    i = self.get_index(varname)
    m = np.linalg.inv(self.matrix)
    m = np.take(m,indices=i,axis=1)
    m = np.take(m,indices=i,axis=0)
    return np.sqrt(m)