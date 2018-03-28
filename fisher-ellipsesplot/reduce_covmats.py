import numpy as np

transform = {
  'LCDM': ['omega_b', 'omega_cdm', 'A_s', 'h', 'M_tot', 'n_s'],
  'HDM': ['omega_b', 'omega_cdm', 'A_s', 'h', 'w0_fld', 'M_tot', 'N_star', 'n_s'],
  'WCDM': ['omega_b', 'omega_cdm', 'A_s', 'h', 'w0_fld', 'M_tot', 'n_s'],
}

translation = {
  'omega_b': 'omegab',
  'omega_cdm': 'omegac',
  'A_s': '10^9As',
  'h': 'h',
  'M_tot': 'M_nu',
  'n_s': 'ns',
  'w0_fld': 'w0',
  'N_star': 'N*',
}

for m, out_params in transform.items():
  with open('CMs/{}.covmat'.format(m)) as f:
    line = f.readline().split()[1:]
    indices = []
    for param in out_params:
      needle = '{},'.format(param)
      indices.append(line.index(needle))

    data = np.loadtxt('CMs/{}.covmat'.format(m))

    cm = np.take(np.take(data, indices, axis=1), indices, axis=0)

    header_params = [translation[key] for key in out_params]
    as_index = out_params.index('A_s')
    scale = np.identity(len(out_params))
    scale[as_index, as_index] = 1e9
    cm = np.dot(np.dot(scale, cm), scale)

    np.savetxt('CMs/{}_reduced.covmat'.format(m), cm, header='\t'.join(header_params))

