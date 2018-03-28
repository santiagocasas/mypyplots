import numpy as np
from scipy.optimize import brentq

# HDM consistency relations
def xi(w0, Omega):
  xi = 16.0 / 3.0 / (1+w0) * (1.0 / np.sqrt(Omega) - 0.5 * (1.0/Omega - 1.0) * np.log((1.0+np.sqrt(Omega))/(1.0-np.sqrt(Omega))) )**2 - 6.0
  return 1.0 / xi

def n_s(w0, N, Omega):
  x = xi(w0, Omega)
  return 1.0 - 8.0 * x / np.tanh(4*x*N)

def alpha_s(w0, N, Omega):
  x = xi(w0, Omega)
  return -32.0 * x**2 / (np.sinh(4*x*N))**2


# data from Planck XX - Table 3
# or Planck XIII Table 4  TTTEEE+lowP
# OMEGA_B_PLANCK = 0.02225
# OMEGA_C_PLANCK = 0.1198
# THETA_MC_PLANCK = 1.04077
# TAU_PLANCK = 0.079
# LN_10AS_PLANCK = 3.094
# AS_10_9_PLANCK = np.exp(LN_10AS_PLANCK) / 10.0
# N_S_PLANCK = 0.9645
# N_S_PLANCK_SIG = 0.0049
# N_S_PLANCK_1SIG = N_S_PLANCK - 1.0 * N_S_PLANCK_SIG
# N_S_PLANCK_2SIG = N_S_PLANCK - 2.0 * N_S_PLANCK_SIG
# H_PLANCK = 0.6727
# OMEGA_M = 0.3156


# SUM_MASS_NUS = 0.06
# N_EFOLDS_INFLATION = 60

# HDM_WO_CENTRAL = brentq(lambda(w): n_s(w, N_EFOLDS_INFLATION, 1.0 - OMEGA_M) - N_S_PLANCK, -0.9999, -0.9)
# HDM_W0_1SIG = brentq(lambda(w): n_s(w, N_EFOLDS_INFLATION, 1.0 - OMEGA_M) - N_S_PLANCK_1SIG, -0.9999, -0.9)
# HDM_W0_2SIG = brentq(lambda(w): n_s(w, N_EFOLDS_INFLATION, 1.0 - OMEGA_M) - N_S_PLANCK_2SIG, -0.9999, -0.9)

# these are the values from the final run
HDM_FIDUCIAL = {
    'omega_b': 0.022275,
    'omega_cdm': 0.11846,
    'h': 0.67889,
    '10^9A_s': np.exp(3.0661) / 10.0,
    'tau_reio': 0.06714,
    'M_tot': 0.062274,
    'w0_fld': -0.99697,
    'N_star': 60.661,
    'n_s': 0.96544, # this is not exactly on the consistency relation, calculate instead
    'r': 0.0031321
}
HDM_FIDUCIAL_W0_SIG = 0.000887
# HDM_FIDUCIAL_NS_SIG = 0.00246

N_EFOLDS_HDM = 60.661

LCDM_FIDUCIAL = {
    'omega_b': 0.022342,
    'omega_cdm': 0.11753,
    'h': 0.68212,
    '10^9A_s': np.exp(3.0794) / 10.0,
    'tau_reio': 0.075253,
    'n_s': 0.97059,
    'r': 0.031786,
    'M_tot': 0.07922,
}
LCDM_FIDUCIAL_NS_SIG = 0.00485

#print "== cosmology =="
# print "n_s 1 sigma = {}   corresponds to w0={}".format(n_s(HDM_W0_1SIG, N_EFOLDS_INFLATION, 1.0 - OMEGA_M), HDM_W0_1SIG)
# print "n_s 2 sigma = {}   corresponds to w0={}".format(n_s(HDM_W0_2SIG, N_EFOLDS_INFLATION, 1.0 - OMEGA_M), HDM_W0_2SIG)

  
