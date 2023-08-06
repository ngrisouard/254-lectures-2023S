""" This script will create bifurcation plots for L12 of PHY254
Nicolas Grisouard, University of Toronto, August 2021
"""
import ddp_routines as ddp
from scipy.constants import pi

ftsz = 12  # font size
m = 1.0  # [kg] mass
g = 9.8  # [m/s] gravity
ell = 1.0  # [m] pendulum length
omega0 = (g/ell)**.5

omegad = 2*omega0/3  # driving frequency
Td = 2*pi/omegad  # driving period
gamma = omega0/4  # damping
ntTd = 500  # number of time steps per driving period
nbetas = 500  # number of values of beta to compute for each figure

theta0, dottheta0 = -pi/2, 0.  # Note: it is now useless

# beta_min, beta_max = 0.9, 1.5
# ddp.plot_bifurcation(beta_min, beta_max, nbetas, omegad, omega0, gamma,
#                      theta0, dottheta0, ntTd, ftsz)

# beta_min, beta_max = 1.06, 1.086
# ddp.plot_bifurcation(beta_min, beta_max, nbetas, omegad, omega0, gamma,
#                      theta0, dottheta0, ntTd, ftsz)

beta_min, beta_max = 1.076, 1.083
ddp.plot_bifurcation(beta_min, beta_max, nbetas, omegad, omega0, gamma,
                     theta0, dottheta0, ntTd, ftsz)
