""" This script will create sets of time series, phase and Fourier plots for
a given realization of the DDP
Nicolas Grisouard, University of Toronto, August 2021
"""
import ddp_routines as ddp
from scipy.constants import pi
import numpy as np
import matplotlib.pyplot as plt

ftsz = 12  # font size
m = 1.0  # [kg] mass
g = 9.8  # [m/s] gravity
ell = 1.0  # [m] pendulum length
omega0 = (g/ell)**.5

omegad = 2*omega0/3  # driving frequency
Td = 2*pi/omegad  # driving period
gamma = omega0/4  # damping
ntTd = 500  # number of time steps per driving period
nbetas = 400  # number of values of beta to compute for each figure

dt = Td/ntTd  # time step
num_cycles = 40  # number of dirving periods
t_end = num_cycles*Td  # final time
nt = ntTd*num_cycles  # number of time steps
time = np.arange(0., t_end, dt)  # initialize the time array

theta0, dottheta0 = -pi/2, 0.

beta = 1.07  # have to run a new one
theta, dottheta = ddp.generate_time_series(theta0, dottheta0, omegad, omega0,
                                           beta, gamma, time)

ddp.plot_TS(theta, dottheta, omegad, omega0, beta,
            gamma, time/Td, ftsz, tmin=24.)

ddp.plot_spectrum(theta, omegad, omega0, beta, gamma, time, 24., ftsz)

theta_wrapped = ddp.wrap_theta(theta)
ddp.plot_phase(theta_wrapped, dottheta, omegad, omega0, beta, gamma, time,
               ftsz, nconts=32)
