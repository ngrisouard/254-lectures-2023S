"""
Routines to support the damped-driven pendulum (DDP) simulations
Nicolas Grisouard, University of Toronto, August 2021
"""
import numpy as np
from scipy.constants import pi
import matplotlib.pyplot as plt


def ddtheta(th, dth, wd, w0, g, t, b):
    """ returns theta acceleration accounting for natural frequency w0,
    forcing frequency w and amplitude b, damping g, system state th and dth """
    return -w0**2 * np.sin(th) - 2*g*dth + b*w0**2*np.cos(wd*t)


def generate_time_series(th0, dth0, wd, w0, b, g, t):
    """ Generate theta and dtheta/dt with the Euler-Cromer method
    th0 = initial angle,
    dth0 = initial angular velocity,
    wd = driving frequency
    w0 = natural oscillation frequency
    t = time array
    b = driving amplitude
    g = damping coefficient """
    dt = t[1] - t[0]  # time step
    th = 0*t  # initialize theta array
    th[0] = th0  # initial angle
    dth = 0*t  # initialize dtheta/dt
    dth[0] = dth0  # initial angular velocity
    for i in range(0, len(t)-1):
        # First, update dtheta/dt
        dth[i+1] = dth[i] + dt*ddtheta(th[i], dth[i], wd, w0, g, t[i], b)
        # Second, use updated dtheta/dt to update theta
        th[i+1] = th[i] + dt*dth[i+1]

    return th, dth


def plot_TS(th, dth, wd, w0, b, g, t, ftsz, tmin):
    """plot time series of th. t must be time/Td, ftsz is font size,
    tmin is the minimum (might need to ged rid of initial transient) """
    t_plot = t[t >= tmin]  # new array, only for times t>=tmin
    th_plot = th[t >= tmin]/pi  # idem for theta; plus, we scale by pi
    # idem for dtheta/dt; plus, we scale by wd*pi
    dth_plot = dth[t >= tmin]/pi/wd

    plt.subplot(2, 1, 1)
    plt.plot(t_plot, th_plot)
    plt.grid()
    plt.xlabel(r"$t/T_d$", fontsize=ftsz)
    plt.ylabel(r"$\theta/\pi$", fontsize=ftsz)
    plt.xlim(t_plot[0], t_plot[-1])

    plt.subplot(2, 1, 2)
    plt.plot(t_plot, dth_plot)
    plt.grid()
    plt.xlabel(r"$t/T_d$", fontsize=ftsz)
    plt.ylabel(r"$\dot\theta/(\pi\omega_d)$", fontsize=ftsz)
    plt.xlim(t_plot[0], t_plot[-1])

    plt.suptitle(r"Time series for " +
                 r"$\omega_0/\omega_d = {0:.1f}$, ".format(w0/wd) +
                 r"$Q = {0:.0f}$, $\beta = {1:.3f}$ rad".format(w0/g, b),
                 fontsize=ftsz)
    plt.tight_layout()
    plt.show()
    return


def plot_phase(th, dth, wd, w0, b, g, t, ftsz, nconts):
    """ draw phase plot; th is theta, dth is dtheta/dt, nconts is number of
    contours for U """
    # First, we plot the potential energy in the background
    thm = min(min(th)/pi, -1)  # lower x limit;
    thM = max(max(th)/pi, +1)  # upper x-limit
    dthM = 1.1*max(max(abs(dth))/pi/wd, +1)  # y-limit; symmetric around zero
    th_lin = np.linspace(thm, thM, 512)
    dth_lin = np.linspace(-dthM, dthM, 64)
    thg, dthg = np.meshgrid(th_lin*pi, dth_lin*pi*wd)
    # not the right units; only want the shape
    Ug = 0.5 * dthg**2 + w0**2 * (1-np.cos(thg))
    # alpha controls transparency
    plt.contourf(thg/pi, dthg/pi/wd, Ug, nconts, alpha=0.1)
    plt.contour(thg/pi, dthg/pi/wd, Ug, [2*w0**2],  # Separatrix
                colors='r', linestyles='--', linewidths=0.5)
    
    # Then we plot the phase trajectory
    # plt.plot(th/pi, dth/pi/wd, 'r.', markersize=2)
    # c=t means color-coded by t; s is size
    plt.scatter(th/pi, dth/pi/wd, c=t, cmap='copper_r', s=1.0)
    plt.ylabel(r'$\dot\theta/(\pi\omega_d)$ [rad/s]')
    plt.xlabel(r'$\theta/\pi$')
    plt.title(r"$\phi$-plot for $\omega_0/\omega_d = {0:.1f}$, ".format(w0/wd)
              + r"$Q = {0:.0f}$, $\beta = {1:.3f}$ rad".format(w0/g, b),
              fontsize=ftsz)
    # plt.axis(xmin=-1, xmax=1)
    plt.tight_layout()
    plt.show()
    return


def plot_spectrum(th, wd, w0, b, g, t, tmin, ftsz):
    """ This one plots the spectrum. Best results will be achieved if
    tmin is as small as possible, but large enough that the transient
    is excluded. It should also be an integer number of periods for the
    function to be periodic. Otherwise, big wiggles. """
    from numpy.fft import rfft, rfftfreq  # the FFT for real entries

    th_for_Fourier = th[t >= tmin]
    An = 2*rfft(th_for_Fourier, norm='forward')
    # To make the Fourier plot interpretable, numpy.fft has a function to
    # create an x-axis made out of omegas instead of mode numbers. See below.
    w_array = rfftfreq(len(th_for_Fourier), t[1]-t[0])*2*pi

    An = An[w_array <= 10*wd]  # no need to plot higher frequencies
    w_array = w_array[w_array <= 10*wd]  # no need to plot higher frequencies

    plt.semilogy(w_array/wd, abs(An), label='$A_n$')
    plt.xlabel(r'$\omega/\omega_d$', fontsize=ftsz)
    plt.ylabel(r'$|A_n|$', fontsize=ftsz)
    plt.title(r"Fourier coeffs for " +
              r"$\omega_0/\omega_d = {0:.1f}$, ".format(w0/wd) +
              r"$Q = {0:.0f}$, $\beta = {1:.3f}$".format(w0/g, b),
              fontsize=ftsz)
    plt.xlim(0., 10.)
    plt.axvline(1., ls=":", lw=1., c='k')
    plt.grid()
    plt.tight_layout()
    plt.show()

    return


def Poincare_section(n_cyc, n_per_Td, it_min, th, dth, t, wd, w0, g, b, ftsz):
    """ draw Poincare section; ncyc is how many forcing cycles we computed,
    n_per_Td is how many iterations per driving cycle, it_min is where we start
    (to skip transient), th is theta, dth is dtheta/dt  """

    # We select which points we need: start at it_min, skip every n_per_cyc
    PC_t = t[it_min::n_per_Td]*2*pi/wd  # the time sub-array
    Pc_th = th[it_min::n_per_Td]  # theta sub-array
    PC_dth = dth[it_min::n_per_Td]  # dtheta/dt sub-array

    Pc_th_wrapped = wrap_theta(Pc_th)
    plot_phase(Pc_th_wrapped, PC_dth, wd, w0, b, g, PC_t, ftsz, nconts=32)
    return


def wrap_theta(th):
    """ map all theta values onto the branch -pi<theta<pi """
    th_wrapped = th[:] % (2*pi)  # First, move everything between [0, 2*pi]
    # second, move everything between [pi, 2*pi] in [-pi, 0]
    for i in range(len(th)):
        if th_wrapped[i] > pi:
            th_wrapped[i] = th_wrapped[i] - 2*pi
    return th_wrapped


def plot_bifurcation(bmin, bmax, nbs, wd, w0, g, th0, dth0, n_per_Td, ftsz):
    """ Bifurcation plot for the DDP; produces two: one for velocity,
    one for angle """

    # Redefine time quantities because bifurcation diagrams need more points
    Td = 2*pi/wd  # driving period
    dt = Td/n_per_Td  # time step
    num_cycles = 300  # number of driving periods
    t_end = num_cycles*Td  # final time

    t = np.arange(0., t_end, dt)  # initialize the time array

    it_min = 100*n_per_Td  # to eliminate transients

    # we will calculate time series for each of these
    beta_array = np.linspace(bmin, bmax, nbs)
    # this empty array will be the beta values of each dot
    x_values = np.empty(0)

    f1, ax1 = plt.subplots(1, 1, dpi=300)  # figure for theta
    f2, ax2 = plt.subplots(1, 1, dpi=300)  # figure for dtheta/dt

    for beta in beta_array:  # scan through values of beta
        thvals, dthvals = generate_time_series(th0, dth0, wd, w0, beta, g, t)

        # Select which points we need, like in the Poincare section
        # start at it_min, skip every n_per_cyc
        Pc_thvals = wrap_theta(thvals[it_min::n_per_Td])/pi
        Pc_dthvals = dthvals[it_min::n_per_Td]/pi/wd

        # We now have computed the values for this beta, we add to the plot
        # Below is an array of same length as Pc_vals full of beta
        x_values = np.full(len(Pc_thvals), beta)
        ax1.plot(x_values, Pc_thvals, 'b.', markersize=.5)
        ax2.plot(x_values, Pc_dthvals, 'g.', markersize=.5)

    # Finish the plot for theta
    ax1.set_xlabel(r'$\beta$ [rad]', fontsize=ftsz)
    ax1.set_ylabel(r'$\theta/\pi$', fontsize=ftsz)
    ax1.set_title(r"Bifurcation plot for " +
                  r"$\omega_0/\omega_d = {0:.1f}$, ".format(w0/wd) +
                  r"$Q = {0:.0f}$".format(w0/g), fontsize=ftsz)

    if abs(bmax-1.083) <= 0.0001:
        ax1.set_ylim(-0.04, -0.0175)
    ax1.grid()
    f1.tight_layout()
    f1.savefig('ddp_bifurcation_bmin{0:.0f}_bmax{1:.0f}_theta.png'.format(
        1000*bmin, 1000*bmax))

    # Finish the plot for dtheta
    ax2.set_xlabel(r'$\beta$ [rad]', fontsize=ftsz)
    ax2.set_ylabel(r'$\dot\theta/(\pi\omega_d)$', fontsize=ftsz)
    ax2.set_title(r"Bifurcation plot for " +
                  r"$\omega_0/\omega_d = {0:.1f}$, ".format(w0/wd) +
                  r"$Q = {0:.0f}$".format(w0/g), fontsize=ftsz)
    if abs(bmax-1.083) <= 0.0001:
        ax2.set_ylim(0.93, 0.9675)
    ax2.grid()
    f2.tight_layout()
    f2.savefig('ddp_bifurcation_bmin{0:.0f}_bmax{1:.0f}_dtheta.png'.format(
        1000*bmin, 1000*bmax))

    # plt.show()  # showing would pause the calculation
    return
