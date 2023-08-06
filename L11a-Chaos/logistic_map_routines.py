"""
Routines to support the logistic map simulations
Nicolas Grisouard, University of Toronto, August 2021
"""
import numpy as np
import matplotlib.pyplot as plt


def iterate_logistic_map(mu, x0):
    xns = [x0]  # one element for now
    old_x = 0.
    x = x0
    counter = 0  # counting how many steps
    while counter < 400:
        counter += 1
        y = mu*x*(1-x)  # do step 1
        old_x = x
        x = y  # do step 2
        xns.append(y)  # add the next value to the list

    return xns


def plot_logistic_map(mu, x0, skip_n_ites=0):
    """ x0 is the first value, skip_n_ites (optional) allows one to skip the
    transient """
    x = np.linspace(0, 1, 100)  # array of x_n's to plot x_{n+1}
    y = mu*x*(1.0-x)  # array of x_{n+1}'s

    # Plot the logistic map
    plt.subplot(1, 2, 1)
    plt.plot(x, y, 'b')
    plt.plot(x, x, 'k')

    # now do the steps
    list_of_xn = iterate_logistic_map(mu, x0)
    for n in range(skip_n_ites, len(list_of_xn)-1):
        xn = list_of_xn[n]
        yn = list_of_xn[n+1]
        if n == 0:
            plt.plot([xn, xn], [0, yn], 'r--')  # plot step 1
        else:
            plt.plot([xn, xn], [xn, yn], 'r--')  # plot step 1
        plt.plot([xn, yn], [yn, yn], 'g-.')  # plot step 2

    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.xlabel('$x_n$')
    plt.ylabel('$x_{n+1}$')

    plt.subplot(1, 2, 2)
    plt.plot(list_of_xn)
    plt.xlim(skip_n_ites, len(list_of_xn)-1)
    plt.xlabel('$n$')
    plt.ylabel('$x_n$')
    plt.grid()

    plt.suptitle('Logistic Map for $\mu = {0:.1f}$'.format(mu))

    plt.tight_layout()
    plt.show()
    return


def bifurcation_logistic(mumin, mumax, skip_n_ites, ymin=0, ymax=1):
    """ Bifurcation plot """
    x0 = 0.8  # no need to vary it
    mus = np.linspace(mumin, mumax, 400)
    for mu in mus:
        list_of_xn = iterate_logistic_map(mu, x0)
        reduced_list_of_xn = list_of_xn[skip_n_ites:]
        x_values = np.full(len(reduced_list_of_xn), mu)  # array full of mu
        plt.plot(x_values, reduced_list_of_xn, 'b.', markersize=1.)

    plt.xlabel('$\mu$')
    plt.ylabel('Stationary values of $x_n$')
    plt.title('Bifurcation diagram for the logistic map')
    plt.xlim(mumin, mumax)
    plt.ylim(ymin, ymax)
    plt.grid()
    plt.tight_layout()
    plt.show()
    return


def lyapunov_logistic(mu, x01, epsilon):
    """ iterate the two maps in parallel and compute their difference at each
    step. """
    def map_it(x):
        return mu*x*(1-x)

    npts = 60
    difference = np.empty(npts)
    x1 = x01
    x2 = x01 + epsilon
    difference[0] = x2 - x1

    for n in range(1, npts):  # should be enough
        y1, y2 = map_it(x1), map_it(x2)
        difference[n] = y2-y1
        x1, x2 = y1, y2

    return difference


def plot_two_divergences(mu1, mu2, x01, epsilon):
    """ Run two divergences side-by-side """
    diff1 = lyapunov_logistic(mu1, x01, epsilon)  # first mu
    diff2 = lyapunov_logistic(mu2, x01, epsilon)  # second mu

    plt.subplot(2, 1, 1)
    plt.semilogy(abs(diff1))
    plt.xlabel('$n$')
    plt.ylabel('$|x_n^{(2)}- x_n^{(1)}|$')
    plt.title(r"$\mu = {}$".format(mu1))
    plt.grid()

    plt.subplot(2, 1, 2)
    plt.semilogy(abs(diff2))
    plt.xlabel('$n$')
    plt.ylabel('$|x_n^{(2)}- x_n^{(1)}|$')
    plt.title(r"$\mu = {}$".format(mu2))
    plt.grid()

    plt.tight_layout()
    return
