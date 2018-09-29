"""
This file implements the code used to generate figures for the paper.
"""

from itertools import zip_longest
from os import makedirs
from os.path import join

import numpy as np
from scipy.linalg import lstsq
from scipy.special import erf
from matplotlib import pyplot as plt
from matplotlib import ticker, rc, cycler

from skg import gauss_fit


SAVE = True
OUTPUT = 'generated/reei'


#############
# Utilities #
#############

rc('text', usetex=True)

def format_plot(fig, ax):
    ax.set_aspect('equal')
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))

    ax.set_prop_cycle(cycler('color', 'k'))

def annotate(ax, text, xy, xytext):    
    ax.annotate(text, xy=xy, xytext=xytext, fontsize=14,
                arrowprops=dict(facecolor='k', arrowstyle='->'))

def write_line(file, indent, line):
    file.write(' ' * indent)
    file.write(line)
    file.write('\n')

plt.ioff()

if SAVE:
    makedirs(OUTPUT, exist_ok=True)


############
# Figure 1 #
############

exact = -0.3, 0.4

x = np.array([-0.992, -0.935, -0.836, -0.404, -0.326,
              -0.042,  0.068,  0.302,  0.439,  0.58])
y = np.array([0.238, 0.262, 0.38, 1.041, 0.922,
              0.755, 0.589, 0.34, 0.193, 0.083])

S = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
S = np.insert(S, 0, 0)

xy = x * y
T = np.cumsum(0.5 * (xy[1:] + xy[:-1]) * np.diff(x))
T = np.insert(T, 0, 0)

A = np.stack((S, T), axis=1)
b = y - y[0]

(A1, B1), *_ = lstsq(A, b, overwrite_a=True, overwrite_b=True)
fit = np.array([-A1 / B1, np.sqrt(-1.0 / B1)])

domain = np.linspace(-1.0, 1.0, 1000)
fig1, ax1 = plt.subplots()
format_plot(fig1, ax1)

ax1.text(-0.2, 1.1, '$f(x)$', fontdict={'size': 14})
ax1.text(1.1, -0.1, '$x$', fontdict={'size': 14})

# f(x) = 1/(sigma sqrt(2 * pi) * exp(-1/2 * ((x - mu) / sigma)**2))
ax1.plot(x, y, 'k+', markersize=10)
ax1.plot(domain, gauss_fit.model(domain, *exact), '--')
ax1.plot(domain, gauss_fit.model(domain, *fit), '-')
annotate(ax1, '$f_k$', xy=(x[7], y[7]), xytext=(0.42, 0.45))

# S(x) = mu * erf((x - mu) / (sqrt(2) * sigma)) / (2 * sigma**2)
def Smodel(x, mu, sigma):
    return 0.5 * erf((x - mu) / (np.sqrt(2) * sigma))
ax1.plot(x, S, 's-', markersize=7, markerfacecolor='none')
ax1.plot(domain, Smodel(domain, *exact) - Smodel(x[0], *exact), '--')
#ax1.plot(domain, Sk(domain, *fit) - Sk(x[0], *fit), '-')
annotate(ax1, '$S_k$', xy=(x[7], S[7]), xytext=(0.42, 0.72))

# T(x) = mu * S(x) - sigma**2 * f(x)
def Tmodel(x, mu, sigma):
    return mu * Smodel(x, mu, sigma) - sigma**2 * gauss_fit.model(x, mu, sigma)
ax1.plot(x, T, 'D-', markersize=7, markerfacecolor='none')
ax1.plot(domain, Tmodel(domain, *exact) - Tmodel(x[0], *exact), '--')
#ax1.plot(domain, Tk(domain, *fit) - Tk(x[0], *fit), '-')
annotate(ax1, '$T_k$', xy=(x[7], T[7]), xytext=(0.42, -0.18))

if SAVE:
    fig1.savefig(join(OUTPUT, 'gauss-plot.png'), figsize=(6, 4.5), dpi=300,
                 bbox_inches='tight')
    def ex_fmt(name, value):
        return ':math:`{}` = {:< 0.6g}'.format(name, value)
    extra = ['', ex_fmt(r'\sigma_e', exact[1]), ex_fmt(r'\mu_e', exact[0]), '',
             ex_fmt(r'\sigma_1', fit[1]), ex_fmt(r'\mu_1', fit[0])]
    with open(join(OUTPUT, 'gauss-data.rst'), 'w') as table:
        write_line(table, 0,
                   '+-----------+-------------+-------------+-------------+'
                   '-------------+-------------------------------+')
        write_line(table, 0,
                   '| :math:`k` | :math:`x_k` | :math:`f_k` | :math:`S_k` |'
                   ' :math:`T_k` |                               |')
        write_line(table, 0,
                   '+===========+=============+=============+=============+'
                   '=============+===============================+')
        for k, items in enumerate(zip_longest(x, y, S, T, extra, fillvalue=''),
                                  start=1):
            write_line(table, 0,
                       '| {:<9d} | {:< 11.3g} | {:< 11.3g} | {:< 11.6g} |'
                       ' {:< 11.6g} | {:<29s} |'.format(k, *items))
            write_line(table, 0,
                       '+-----------+-------------+-------------+-----------'
                       '--+-------------+-------------------------------+')


if not SAVE:
    plt.show()
