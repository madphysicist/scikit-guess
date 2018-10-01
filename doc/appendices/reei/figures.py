"""
This file implements the code used to generate figures for the paper.

When imported, the figures and tables are saved to files. When run as a
script, they are displayed as Matplotlib interactive figures and
terminal output, respectively.
"""

from itertools import chain, repeat, starmap
from os import makedirs
from os.path import join

import numpy as np
from scipy.linalg import lstsq
from scipy.special import erf, erfinv
from matplotlib import pyplot as plt
from matplotlib import ticker, rc, cycler

from skg import gauss_cdf_fit
from skg import gauss_pdf_fit


OUTPUT_FOLDER = 'generated/reei'


#############
# Utilities #
#############


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


def save_fig(name, fig):
    fig.savefig(join(OUTPUT_FOLDER, name + '.png'), figsize=(6, 4.5),
                dpi=300, bbox_inches='tight')


def save_table(name, string):
    with open(join(OUTPUT_FOLDER, name + '.rst'), 'w') as file:
        file.write(string)


def gen_table(cols, specs=None, heading=None):
    """
    Generate a sphinx table of the selected data.

    Tables may contain a heading, which will be filled in from the first
    element of each column in `cols`. Headings will always be treated as
    strings.

    Parameters
    ----------
    cols : list[Sequence]
        A matrix of columns containing data. Each element should be a
        sequence of fixed (but not necessarily equal) length.
    specs : list or None
        A sequence of the same length as `cols`, or None to indicate
        default formatting. The elements are format specs that will be
        applied to the elements of the column they correspond to.
    heading : list or None
        A list of heading elements to prepend to the list.

    Return
    ------
    table : str
        A string containing an rsT/sphinx compatible table of the data
        in `cols`. The output always ends with a newline.
    """
    def filler(fill='-'):
        return tr((fill * w for w in widths), '+', fill)

    def tr(items, sep='|', fill=' '):
        mid = '{1}{0}{1}'.format(sep, fill).join(items)
        return '{0}{1}{2}{1}{0}\n'.format(sep, fill, mid)

    def pad(row):
        return starmap('{:<{}}'.format, zip(row, widths))

    def format(item, spec):
        if item in (None, ''):
            return ''
        try:
            if isinstance(item, tuple):
                return spec.format(*item)
            return spec.format(item)
        except (ValueError, TypeError):
            return str(item)

    # First format the content
    if specs is None:
        specs = repeat('{}', len(cols))
    elif isinstance(specs, str):
        specs = repeat(specs, len(cols))
    else:
        specs = ('{}' if spec is None else spec for spec in specs)
    height = len(max(cols, key=len))

    content = [[format(item, spec)
                    for item in chain(col, repeat('', height - len(col)))]
               for col, spec in zip(cols, specs)]

    # Make sure the headings are formatted too when computing widths
    if heading:
        heading = list(map(str, heading))
        widths = [max(len(h), max(map(len, col)))
                      for h, col in zip(heading, content)]
    else:
        widths = [max(map(len, col)) for col in content]

    # Pad all the cells out and transpose the list
    rows = [pad(row) for row in zip(*content)]

    prefix = filler()
    if heading:
        prefix += tr(pad(heading)) + filler('=')
    table = filler().join(map(tr, rows))
    suffix = filler()
    return prefix + table + suffix


####################
# Matplotlib Setup #
####################

rc('text', usetex=True)
plt.ioff()


#############
# Gauss PDF #
#############

def gauss_pdf():
    """
    Generates a figure and table for the Gauss PDF data in the paper.
    """
    def Smodel(x, mu, sigma):
        return 0.5 * erf((x - mu) / (np.sqrt(2) * sigma))

    def Tmodel(x, mu, sigma):
        return mu * Smodel(x, mu, sigma) - \
            sigma**2 * gauss_pdf_fit.model(x, mu, sigma)

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
    fig, ax = plt.subplots()
    format_plot(fig, ax)

    ax.text(-0.2, 1.1, '$f(x)$', fontdict={'size': 14})
    ax.text(1.1, -0.1, '$x$', fontdict={'size': 14})

    # f(x) = 1/(sigma sqrt(2 * pi) * exp(-1/2 * ((x - mu) / sigma)**2))
    ax.plot(x, y, '+', markersize=10)
    ax.plot(domain, gauss_pdf_fit.model(domain, *exact), '--')
    ax.plot(domain, gauss_pdf_fit.model(domain, *fit), '-')
    annotate(ax, '$f_k$', xy=(x[7], y[7]), xytext=(0.42, 0.45))

    # S(x) = mu * erf((x - mu) / (sqrt(2) * sigma)) / (2 * sigma**2)
    ax.plot(x, S, 's-', markersize=7, markerfacecolor='none')
    ax.plot(domain, Smodel(domain, *exact) - Smodel(x[0], *exact), '--')
    #ax.plot(domain, Sk(domain, *fit) - Sk(x[0], *fit), '-')
    annotate(ax, '$S_k$', xy=(x[7], S[7]), xytext=(0.42, 0.72))

    # T(x) = mu * S(x) - sigma**2 * f(x)
    ax.plot(x, T, 'D-', markersize=7, markerfacecolor='none')
    ax.plot(domain, Tmodel(domain, *exact) - Tmodel(x[0], *exact), '--')
    #ax.plot(domain, Tk(domain, *fit) - Tk(x[0], *fit), '-')
    annotate(ax, '$T_k$', xy=(x[7], T[7]), xytext=(0.42, -0.18))

    extra = ['', ('sigma_e', exact[1]), ('mu_e', exact[0]),
             '', ('sigma_1', fit[1]), ('mu_1', fit[0])]
    return fig, gen_table(
        cols=[np.arange(x.size) + 1, x, y, S, T, extra], specs=[
            '{:d}', ' {: 0.3g}', '{: 0.3g}', '{: 0.6g}', '{: 0.6g}',
            ':math:`\\{}` = {:< 0.6g}'
        ], heading=[
            ':math:`k`', ':math:`x_k`', ':math:`f_k`',
            ':math:`arfErf(2 F_k - 1)`', ''
        ]
    )


#############
# Gauss CDF #
#############

def gauss_cdf():
    """
    Generates a figure and table for the Gauss CDF data in the paper.
    """
    exact = 0.3, 0.4

    x = np.array([-0.914, -0.556, -0.49, -0.195, 0.019,
                   0.045,  0.587,  0.764,  0.81, 0.884])
    y = np.array([0.001, 0.017, 0.021, 0.097, 0.258,
                  0.258, 0.704, 0.911, 0.911, 0.979])

    a = np.stack((x, np.ones_like(x)), axis=1)
    b = erfinv(2 * y - 1)

    (A, B), *_ = lstsq(a, b, overwrite_a=True, overwrite_b=True)
    fit = np.array([-B / A, 1 / (np.sqrt(2.0) * A)])

    domain = np.linspace(-1.0, 1.0, 1000)
    fig, ax = plt.subplots()
    format_plot(fig, ax)

    ax.text(-0.2, 1.1, '$F(x)$', fontdict={'size': 14})
    ax.text(1.1, -0.1, '$x$', fontdict={'size': 14})

    # F(x) = 1/2 (1 + erf((x - mu) / (sqrt(2) * sigma))
    ax.plot(x, y, '+', markersize=10)
    ax.plot(domain, gauss_cdf_fit.model(domain, *exact), '--')
    ax.plot(domain, gauss_cdf_fit.model(domain, *fit), '-')

    # asymptote
    ax.plot([0, 1], [1, 1], 'k-', linewidth=0.5)

    # y(x) = erfinv(2F(x) - 1)
    b = erfinv(2.0 * y - 1)

    extra = ['', ('sigma_e', exact[1]), ('mu_e', exact[0]),
             '', ('sigma_1', fit[1]), ('mu_1', fit[0])]
    return fig, gen_table(
        cols=[np.arange(x.size) + 1, x, y, b, extra], specs=[
            '{:d}', '{: 0.3g}', '{: 0.3g}', '{: 0.6g}',
            ':math:`\\{}` = {: 0.6g}'
        ], heading=[
            ':math:`k`', ':math:`x_k`', ':math:`F_k`',
            ':math:`argErf(2 F_k - 1)`', ''
        ]
    )

############
# Erf Test #
############

def erf_test():
    x = [0.001, 0.1, 1.0, 2.0, 2.699, 2.701, 4.0, 5.0]
    y = [0.001128378791, 0.112462916, 0.8427007929, 0.995322265,
         0.9998648953, 0.998664351, 0.9999999846, 0.9999999999984]
    return None, gen_table(
        cols=[x, y], specs=['{:0.4g}', '{:0.13g}'],
        heading=[':math:`x = argErf(y)`', ':math:`Erf(x) = y`']
    )


if __name__ != '__main__':
    makedirs(OUTPUT_FOLDER, exist_ok=True)

for func in [gauss_pdf, gauss_cdf, erf_test]:
    name = func.__name__.replace('_', '-')
    title = name.replace('-', ' ').upper()
    figure, table = func()
    if __name__ == '__main__':
        if table:
            print(title, table, sep='\n\n')
        if figure:
            figure.suptitle(title)
    else:
        if figure:
            save_fig('{}-plot'.format(name), figure)
        if table:
            save_table('{}-data'.format(name), table)

if __name__ == '__main__':
    plt.show()
