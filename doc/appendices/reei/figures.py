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
from matplotlib.axis import Ticker
from matplotlib.ticker import (
    FuncFormatter, FixedLocator, MultipleLocator, NullFormatter, NullLocator,
    ScalarFormatter
)
from matplotlib.transforms import Affine2D
from matplotlib import rc, cycler

from skg import exp_fit, gauss_cdf_fit, gauss_pdf_fit, weibull_cdf_fit


OUTPUT_FOLDER = 'generated/reei'


#############
# Utilities #
#############


def format_plot(aspect='equal', x_zero=True, y_zero=True,
                majx=0.5, majy=0.5, minx=0.1, miny=0.1, figsize=(6.0, 4.5)):
    fig, ax = plt.subplots(figsize=figsize)
    if aspect is not None:
        ax.set_aspect(aspect)
    if y_zero:
        ax.spines['left'].set_position('zero')
    if x_zero:
        ax.spines['bottom'].set_position('zero')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.xaxis.set_minor_locator(MultipleLocator(minx))
    ax.xaxis.set_major_locator(MultipleLocator(majx))
    ax.yaxis.set_minor_locator(MultipleLocator(miny))
    ax.yaxis.set_major_locator(MultipleLocator(majy))

    ax.set_prop_cycle(cycler('color', 'k'))
    return fig, ax


def fix_plot_zeros(ax, offset=25):
    """
    Remove the zero on the y-axis, and shift the x-axiz zero label to
    the left.
    """
    ax.figure.canvas.draw()

    class XFormatter(ScalarFormatter):
        def __init__(self, z='$0$'):
            super().__init__()
            self.z = z
        def __call__(self, x, pos=None):
            if x == 0:
                return self.z
            return super().__call__(x, pos)

    base_trans = ax.get_xticklabels()[0].get_transform()
    def movelabel(evt=None):
        for tick in ax.xaxis.get_major_ticks():
            if tick.get_loc() == 0:
                if isinstance(offset, str):
                    ha = offset
                    trans = base_trans
                else:
                    ha = 'right'
                    trans = base_trans + Affine2D().translate(-offset, 0.0)
            else:
                trans = base_trans
                ha = 'center'
            plt.setp(tick.label, transform=trans, ha=ha)

    ax.xaxis.set_major_formatter(XFormatter(z='0'))
    ax.yaxis.set_major_formatter(XFormatter(z=''))
    movelabel()
    ax.callbacks.connect('xlim_changed', movelabel)


def xlabel(ax, label, adjust=50, xloc=1.0):
    ax.set_xlabel(label, x=xloc, ha='center', va='top', labelpad=0, fontsize=12)
    ax.xaxis.get_label().set_transform(
        ax.xaxis.get_label().get_transform() + Affine2D().translate(0, adjust)
    )


def ylabel(ax, label, adjust=50, yloc=1.0, va='center'):
    ax.set_ylabel(label, rotation=0, y=yloc, ha='right', va=va,
                  labelpad=0, fontsize=12)
    ax.yaxis.get_label().set_transform(
        ax.yaxis.get_label().get_transform() +
        Affine2D().translate(adjust, 0.0)
    )


def annotate(ax, text, xy, xytext):    
    ax.annotate(text, xy=xy, xytext=xytext, fontsize=14,
                arrowprops=dict(facecolor='k', arrowstyle='->'))


def save_fig(name, fig):
    fig.savefig(join(OUTPUT_FOLDER, name + '.png'),
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
    fig, ax = format_plot()
    xlabel(ax, '$x$')
    ylabel(ax, '$f(x)$', adjust=55)

    # f(x) = 1/(sigma sqrt(2 * pi) * exp(-1/2 * ((x - mu) / sigma)**2))
    ax.plot(x, y, '+', markersize=8)
    ax.plot(domain, gauss_pdf_fit.model(domain, *exact), '--', lw=0.5)
    ax.plot(domain, gauss_pdf_fit.model(domain, *fit), '-', lw=0.5)
    annotate(ax, '$f_k$', xy=(x[7], y[7]), xytext=(0.42, 0.45))

    # S(x) = mu * erf((x - mu) / (sqrt(2) * sigma)) / (2 * sigma**2)
    ax.plot(x, S, 's-', markersize=7, markerfacecolor='none', lw=0.5)
    ax.plot(domain, Smodel(domain, *exact) - Smodel(x[0], *exact), '--', lw=0.5)
    #ax.plot(domain, Sk(domain, *fit) - Sk(x[0], *fit), '-')
    annotate(ax, '$S_k$', xy=(x[7], S[7]), xytext=(0.42, 0.72))

    # T(x) = mu * S(x) - sigma**2 * f(x)
    ax.plot(x, T, 'D-', markersize=7, markerfacecolor='none', lw=0.5)
    ax.plot(domain, Tmodel(domain, *exact) - Tmodel(x[0], *exact), '--', lw=0.5)
    #ax.plot(domain, Tk(domain, *fit) - Tk(x[0], *fit), '-')
    annotate(ax, '$T_k$', xy=(x[7], T[7]), xytext=(0.42, -0.18))

    fix_plot_zeros(ax)

    extra = ['', ('sigma_e', exact[1]), ('mu_e', exact[0]),
             '', ('sigma_1', fit[1]),   ('mu_1', fit[0])]
    return fig, gen_table(
        cols=[np.arange(x.size) + 1, x, y, S, T, extra], specs=[
            '{:d}', ' {: 0.3g}', '{: 0.3g}', '{: 0.6g}', '{: 0.6g}',
            ':math:`\\{}` = {:< 0.6g}'
        ], heading=[
            ':math:`k`', ':math:`x_k`', ':math:`f_k`',
            ':math:`S_k`', ':math:`T_k`', ''
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
    fig, ax = format_plot()
    ax.set_ylim(0, 1.1)
    xlabel(ax, '$x$')
    ylabel(ax, '$F(x)$', adjust=60)

    # F(x) = 1/2 (1 + erf((x - mu) / (sqrt(2) * sigma))
    ax.plot(x, y, '+', markersize=6)
    ax.plot(domain, gauss_cdf_fit.model(domain, *exact), '--', lw=0.5)
    ax.plot(domain, gauss_cdf_fit.model(domain, *fit), '-', lw=0.5)

    # asymptote
    ax.plot([0, 1], [1, 1], 'k-', linewidth=0.5)

    # y(x) = erfinv(2F(x) - 1)
    b = erfinv(2.0 * y - 1)

    fix_plot_zeros(ax, offset='center')

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
    """
    Generates a table of test data for the Erf and argErf listings in
    the paper.

    Does not generate a figure.
    """
    x = [0.001, 0.1, 1.0, 2.0, 2.699, 2.701, 4.0, 5.0]
    y = [0.001128378791, 0.112462916, 0.8427007929, 0.995322265,
         0.9998648953, 0.998664351, 0.9999999846, 0.9999999999984]
    return None, gen_table(
        cols=[x, y], specs=['{:0.4g}', '{:0.13g}'],
        heading=[':math:`x = argErf(y)`', ':math:`Erf(x) = y`']
    )


#######
# Exp #
#######

def exp():
    """
    Generates a figure and table for the exponential data in the paper.
    """
    exact = 0.3, 0.6, 1.7

    x = np.array([
        -0.99, -0.945, -0.874, -0.859, -0.64,
        -0.573, -0.433, -0.042, -0.007, 0.054,
        0.088, 0.222, 0.401, 0.465, 0.633,
        0.637, 0.735, 0.762, 0.791, 0.981,
    ])
    y = np.array([
        0.418, 0.412, 0.452, 0.48, 0.453,
        0.501, 0.619, 0.9, 0.911, 0.966,
        0.966, 1.123, 1.414, 1.683, 2.101,
        1.94, 2.473, 2.276, 2.352, 3.544,
    ])

    s = np.concatenate(([0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))))

    M = np.stack((x - x[0], s), axis=1)
    Y = y - y[0]

    (A, B), *_ = lstsq(M, Y, overwrite_a=True, overwrite_b=True)

    a, c = -A / B, B

    m = np.stack((np.ones_like(x), np.exp(c * x)), axis=1)

    (a, b), *_ = lstsq(m, y, overwrite_a=True, overwrite_b=False)

    fit = np.array([a, b, c])

    domain = np.linspace(-1.0, 1.0, 1000)
    fig, ax = format_plot(aspect=0.5, majy=1, miny=0.5)
    ax.set_ylim(0, 4)
    xlabel(ax, '$x$')
    ylabel(ax, '$y$', adjust=15)

    # y(x) = a + b * exp(c * x)
    ax.plot(x, y, '+', markersize=6)
    ax.plot(domain, exp_fit.model(domain, *exact), '--', lw=0.5)
    ax.plot(domain, exp_fit.model(domain, *fit), '-', lw=0.5)

    fix_plot_zeros(ax, offset='center')
    next(tick for tick in ax.yaxis.get_major_ticks()
         if tick.get_loc() == 4.0).set_visible(False)

    extra = [
        '', ('a_e', exact[0], 1), ('b_e', exact[1], 1), ('c_e', exact[2], 2),
        '', ('a_2', fit[0], 6), ('b_2', fit[1], 6), ('c_2', fit[2], 7)
    ]
    return fig, gen_table(
        cols=[
            np.arange(x.size) + 1, x, list(zip(y, 3 + (y >= 1))), s, extra
        ], specs=[
            '{:d}', '{: 0.3g}', '{: 0.{}g}', '{: 0.6g}',
            ':math:`{}` = {: 0.{}g}'
        ], heading=[
            ':math:`k`', ':math:`x_k`', ':math:`y_k`', ':math:`S_k`', ''
        ]
    )


###############
# Weibull CDF #
###############

def weibull_cdf():
    """
    Generates a figure and table for the Weibull CDF data in the paper.
    """
    exact = 2.4, 1.6, 0.8

    t = np.array([
        1.202, 1.397, 1.537, 1.57, 1.768, 1.856, 1.87, 1.889, 1.918, 2.098,
        2.212, 2.349, 2.453, 2.557, 2.596, 2.602, 2.678, 2.706, 3.089, 3.441
    ])
    F = np.array([
        0.033, 0.082, 0.131, 0.181, 0.23, 0.279, 0.328, 0.377, 0.426, 0.475,
        0.524, 0.573, 0.622, 0.671, 0.72, 0.769, 0.818, 0.867, 0.916, 0.965
    ])

    x = np.log(-np.log(1.0 - F))
    y = np.log(t)
    s = np.concatenate(([0], np.cumsum(0.5 * (t[1:] + t[:-1]) * np.diff(x))))

    a, b, c = exp_fit(x, t, sorted=True)

    fit = np.array([1.0 / c, b, a])

    def xtrans(t):
        return np.log(t)
    def ytrans(F):
        return np.log(-np.log(np.subtract(1.0, F)))

    domain_e = np.linspace(1.125, 3.5, 1000)
    domain_f = np.linspace(1.15, 3.5, 1000)

    fig, ax = format_plot(figsize=(4.5, 6.0), aspect=None, x_zero=False,
                          majy=1.0, miny=1.0, majx=0.5, minx=0.1)
    ax.spines['bottom'].set_bounds(0, 1.8)
    ax.spines['bottom'].set_position(('data', -4))
    ax.set_xlim(0, 1.5)
    ax.spines['left'].set_bounds(-4, 1.55)
    ax.set_ylim(*ytrans([0.01, 0.995]))
    xlabel(ax, '$ln(t)$', xloc=1.15, adjust=40)
    ylabel(ax, '$ln(-ln(1-F))$', va='bottom', adjust=325, yloc=0.99)

    # Add parasite axes
    fig.subplots_adjust(left=0.25, bottom=0.25)
    par_x = ax.twinx()
    par_y = ax.twiny()

    # Configure the x-axis

    # Turn off spines and y-axis
    plt.setp(list(par_x.spines.values()) + [par_x.yaxis], visible=False)
    # Turn on and offset lower spine
    plt.setp(par_x.spines['bottom'], position=('axes', -0.01), visible=True)
    plt.setp(par_x.xaxis, ticks_position='bottom', visible=True)
    # Set tickers and ticks
    par_x.xaxis.major = Ticker()
    par_x.xaxis.minor = Ticker()
    plt.setp(par_x.xaxis,
        major_locator=FixedLocator(xtrans([1, 2, 3, 4])),
        major_formatter=FuncFormatter(
            lambda x, pos=None: '${:0.2g}$'.format(np.exp(x))
        ),
        minor_locator=NullLocator(),
        minor_formatter=NullFormatter()
    )
    # Set axis bounds
    par_x.spines['bottom'].set_bounds(*xtrans([1, 4.8]))
    # Add label
    xlabel(par_x, '$t$', xloc=1.015)

    # Configure the y-axis

    # Turn off spines and y-axis
    plt.setp(list(par_y.spines.values()) + [par_y.xaxis], visible=False)
    # Turn on and offset left spine
    plt.setp(par_y.spines['left'], position=('axes', -0.15), visible=True)
    plt.setp(par_y.yaxis, ticks_position='left', visible=True)
    # Set tickers and ticks
    par_y.yaxis.major = Ticker()
    par_y.yaxis.minor = Ticker()
    plt.setp(par_y.yaxis,
        major_locator=FixedLocator(
            ytrans([0.05, 0.1, 0.5, 0.9, 0.95, 0.99])
        ),
        major_formatter=FuncFormatter(
            lambda x, pos=None: '${:0.2g}$'.format(-np.expm1(-np.exp(x)))
        ),
        minor_locator=FixedLocator(ytrans(np.concatenate((
            np.arange(0.01, 0.1, 0.01),
            np.arange(0.1, 1.0, 0.1),
            np.arange(0.9, 1.0, 0.01)
        )))),
        minor_formatter=NullFormatter()
    )
    # Set axis bounds
    par_y.spines['left'].set_bounds(*ytrans([0.01, 0.995]))
    # Add label
    ylabel(par_y, '$F$', adjust=100, va='bottom')

    # Draw plots
    ax.plot(y, x, '+', markersize=6)
    ax.plot(xtrans(domain_e), ytrans(weibull_cdf_fit.model(domain_e, *exact)),
            '--', lw=0.5)
    ax.plot(xtrans(domain_f), ytrans(weibull_cdf_fit.model(domain_f, *fit)),
            '-', lw=0.5)

    extra = [
        '', ('alpha_e', exact[0], 2), ('beta_e', exact[1], 2),
            ('mu_e', exact[2], 1),
        '', ('alpha_c', fit[0], 6), ('beta_c', fit[1], 6), ('mu_c', fit[2], 6)
    ]
    return fig, gen_table(
        cols=[np.arange(x.size) + 1, t, F, x, s, extra], specs=[
            '{:d}', '{: 0.4g}', '{: 0.3g}', '{: 0.6g}', '{: 0.6g}',
            ':math:`\\{}` = {: 0.{}g}'
        ], heading=[
            ':math:`k`', ':math:`t_k`', ':math:`F_k`',
            ':math:`ln(-ln(1 - F_k))`', ':math:`S_k`', ''
        ]
    )


if __name__ != '__main__':
    makedirs(OUTPUT_FOLDER, exist_ok=True)

for func in [gauss_pdf, gauss_cdf, erf_test, exp, weibull_cdf]:
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
