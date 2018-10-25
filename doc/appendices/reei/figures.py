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


def annotate(ax, text, xy, xytext, fs=14, color='k'):
    ax.annotate(text, xy=xy, xytext=xytext, fontsize=fs, color=color,
                arrowprops=dict(color=color, arrowstyle='->'))


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
    ax.plot(domain, Smodel(domain, *exact) - Smodel(x[0], *exact),
            '--', lw=0.5)
    #ax.plot(domain, Sk(domain, *fit) - Sk(x[0], *fit), '-')
    annotate(ax, '$S_k$', xy=(x[7], S[7]), xytext=(0.42, 0.72))

    # T(x) = mu * S(x) - sigma**2 * f(x)
    ax.plot(x, T, 'D-', markersize=7, markerfacecolor='none', lw=0.5)
    ax.plot(domain, Tmodel(domain, *exact) - Tmodel(x[0], *exact),
            '--', lw=0.5)
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


#############
# Sinusoids #
#############

class Sinusoid:
    """
    Keeps track of the shared data for all the sinusoidal regressions
    described in the paper.
    """

    # a, b, c, omega
    exact = (-0.4, 1.3, -0.6, 2)

    x = np.array([
        -1.983, -1.948, -1.837, -1.827, -1.663, -0.815, -0.778, -0.754,
        -0.518, 0.322, 0.418, 0.781, 0.931, 1.51, 1.607,
    ])
    y = np.array([
        0.936, 0.81, 0.716, 0.906, 0.247, -1.513, -1.901, -1.565,
        -1.896, 0.051, 0.021, 1.069, 0.862, 0.183, 0.311,
    ])

    domain = np.linspace(-2.0, 2.0, 1000)

    def __init__(self):
        exact = self.model(self.x, *self.exact)
        self.rho0 = np.hypot(self.exact[1], self.exact[2])
        self.phi0 = np.arctan2(self.exact[2], self.exact[1])
        self.rms0 = np.std(self.y - exact)

        # a, rho, omega, phi
        self.exact2 = (self.exact[0], self.rho0, self.exact[3], self.phi0)

    @staticmethod
    def model(x, a, b, c, omega):
        """
        The model for a simple sinusoid.
        """
        t = omega * x
        return a + b * np.sin(t) + c * np.cos(t)

    @staticmethod
    def model2(x, a, rho, omega, phi):
        """
        Alternative phrasing for the same model.
        """
        return a + rho * np.sin(omega * x + phi)

    @staticmethod
    def S(x, a, b, c, omega):
        """
        First order integral of the sinusoid.
        """
        t = omega * x
        s = a * x - b / omega * np.cos(t) + c / omega * np.sin(t)
        s -= s[0]
        return s

    @staticmethod
    def SS(x, a, b, c, omega):
        """
        Second order integral of the sinusoid.
        """
        t = omega * x
        C = a * x[0] - b / omega * np.cos(t[0]) + c / omega * np.sin(t[0])
        ss = 0.5 * a * x**2 - C * x - \
                b / omega**2 * np.sin(t) - c / omega**2 * np.cos(t)
        ss -= ss[0]
        return ss

    @staticmethod
    def derivative(x, a, b, c, omega):
        """
        The derivative of the sinusoid.
        """
        t = omega * x
        return b * omega * np.cos(t) - c * omega * np.sin(t)

    @classmethod
    def fit1(cls, x, y):
        """
        Integral-only fitting function.
        """
        x, y = map(np.asfarray, (x, y))

        d = 0.5 * np.diff(x)

        S = np.cumsum((y[1:] + y[:-1]) * d)
        S = np.insert(S, 0, 0)

        SS = np.cumsum((S[1:] + S[:-1]) * d)
        SS = np.insert(SS, 0, 0)

        M = np.stack((SS, x**2, x, np.ones_like(x)), axis=1)
        (A1, B1, C1, D1), *_ = lstsq(M, y, overwrite_a=True, overwrite_b=False)

        x0 = cls.x[0]
        omega = np.sqrt(-A1) if A1 <= 0 else np.nan
        a = 2 * B1 / omega**2
        p = B1 * x0**2 + C1 * x0 + D1 - a
        q = (C1 + 2 * B1 * x0) / omega
        t = omega * x0
        b = p * np.sin(t) + q * np.cos(t)
        c = p * np.cos(t) - q * np.sin(t)
        return (a, b, c, omega)

    @classmethod
    def fit2(cls, x, y):
        """
        Second order fit using inverse tangent.
        """
        a1, b1, c1, omega1 = cls.fit1(x, y)
        if np.isnan(a1):
            return (np.nan,) * 4

        rho1 = np.hypot(c1, b1)
        phi1 = np.arctan2(c1, b1)

        Phi = cls.atan_x(y - a1, rho1)
        kk = np.round((omega1 * x + phi1) / np.pi)

        theta = (-1)**kk * Phi + np.pi * kk
        M = np.stack((x, np.ones_like(x)), axis=-1)

        (omega2, phi2), *_ = lstsq(M, theta,
                                   overwrite_a=True, overwrite_b=True)
        a2 = a1
        b2 = rho1 * np.cos(phi2)
        c2 = rho1 * np.sin(phi2)

        return (a2, b2, c2, omega2)

    @classmethod
    def fit3(cls, x, y):
        """
        Third order fit using the classical approach.
        """
        a2, b2, c2, omega2 = cls.fit2(x, y)
        t = omega2 * x
        M = np.stack((np.ones_like(x), np.sin(t), np.cos(t)), axis=-1)
        fit3, *_ = lstsq(M, y, overwrite_a=True)
        return tuple(fit3) + (omega2,)

    @staticmethod
    def atan_x(f, rho):
        r"""
        Computes :math:`arctan \left(\frac{f}{\sqrt{\rho^2 - f^2}} \right)`
        for :math:`\rho^2 > f^2`. In all other cases, returns
        :math:`\frac{\pi}{2} with the same sign as :math:`f`.
        """
        dis = rho**2 - f**2
        m = (dis >= 0)
        out = np.arctan(np.divide(f, np.sqrt(dis, where=m), where=m), where=m)
        np.copysign(np.pi / 2.0, f, where=~m, out=out)

        return out

    @classmethod
    def phi(cls, x, a, rho, omega, phi):
        r"""
        Computes

        .. math::

           \Phi(x) = arctan \left( \frac{f(x) - a}
               {\sqrt{\rho^2 - \left( f(x) - a\right)^2}} \right)

        Based on :math:`f(x) = a + \rho \; sin(\omega \; x + \varphi)`.
        """
        f = cls.model2(x, 0.0, rho, omega, phi)
        return cls.atan_x(f, rho)

    @classmethod
    def phi2(cls, x, y, a, rho):
        r"""
        Computes

        .. math::

           \Phi(x) = arctan \left( \frac{y - a}
               {\sqrt{\rho^2 - \left( y - a\right)^2}} \right)
        """
        f = y - a
        return cls.atan_x(f, rho)

    def gen_omega_cdf(self, tx, fit, ns=None, x_rand=True, y_sigma=0.0,
                      omega=1, samp=10000, an=None):
        """
        Generate a figure and table with CDFs for each number of
        points-per-cycle in `ns`.
        """
        def unpack_tx(tx):
            if isinstance(tx, tuple):
                return tx
            return tx, 0.8

        if ns is None:
            ns = [8, 10, 12, 15, 20, 50]
        ns = np.array(ns)

        omega_e = 2.0 * np.pi
        wm = []

        fig, ax = format_plot(aspect=0.15, x_zero=False, majx=0.1, minx=0.1)
        ax.spines['left'].set_position(('data', 1))
        ax.set_xlim(0.89, 1.35)
        ax.set_ylim(0, 1.099)
        xlabel(ax, r'$\frac{{\omega_{}}}{{\omega_e}}$'.format(omega))
        ylabel(ax, '$P$')

        if an is None:
            ax.text(0.5 * sum(tx[-2:]), 0.9, '$n_p =$',
                    fontsize=8, va='top', ha='left')
        else:
            annotate(ax, '$n_p$', unpack_tx(tx[0]), an, fs=8)

        for p, t in zip(ns, tx):
            x = np.random.rand(samp, p)
            x.sort(axis=1)
            y = np.sin(omega_e * x) + np.random.normal(loc=0.0, scale=y_sigma,
                                                       size=x.shape)
            # TODO: This is a good usecase for multiple simultaneous fits.
            ratios = [fit(*k)[-1] for k in zip(x, y)]
            ratios = np.array([f / omega_e for f in ratios if not np.isnan(f)])

            pdf, bins = np.histogram(ratios, bins=samp // 20, density=True)
            pdf *= np.diff(bins)
            bins = 0.5 * (bins[1:] + bins[:-1])
            cdf = np.cumsum(pdf)
            wm.append(np.interp(0.5, cdf, bins))

            ax.plot(bins, cdf, lw=0.5)
            t, u = unpack_tx(t)
            ax.text(t, u, '${}$'.format(p), fontsize=8, va='top', ha='left')

        wm = np.array(wm)

        ax.plot([1, 1.2], [0.5, 0.5], ':', lw=0.5)
        ax.plot([1, 1.4], [1, 1], lw=0.75)
        ax.text(1.21, 0.5,
                r'$P \left( \frac{{\omega_{0}}}{{\omega_e}} < '
                r'\frac{{\omega_{{{0}m}}}}{{\omega_e}} \right) = '
                '0.5$'.format(omega),
                va='center', ha='left')

        fix_plot_zeros(ax)

        return fig, gen_table(
            cols=[ns[:len(wm)], wm], specs=['{:d}', '{:0.3f}'],
            heading=[
                ':math:`n_p`',
                r':math:`\frac{{\omega_{{{}m}}}}{{\omega_e}}`'.format(omega)
            ]
        )

    def sin_exact(self):
        """
        Generates a figure and table for the exact sinusoidal data in
        the paper.
        """
        fig, ax = format_plot(majx=1, minx=1, majy=1, miny=1)
        xlabel(ax, '$x$')
        ylabel(ax, '$y$')

        ax.plot(self.x, self.y, '+', markersize=6)
        ax.plot(self.domain, self.model(self.domain, *self.exact),
                '--', lw=0.5)

        fix_plot_zeros(ax)

        extra = [
            '', ('\\omega_e', self.exact[-1], 0), ('a_e', self.exact[0], 1),
            ('b_e', self.exact[1], 1), ('c_e', self.exact[2], 1),
            ('\\rho_e', self.rho0, 6), ('\\sigma_e', self.rms0, 4),
        ]
        return fig, gen_table(
            cols=[np.arange(self.x.size) + 1, self.x, self.y, extra],
            specs=['{:d}', '{:0.3f}', '{:0.3f}', ':math:`{}` = {:.{}f}'],
            heading=[':math:`k`', ':math:`x_k`', ':math:`y_k`', '']
        )

    def sin_nomega(self):
        """
        Generates a figure and table for the sinusoid with known
        frequency in the paper.
        """
        fig, ax = format_plot(majx=1, minx=1, majy=1, miny=1)
        xlabel(ax, '$x$')
        ylabel(ax, '$y$')

        omega_e = self.exact[-1]
        t = omega_e * self.x
        M = np.stack((np.ones_like(self.x), np.sin(t), np.cos(t)), axis=1)
        p = self.y

        (a, b, c), *_ = lstsq(M, p, overwrite_a=True, overwrite_b=False)
        fit = np.array([a, b, c, omega_e])

        rho = np.hypot(b, c)
        rms = np.std(self.y - self.model(self.x, *fit))

        ax.plot(self.x, self.y, '+', markersize=6)
        ax.plot(self.domain, self.model(self.domain, *self.exact),
                '--', lw=0.5)
        ax.plot(self.domain, self.model(self.domain, *fit), '-', lw=0.5)

        fix_plot_zeros(ax)

        extra = [
            ('\\omega_e', omega_e, 0), ('a_0', fit[0], 6),
            ('b_0', fit[1], 6), ('c_0', fit[2], 6), ('\\rho_0', rho, 6),
            ('\\sigma_0', rms, 6),
        ]
        return fig, gen_table(
            cols=[extra], specs=[':math:`{}` = {:.{}f}'], heading=None
        )

    def sin_int(self):
        """
        Generates a figure and table for the integral-only sinudoidal
        data in the paper.
        """
        d = 0.5 * np.diff(self.x)

        S = np.cumsum((self.y[1:] + self.y[:-1]) * d)
        S = np.insert(S, 0, 0)

        SS = np.cumsum((S[1:] + S[:-1]) * d)
        SS = np.insert(SS, 0, 0)

        diff = np.diff(self.y) / np.diff(self.x)
        diff_x = 0.5 * (self.x[1:] + self.x[:-1])

        fig, ax = format_plot(aspect=0.4,
                              majx=1.0, minx=1.0, majy=1.0, miny=1.0)
        ax.set_ylim(-3.5, 4.1)
        xlabel(ax, '$x$')

        # Plot the derivative
        ax.plot(self.domain, self.derivative(self.domain, *self.exact),
                '--', lw=0.75)
        ax.plot(diff_x, diff, 's-', ms=1.5, lw=-0.75)
        annotate(ax, '$f\'(x)$', xy=(-1.1, self.derivative(-1.1, *self.exact)),
                 xytext=(-1.5, -2.0), fs=12)
        annotate(ax, '$f\'_k$', xy=(diff_x[1], diff[1]), xytext=(-1.4, -1.5),
                 fs=12)

        # Plot the first-order integral
        ax.plot(self.domain, self.S(self.domain, *self.exact), 'b--', lw=0.75)
        ax.plot(self.x, S, 'bs-', ms=1.5, lw=0.75)
        annotate(ax, '$S_k$', xy=(self.x[11], S[11]), xytext=(0.5, -0.8),
                 fs=12)

        # Plot the second-order integral
        ax.plot(self.domain, self.SS(self.domain, *self.exact), 'r--', lw=0.75)
        ax.plot(self.x, SS, 'rs-', ms=1.5, lw=0.75)
        annotate(ax, '$SS_k$', xy=(self.x[11], SS[11]), xytext=(0.45, -2.8),
                 fs=12)

        fix_plot_zeros(ax)

        return fig, gen_table(
            cols=[np.arange(self.x.size) + 1, S, SS],
            specs=['{:d}', '{:0.6g}', '{:0.6g}'],
            heading=[':math:`k`', ':math:`S_k`', ':math:`SS_k`']
        )

    def sin_eq_nd(self):
        r"""
        Generates a figure and table showing the effects of :math:`n_k`
        on :math:`\omega_1 / \omega_e` for the equidistant,
        non-dispersive sinusoid in the paper.
        """
        omega_e = 2.0 * np.pi
        n = np.arange(5, 21)

        def omega(p):
            x = np.linspace(0.0, 1.0, p)
            y = np.sin(omega_e * x)
            return self.fit1(x, y)[-1]

        ratios = np.array([omega(p) / omega_e for p in n])

        fig, ax = format_plot(aspect=40, x_zero=False, y_zero=False,
                              majx=5, minx=1, majy=0.1, miny=0.1)
        ax.set_xlim(4.1, 21.9)
        ax.set_ylim(1, 1.33)
        xlabel(ax, '$n_p$')
        ylabel(ax, r'$\frac{\omega_1}{\omega_e}$')

        ax.plot(n, ratios, 's', ms=1.8, markerfacecolor='none')

        return fig, gen_table(
            cols=[n, ratios], specs=['{:d}', '{:0.4g}'],
            heading=[':math:`n_p`', r':math:`\frac{\omega_1}{\omega_e}`']
        )

    def sin_rand_nd(self):
        r"""
        Generates a figure and table showing the effects of :math:`n_k`
        on :math:`\omega_1 / \omega_e` for the random, non-dispersive
        sinusoid in the paper.
        """
        return self.gen_omega_cdf(
            fit=self.fit1, tx=[1.185, 1.125, 1.09, 1.06, 1.03, 1.01]
        )

    def sin_rand_d(self):
        r"""
        Generates a figure and table showing the effects of :math:`n_k`
        on :math:`\omega_1 / \omega_e` for the random, dispersive
        sinusoid in the paper.
        """
        return self.gen_omega_cdf(
            fit=self.fit1, tx=[1.21, 1.14, 1.1, 1.07, 1.04, 1.02], y_sigma=0.1
        )

    def sin_saw(self):
        """
        Generates a figure showing the transformation of a sawtooth
        function into a line.
        """
        fit = self.fit1(self.x, self.y)

        rho1 = np.hypot(fit[1], fit[2])
        phi1 = np.arctan2(fit[2], fit[1])

        phi_exact = self.phi(self.domain, *self.exact2)
        phi_data = self.phi2(self.x, self.y, fit[0], rho1)

        kk_exact = np.round((self.exact[-1] * self.domain + self.phi0) / np.pi)
        kk_data = np.round((fit[-1] * self.x + phi1) / np.pi)

        theta_exact = (-1)**kk_exact * phi_exact + np.pi * kk_exact
        theta_data = (-1)**kk_data * phi_data + np.pi * kk_data

        fig, ax = format_plot(aspect=0.4, majx=1, minx=1, majy=1, miny=1)
        xlabel(ax, '$x_k$')

        ax.plot(self.domain, phi_exact, '--', lw=0.75)
        ax.plot(self.domain, theta_exact, '-', lw=0.4)
        ax.plot(self.x, phi_data, '+', markersize=6)
        ax.plot(self.x, theta_data, 's', lw=0.5, ms=4, markerfacecolor='none')

        annotate(ax, r'$(\Phi)$', (-1.75, 1), (-1.5, 1.6))
        annotate(ax, r'$(\theta)$', (-1.85, -3.7), (-2, -3))

        fix_plot_zeros(ax)

        return fig, gen_table(
            cols=[
                np.arange(kk_data.size) + 1, phi_data,
                kk_data.astype(np.int), theta_data
            ],
            specs=['{:d}', '{:0.6g}', '{:d}', '{:0.6g}'],
            heading=[
                ':math:`k`', r':math:`\Phi_k`', ':math:`K_k`',
                r':math:`\theta_k`'
            ]
        )

    def sin_rand_nd2(self):
        r"""
        Generates a figure and table showing the effects of :math:`n_k`
        on :math:`\omega_2 / \omega_e` for the random, non-dispersive
        sinusoid in the paper.
        """
        return self.gen_omega_cdf(
            fit=self.fit2, omega=2, tx=[
                (1.06, 0.8), (1.04, 0.82), (1.03, 0.84), (1.02, 0.86),
                (1.01, 0.88)
            ], an=(1.08, 0.75)
        )

    def sin_rand_d2(self):
        r"""
        Generates a figure and table showing the effects of :math:`n_k`
        on :math:`\omega_2 / \omega_e` for the random, dispersive
        sinusoid in the paper.
        """
        return self.gen_omega_cdf(
            fit=self.fit2, omega=2, tx=[1.21, 1.14, 1.1, 1.07, 1.04, 1.01],
            y_sigma=0.1
        )

    def sin_final(self):
        """
        Generates a figure showing the final, and all intermediate
        steps, in the optimization.
        """
        fit1 = self.fit1(self.x, self.y)
        fit2 = self.fit2(self.x, self.y)
        fit3 = self.fit3(self.x, self.y)

        fig, ax = format_plot(majx=1, minx=1, majy=1, miny=1)
        xlabel(ax, '$x$')
        ylabel(ax, '$y$')

        ax.plot(self.x, self.y, '+', markersize=6)
        ax.plot(self.domain, self.model(self.domain, *self.exact),
                '--', lw=0.5)
        ax.plot(self.domain, self.model(self.domain, *fit1), 'b-', lw=0.5)
        ax.plot(self.domain, self.model(self.domain, *fit2), 'r-', lw=0.5)
        ax.plot(self.domain, self.model(self.domain, *fit3), 'k-', lw=0.5)

        annotate(ax, '$(1)$', (0.45, self.model(0.45, *fit1)), (0.2, 1.22),
                 fs=8, color='b')
        annotate(ax, '$(2)$', (0.9, self.model(0.9, *fit2)), (0.85, 1.22),
                 fs=8, color='r')
        annotate(ax, '$(3)$', (1.15, self.model(1.15, *fit3)), (1.5, 1.22),
                 fs=8, color='k')

        fix_plot_zeros(ax)

        def rearrange(fit):
            a, b, c, omega = fit
            return [omega, a, b, c, np.hypot(c, b), np.arctan2(c, b)]

        labels = [
            r':math:`\omega`', ':math:`a`', ':math:`b`', ':math:`c`',
            r':math:`\rho`', r':math:`\varphi`',
        ]
        return fig, gen_table(
            cols=[labels, rearrange(fit1), rearrange(fit2), rearrange(fit3)],
            specs=['{}', '{:0.6g}', '{:0.6g}', '{:0.6g}'],
            heading=['', ':math:`(1)`', ':math:`(2)`', ':math:`(3)`']
        )


sinusoid = Sinusoid()
func_list = [
    gauss_pdf, gauss_cdf, erf_test, exp, weibull_cdf,
    sinusoid.sin_exact, sinusoid.sin_nomega, sinusoid.sin_int,
    sinusoid.sin_eq_nd, sinusoid.sin_rand_nd, sinusoid.sin_rand_d,
    sinusoid.sin_saw, sinusoid.sin_rand_nd2, sinusoid.sin_rand_d2,
    sinusoid.sin_final,
]


if __name__ != '__main__':
    makedirs(OUTPUT_FOLDER, exist_ok=True)

for func in func_list:
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
