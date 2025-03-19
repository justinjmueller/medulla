import numpy as np
import matplotlib.pyplot as plt

def mark_pot(ax, exposure, horizontal=False, vadj=0) -> None:
    """
    Add the POT information to the plot. The POT information will
    either be added to the top right corner of the plot along the
    horizontal axis or the top left corner of the plot along the
    vertical axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to add the POT information to.
    exposure : float
        The POT value to use in the marking.
    horizontal : bool, optional
        A flag to indicate if the POT information should be added
        along the horizontal axis. The default is True.
    vadj : float, optional
        The vertical adjustment to use when adding the POT information.
        The default is 0.

    Returns
    -------
    None.
    """
    mag = int(np.floor(np.log10(exposure)))
    usepot = exposure/10**mag
    s = f'{usepot:.2f}'+f'$\\times 10^{{{mag}}}$ POT'
    xrange = ax.get_xlim()
    yrange = ax.get_ylim()
    if horizontal:
        usey = yrange[1] + 0.01*(yrange[1] - yrange[0]) + vadj*(yrange[1] - yrange[0])
        usex = xrange[1] - 0.02*(xrange[1] - xrange[0])
        ax.text(x=usex, y=usey, s=s, fontsize=13, color='black', horizontalalignment='right')
    else:
        usey = yrange[1] + 0.02*(yrange[1] - yrange[0])
        usex = xrange[0] - 0.12*(xrange[1] - xrange[0])
        ax.text(x=usex, y=usey, s=s, fontsize=13, color='black', verticalalignment='top', rotation=90)

def mark_preliminary(ax, label, vadj=0, hadj=0) -> None:
    """
    Add a preliminary label to the plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to add the preliminary label to.
    label : str
        The label to add to the plot to indicate that the plot is
        preliminary.
    vadj : float, optional
        The vertical adjustment to use when adding the label. The
        default is 0.
    hadj : float, optional
        The horizontal adjustment to use when adding the label. The
        default is 0.

    Returns
    -------
    None.
    """
    yrange = ax.get_ylim()
    usey = yrange[1] + 0.01*(yrange[1] - yrange[0]) + vadj*(yrange[1] - yrange[0])
    xrange = ax.get_xlim()
    usex = xrange[0] + 0.025*(xrange[1] - xrange[0]) + hadj*(xrange[1] - xrange[0])
    ax.text(x=usex, y=usey, s=label, fontsize=14, color='#d67a11')

def draw_error_boxes(ax, x, y, xerr, yerr, **kwargs):
    """
    Adds error boxes to the input axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to which the error boxes are to be added.
    x : numpy.array
        The x-coordinates of the error boxes.
    y : numpy.array
        The y-coordinates of the error boxes.
    xerr : numpy.array
        The x-error values of the error boxes.
    yerr : numpy.array
        The y-error values of the error boxes.
    kwargs: dict
        Keyword arguments to be passed to the errorbar function.

    Returns
    -------
    None.
    """
    boxes = [Rectangle((x[i] - xerr[i], y[i] - yerr[i]), 2 * np.abs(xerr[i]), 2 * yerr[i]) for i in range(len(x))]
    pc = PatchCollection(boxes, **kwargs)
    ax.add_collection(pc)