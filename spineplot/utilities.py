import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

def mark_pot(ax, exposure, horizontal=False, vadj=0):
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
    matplotlib.text.Text or None
        The POT text (None if the exposure is not usable).
    """
    if exposure is None or not np.isfinite(exposure) or exposure <= 0:
        return None
    mag = int(np.floor(np.log10(exposure)))
    usepot = exposure/10**mag
    s = f'{usepot:.2f}'+f'$\\times 10^{{{mag}}}$ POT'
    xrange = ax.get_xlim()
    yrange = ax.get_ylim()
    if horizontal:
        usey = 1.01 + vadj
        return ax.text(x=1.0 - 0.025, y=usey, s=s, fontsize=13, color='black',
                       transform=ax.transAxes,
                       horizontalalignment='right', verticalalignment='bottom',
                       clip_on=False)
    else:
        usey = yrange[1] + 0.02*(yrange[1] - yrange[0])
        usex = xrange[0] - 0.12*(xrange[1] - xrange[0])
        return ax.text(x=usex, y=usey, s=s, fontsize=13, color='black', verticalalignment='top', rotation=90)

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
        The vertical adjustment in axes coordinates. The
        default is 0.
    hadj : float, optional
        The horizontal adjustment in axes coordinates. The
        default is 0.

    Returns
    -------
    None.
    """
    usex = 0.025 + hadj
    usey = 1.01 + vadj
    ax.text(
        x=usex,
        y=usey,
        s=label,
        fontsize=14,
        color='#d67a11',
        transform=ax.transAxes,
        horizontalalignment='left',
        verticalalignment='bottom',
        clip_on=False,
    )

def _axes_bbox(ax, artist):
    """
    Returns the bounding box of an artist in axes coordinates.
    """
    renderer = ax.figure.canvas.get_renderer()
    return artist.get_window_extent(renderer=renderer).transformed(ax.transAxes.inverted())

_LEGEND_OPTS = dict(loc='upper right', fontsize=10, handlelength=1.5,
                    handletextpad=0.4, columnspacing=1.0, labelspacing=0.3,
                    borderaxespad=0.3)

def draw_legend(ax, handles=None, labels=None, logy=False, headroom=True,
                max_ncol=4):
    """
    Draw a compact legend in the upper right corner of the axes and,
    if requested, raise the upper y limit so that all data sits below
    it. The number of columns is the one giving the shortest legend
    that still fits inside the axes width.

    The inputs are stored on the axis so the layout can be redone with
    `relayout_figure` after the figure is resized (the legend size is
    fixed in points, so its share of the axes depends on figure size).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to add the legend to.
    handles, labels : list, optional
        The legend handles and labels. Taken from the axis if omitted.
    logy : bool, optional
        Whether the y-axis is logarithmic. The default is False.
    headroom : bool, optional
        Raise the upper y limit to clear the legend. The default is True.
    max_ncol : int, optional
        The maximum number of columns. The default is 4.

    Returns
    -------
    matplotlib.legend.Legend or None
        The legend (None if there is nothing to show).
    """
    if handles is None or labels is None:
        handles, labels = ax.get_legend_handles_labels()
    if len(labels) == 0:
        return None
    ax._spine_legend_layout = dict(handles=list(handles), labels=list(labels),
                                   logy=logy, headroom=headroom,
                                   max_ncol=max_ncol, ylim=ax.get_ylim())
    return _layout_legend(ax)

def _layout_legend(ax):
    """
    Place the legend and set the headroom from the inputs stored by
    `draw_legend`.
    """
    cfg = ax._spine_legend_layout
    handles, labels = cfg['handles'], cfg['labels']
    ax.set_ylim(*cfg['ylim'])
    best, best_height = 1, None
    for ncol in range(1, min(cfg['max_ncol'], len(labels)) + 1):
        bb = _axes_bbox(ax, ax.legend(handles, labels, ncol=ncol, **_LEGEND_OPTS))
        if bb.x0 < 0:
            break
        if best_height is None or bb.height < best_height - 1e-6:
            best, best_height = ncol, bb.height
    legend = ax.legend(handles, labels, ncol=best, **_LEGEND_OPTS)
    if cfg['headroom']:
        make_headroom(ax, [legend], logy=cfg['logy'])
    return legend

def relayout_figure(fig) -> None:
    """
    Redo the legend placement and headroom for every axis of a figure
    that used `draw_legend`. Call this after resizing the figure.
    """
    for ax in fig.axes:
        if getattr(ax, '_spine_legend_layout', None) is not None:
            _layout_legend(ax)

def make_headroom(ax, artists, logy=False, margin=0.03) -> None:
    """
    Raise the upper y limit so that all data on the axis sits below the
    given artists (e.g. the legend) drawn in the top of the axes. The
    limit is only ever raised.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to adjust.
    artists : list
        Artists placed in axes coordinates near the top of the axes.
        None entries are ignored.
    logy : bool, optional
        Whether the y-axis is logarithmic. The default is False.
    margin : float, optional
        Gap in axes fraction between the data and the artists. The
        default is 0.03.

    Returns
    -------
    None.
    """
    artists = [a for a in artists if a is not None]
    if not artists:
        return
    frac = min(_axes_bbox(ax, a).y0 for a in artists) - margin
    if frac <= 0.2:
        frac = 0.2
    top = ax.dataLim.y1
    low, high = ax.get_ylim()
    if not np.isfinite(top):
        return
    if logy:
        if top <= 0 or low <= 0:
            return
        new_high = 10**(np.log10(low) + (np.log10(top) - np.log10(low)) / frac)
    else:
        if top <= low:
            return
        new_high = low + (top - low) / frac
    if new_high > high:
        ax.set_ylim(low, new_high)

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