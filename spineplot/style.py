import matplotlib.pyplot as plt
from os import environ
from os.path import dirname, isabs, isfile, join, normpath

import toml

try:
    from .utilities import mark_preliminary, mark_pot
except ImportError:
    try:
        from utilities import mark_preliminary, mark_pot
    except ImportError:
        mark_preliminary = None
        mark_pot = None

class Style:
    """
    A class designed to encapsulate the configuration of a plotting
    style for a given plot.

    Attributes
    ----------
    _name : str
        The name of the Style object.
    _style : str
        The name of the style sheet to use for the plot.
    _markers : list
        The list of markers to use for the plot.
    _default_figsize : tuple
        The default size of the figure to create.
    _title : str
        The title to place at the top of the plot.
    _mark_pot : bool
        A flag toggling the display of the total POT (exposure) at the
        top of the plot above the axis and below the title.
    _mark_preliminary : str
        A string to be used a label to indicate that the plot is
        preliminary. If None, no label is added.
    _scilimits : tuple, optional
        The scientific limits to use when formatting the axis labels.
        The default is None.
    _plot_kwargs : dict
        A dictionary containing the keyword arguments to be passed
        to the plotting function.
    """
    def __init__(self, name, style_sheet, markers, default_figsize,
                 title, mark_pot=True, mark_pot_horizontal=True,
                 mark_preliminary=True, scilimits=None, plot_kwargs=None) -> None:
        """
        Initializes the Style object with the given kwargs.

        Parameters
        ----------
        name : str
            The name of the style object.
        style_sheet : str
            The name of the style sheet to use for the plot.
        markers : list
            The list of markers to use for the plot.
        default_figsize : tuple
            The default size of the figure to create.
        title : str
            The title to place at the top of the plot.
        mark_pot : bool, optional
            A flag toggling the display of the total POT (exposure) at the
            top of the plot above the axis and below the title. The default
            is True.
        mark_preliminary : str, optional
            A string to be used a label to indicate that the plot is
            preliminary. If None, no label is added. The default is True.
        scilimits : tuple, optional
            The scientific limits to use when formatting the axis labels.
            The default is None.
        plot_kwargs : dict, optional
            A dictionary containing the keyword arguments to be passed
            to the plotting function. The default is None.

        Returns
        -------
        None
        """
        self._name = name
        self._style = f'{environ.get("MEDULLA_PLOT_DIR", "")}/{style_sheet}'
        self._markers = markers
        self._default_figsize = default_figsize
        self._title = None if title == 'none' else title
        self._mark_pot = mark_pot
        self._mark_pot_horizontal = mark_pot_horizontal
        self._mark_preliminary = None if mark_preliminary == 'none' else mark_preliminary
        self._scilimits = scilimits
        self._plot_kwargs = plot_kwargs

    def __enter__(self):
        """
        Sets the style for the plot.

        Returns
        -------
        self : Style
            The created Style object.
        """
        plt.style.use(self._style)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """
        Resets the style for the plot.

        Parameters
        ----------
        exc_type : type
            The exception type.
        exc_val : Exception
            The exception value.
        exc_tb : traceback
            The traceback.

        Returns
        -------
        None
        """
        plt.style.use('default')

    def get_style(self) -> str:
        """
        Returns the value of the style attribute.

        Returns
        -------
        str
            The value of the style attribute.
        """
        return self._style

    def get_color(self, cdx) -> str:
        """
        Returns the color for the given cycle index. This needs to
        correctly loop over the colors in the color cycle in case
        the cycle index is larger than the number of colors in the
        cycle.

        Parameters
        ----------
        cdx : int
            The cycle index of the color to return.

        Returns
        -------
        str
            The color at the given index.
        """
        cdx = cdx % len(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        return plt.rcParams['axes.prop_cycle'].by_key()['color'][cdx]

    def get_marker(self, cdx) -> str:
        """
        Returns the marker for the given cycle index. This needs to
        correctly loop over the markers in the marker cycle in case
        the cycle index is larger than the number of markers in the
        cycle.

        Parameters
        ----------
        cdx : int
            The cycle index of the marker to return.

        Returns
        -------
        str
            The marker at the given index.
        """
        cdx = cdx % len(self._markers)
        return self._markers[cdx]
    
    @property
    def default_figsize(self):
        """
        Returns the value of the default_figsize attribute.

        Returns
        -------
        tuple
            The value of the default_figsize attribute.
        """
        return self._default_figsize

    @property
    def default_title(self) -> str:
        """
        Returns the value of the title attribute.

        Returns
        -------
        str
            The value of the title attribute.
        """
        return self._title
    
    @property
    def mark_pot(self) -> bool:
        """
        Returns the value of the mark_pot attribute.

        Returns
        -------
        bool
            The value of the mark_pot attribute.
        """
        return self._mark_pot

    @property
    def mark_pot_horizontal(self) -> bool:
        """
        Returns the value of the mark_pot_horizontal attribute.

        Returns
        -------
        bool
            The value of the mark_pot_horizontal attribute.
        """
        return self._mark_pot_horizontal
    
    @property
    def mark_preliminary(self) -> str:
        """
        Returns the value of the mark_preliminary attribute.

        Returns
        -------
        str
            The value of the mark_preliminary attribute.
        """
        return self._mark_preliminary
    
    @property
    def scilimits(self) -> tuple:
        """
        Returns the value of the scilimits attribute.

        Returns
        -------
        tuple
            The value of the scilimits attribute.
        """
        return self._scilimits

    @property
    def plot_kwargs(self) -> dict:
        """
        Returns the value of the plot_kwargs attribute.

        Returns
        -------
        dict
            The value of the plot_kwargs attribute.
        """
        return self._plot_kwargs

def load_style_from_toml(style_name, styles_toml_path=None) -> Style:
    """Load a Style object by name from a styles TOML file."""
    if styles_toml_path is None:
        # Default to spineplot/configurations/common/styles.toml.
        styles_toml_path = join(dirname(__file__), 'configurations', 'common', 'styles.toml')

    config = toml.load(styles_toml_path)
    style_entries = config.get('style', [])
    for entry in style_entries:
        if entry.get('name') == style_name:
            style = Style(**entry)

            # If the resolved style path is not usable, resolve relative to spineplot root.
            if not isfile(style.get_style()):
                style_sheet = entry.get('style_sheet', '')
                if style_sheet and not isabs(style_sheet):
                    spineplot_root = dirname(dirname(dirname(styles_toml_path)))
                    style._style = normpath(join(spineplot_root, style_sheet))
            return style

    available = ', '.join([entry.get('name', '<unnamed>') for entry in style_entries])
    raise ValueError(
        f"Style '{style_name}' was not found in '{styles_toml_path}'. "
        f"Available styles: {available}"
    )


def apply_style_to_plot(style_name, ax=None, fig=None, styles_toml_path=None) -> Style:
    """Read a named style from TOML and apply it to the active Matplotlib plot.

    Parameters
    ----------
    style_name : str
        Name of the style block in styles.toml.
    ax : matplotlib.axes.Axes, optional
        Axis to update with title metadata from the style.
    fig : matplotlib.figure.Figure, optional
        Figure to resize with the style default size.
    styles_toml_path : str, optional
        Custom path to styles.toml. Defaults to spineplot's common styles.

    Returns
    -------
    Style
        The loaded Style object.
    """
    style = load_style_from_toml(style_name, styles_toml_path=styles_toml_path)
    plt.style.use(style.get_style())

    if ax is None:
        ax = plt.gca()
    if fig is None and ax is not None:
        fig = ax.figure

    if fig is not None and style.default_figsize is not None:
        fig.set_size_inches(*style.default_figsize)

    if ax is not None and style.default_title is not None:
        ax.set_title(style.default_title)

    # Apply mark_preliminary if available
    if ax is not None and style.mark_preliminary is not None and mark_preliminary is not None:
        mark_preliminary(ax, style.mark_preliminary)

    return style


def make_styled_plot(style_name, figsize=None, figure_size=None, styles_toml_path=None):
    """Convenience: create a figure with style already applied.
    
    Use this to be sure the style is active BEFORE you plot.
    
    Returns
    -------
    fig, ax, style : tuple
        Matplotlib figure, axes, and the loaded Style object.
    
    Example
    -------
    fig, ax, style = make_styled_plot("basic_icarus")
    ax.plot([1, 2, 3], [1, 4, 2])
    plt.show()
    """
    style = load_style_from_toml(style_name, styles_toml_path=styles_toml_path)
    plt.style.use(style.get_style())

    # `figure_size` is an explicit alias for `figsize`.
    if figure_size is not None:
        figsize = figure_size
    
    if figsize is None and style.default_figsize is not None:
        figsize = style.default_figsize
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if style.default_title is not None:
        ax.set_title(style.default_title)
    
    # Apply mark_preliminary if available
    if style.mark_preliminary is not None and mark_preliminary is not None:
        mark_preliminary(ax, style.mark_preliminary)
    
    return fig, ax, style
