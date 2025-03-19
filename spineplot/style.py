import matplotlib.pyplot as plt

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
        self._style = style_sheet
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