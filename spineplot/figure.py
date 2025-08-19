import matplotlib.pyplot as plt

from style import Style

class SpineFigure:
    """
    A class to represent a single figure in a plot. This class serves
    as a base class for any figure that will be plotted with the
    package. All other figure classes will inherit from this class;
    this class will not be used directly and is only meant to provide
    common functionality to all other figure classes.

    Attributes
    ----------
    _figure : matplotlib.figure.Figure
        The parent figure object that will be used to contain all of
        the subplots.
    _figsize : tuple
        The size of the figure to create.
    _style : Style
        The style to use when drawing the figure.
    _title : str
        The title of the figure.
    _axs : list
        A list of axes objects that will be used to plot the figure.
    _artists : list
        A list of SpineArtists to be drawn on the figure.
    _draw_kwargs : list
        A list of keyword arguments (dicts) to pass to the draw method
        of the SpineArtists.

    """
    def __init__(self, figsize, style, title=None):
        """
        Parameters
        ----------
        figsize : tuple
            The size of the figure to create.
        style : Style
            The style to use when drawing the figure.
        title : str, optional
            The title of the figure. The default is None.
        """
        self._figure = None
        self._figsize = figsize
        self._style = style
        self._title = title
        self._axs = []
        self._artists = []
        self._draw_kwargs = []

    def register_spine_artist(self, artist, draw_kwargs):
        """
        Register an artist with the figure. This method is used to
        add an artist to the figure so that it can be displayed when
        the figure is shown.

        Parameters
        ----------
        artist : matplotlib.artist.Artist
            The artist to add to the figure.
        draw_kwargs : dict
            A dictionary of keyword arguments to pass to the draw
            method of the artist.
        """
        self._artists.append(artist)
        self._draw_kwargs.append(draw_kwargs)

    def create(self):
        """
        Create the figure. This method is used to create the figure
        and all of the subplots that will be used to display the data.

        Returns
        -------
        None.
        """
        with self._style as style:
            for axi, ax in enumerate(self._axs):
                self._artists[axi].draw(ax, **self._draw_kwargs[axi], style=style)
            self._figure.suptitle(self._title, y=0.99)
    
    def close(self):
        """
        Close the figure. This method is used to close the figure
        and release the resources associated with it.

        Returns
        -------
        None.
        """
        plt.close(self._figure)

    @property
    def figure(self):
        """
        Returns
        -------
        matplotlib.figure.Figure
            The parent figure object that will be used to contain all
            of the subplots.
        """
        return self._figure

    @property
    def axs(self):
        """
        Returns
        -------
        list
            A list of axes objects that will be used to plot the figure.
        """
        return self._axs

class SimpleFigure(SpineFigure):
    """
    A simple figure with a single axis. This class is used to create
    a single plot with a single set of data. The class is a subclass
    of SpineFigure.
    """
    def __init__(self, figsize, style, title=None):
        """
        Parameters
        ----------
        style : Style
            The style to use when drawing the figure.
        figsize : tuple, optional
            The size of the figure to create. The default is (8, 6).
        """
        super().__init__(figsize=figsize, style=style, title=title)

    def create(self):
        """
        Create the figure. This method is used to create the figure
        and all of the subplots that will be used to display the data.

        Returns
        -------
        None.
        """
        with self._style as style:
            self._figure = plt.figure(figsize=self._figsize)
            self._axs = [self._figure.add_subplot(),]
            super().create()
