import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from spectra import SpineSpectra
from style import Style
from variable import Variable
from utilities import mark_pot, mark_preliminary, draw_error_boxes

class SpineSpectra2D(SpineSpectra):
    """
    A class designed to encapsulate a pair of variables' spectrum for
    an ensemble of samples. The class method add_sample() can be used
    to add a sample to the SpineSpectra2D.

    Attributes
    ----------
    _title : str
        The title of the artist. This will be placed at the top of the
        axis assigned to the artist.
    _xrange : tuple
        The range of the x-axis for the spectrum. This is a tuple of
        the minimum and maximum values for the x-axis. If None, the
        range will be determined by Variable object assigned to the
        x-axis (show_option='2d') or set to (-1,1)
        (show_option='projection').
    _xtitle : str
        The label for the x-axis of the spectrum. If None, the label
        will be determined by the Variable object assigned to the
        x-axis (show_option='2d') or set to '(Y-X)/X'
        (show_option='projection').
    _yrange : tuple
        The range of the y-axis for the spectrum. This is a tuple of
        the minimum and maximum values for the y-axis. If None, the
        range will be determined by Variable object assigned to the
        y-axis (show_option='2d') or set to None
        (show_option='projection').
    _ytitle : str
        The label for the y-axis of the spectrum. If None, the label
        will be determined by the Variable object assigned to the
        y-axis (show_option='2d') or set to 'Entries'
        (show_option='projection').
    _variables : list
        The list of Variable objects for the spectrum.
    _categories : dict
        A dictionary of the categories for the spectrum. This serves as
        a map between the category label in the input TTree and the
        category label for the spectrum (and therefore what is shown
        in a single legend entry).
    _colors : dict
        A dictionary of the colors for the categories in the spectrum.
        This serves as a map between the category label for the
        spectrum (value in the `_categories` dictionary) and the color
        to use for the histogram. The color can be any valid matplotlib
        color string or a cycle indicator (e.g. 'C0', 'C1', etc.).
    _plotdata : dict
        A dictionary of the data for the spectrum. This is a map between
        the category label for the spectrum and the histogram data for
        that category.
    """
    def __init__(self, variables, categories, colors, category_types,
                 title=None, xrange=None, xtitle=None, yrange=None,
                 ytitle=None) -> None:
        """
        Initializes the SpineSpectra2D object.

        Parameters
        ----------
        variables : list
            The list of Variable objects for the spectrum.
        categories : dict
            A dictionary of the categories for the spectrum. This serves
            as a map between the category label in the input TTree and
            the category label for the spectrum (and therefore what is
            shown in a single legend entry).
        colors : dict
            A dictionary of the colors for the categories in the spectrum.
            This serves as a map between the category label for the
            spectrum (value in the `_categories` dictionary) and the color
            to use for the histogram. The color can be any valid matplotlib
            color string or a cycle indicator (e.g. 'C0', 'C1', etc.).
        category_types : dict
            A dictionary of the types for the categories in the spectrum.
            This serves as a map between the category label for the spectrum
            (value in the `_categories` dictionary) and the type of plot to
            use for the histogram. The type should be either 'histogram' or
            'scatter' to correspond to a stacked histogram or scatter plot,
            respectively.
        title : str, optional
            The title of the artist. This will be placed at the top of
            the axis assigned to the artist. The default is None.
        xrange : tuple, optional
            The range of the x-axis for the spectrum. This is a tuple of
            the minimum and maximum values for the x-axis. If None, the
            range will be determined by Variable object assigned to the
            x-axis (show_option='2d') or set to (-1,1)
            (show_option='projection'). The default is None.
        xtitle : str, optional
            The label for the x-axis of the spectrum. If None, the label
            will be determined by the Variable object assigned to the
            x-axis (show_option='2d') or set to '(Y-X)/X'
            (show_option='projection'). The default is None.
        yrange : tuple, optional
            The range of the y-axis for the spectrum. This is a tuple of
            the minimum and maximum values for the y-axis. If None, the
            range will be determined by Variable object assigned to the
            y-axis (show_option='2d') or set to None
            (show_option='projection'). The default is None.
        ytitle : str, optional
            The label for the y-axis of the spectrum. If None, the label
            will be determined by the Variable object assigned to the
            y-axis (show_option='2d') or set to 'Entries'
            (show_option='projection'). The default is None.

        Returns
        -------
        None.
        """
        super().__init__(variables, categories, colors, title,
                         xrange, xtitle, yrange, ytitle)
        self._category_types = category_types
        self._plotdata_diagonal = None
        self._binedges_diagonal = None
        self._plotdata_absdiag = None
        self._binedges_absdiag = None

    def add_sample(self, sample, is_ordinate) -> None:
        """
        Adds a sample to the SpineSpectra2D object. The sample's data
        is extracted per category and stored for later plotting.
        Multiple samples may have overlapping categories, so the data
        is stored in a dictionary with the category as the key.

        Parameters
        ----------
        sample : Sample
            The sample to add to the SpineSpectra2D object.
        is_ordinate : bool
            A flag to indicate if the sample is the ordinate sample.

        Returns
        -------
        None.
        """
        super().add_sample(sample, is_ordinate)

        if self._plotdata is None:
            self._plotdata = {}
            self._binedges = {}
        if self._plotdata_diagonal is None:
            self._plotdata_diagonal = {}
            self._binedges_diagonal = {}
        if self._plotdata_absdiag is None:
            self._plotdata_absdiag = {}
            self._binedges_absdiag = {}
        

        # Check if a mask is present for the variables. If so, we need
        # to combine the masks for the two variables.
        if self._variables[0].mask is not None and self._variables[1].mask is not None:
            joint_mask = f'{self._variables[0].mask} and {self._variables[1].mask}'
        elif self._variables[0].mask is not None:
            joint_mask = self._variables[0].mask
        elif self._variables[1].mask is not None:
            joint_mask = self._variables[1].mask
        else:
            joint_mask = None

        data, weights = sample.get_data([self._variables[0]._key, self._variables[1]._key], joint_mask)        
        for category, values in data.items():
            if category not in self._categories.keys():
                continue
            if self._categories[category] not in self._plotdata:
                self._plotdata[self._categories[category]] = np.zeros((self._variables[0]._nbins, self._variables[1]._nbins))
            xr = self._variables[0]._range if self._xrange is None else self._xrange
            yr = self._variables[1]._range if self._yrange is None else self._yrange
            h = np.histogram2d(values[0], values[1], bins=(self._variables[0]._nbins, self._variables[1]._nbins), range=(xr, yr), weights=weights[category])
            self._plotdata[self._categories[category]] += h[0]
            self._binedges[self._categories[category]] = h[1]

            if self._categories[category] not in self._plotdata_diagonal:
                self._plotdata_diagonal[self._categories[category]] = np.zeros(self._variables[0]._nbins)
            diag = np.divide(values[1] - values[0], values[0])
            xr = (-1, 1) if self._xrange is None else self._xrange
            h = np.histogram(diag, bins=self._variables[0]._nbins, range=xr, weights=weights[category])
            self._plotdata_diagonal[self._categories[category]] += h[0]
            self._binedges_diagonal[self._categories[category]] = h[1]

            if self._categories[category] not in self._plotdata_absdiag:
                self._plotdata_absdiag[self._categories[category]] = np.zeros(self._variables[0]._nbins)
            absdiag = values[1] - values[0]
            h = np.histogram(absdiag, bins=self._variables[0]._nbins, range=(-100,100), weights=weights[category])
            self._plotdata_absdiag[self._categories[category]] += h[0]
            self._binedges_absdiag[self._categories[category]] = h[1]

    def draw(self, ax, style, show_option='2d', draw_identity=True,
             draw_colorbar=True, invert_stack_order=False,
             fit_type=None, logx=False, logy=False, logz=False,
             draw_stat_error=False) -> None:
        """
        Plots the data for the SpineSpectra2D object.

        Parameters
        ----------
        ax : Axes
            The matplotlib Axes object to draw the plot on.
        style : Style
            The Style object to use for the plot.
        show_option : str
            The option to use for the plot. This can be one of a few
            options (default is '2d'):
                '2d'         - Draw a 2D histogram of the data.
                'projection' - Draw a projection of the data about the
                               diagonal
                'absdiff'    - Draw a absolute difference of the data
        draw_identity : bool
            A flag to indicate if the identity line should be drawn on
            the plot. The default is True.
        invert_stack_order : bool
            A flag to indicate if the stack order in the legend should
            be inverted. The default is False.
        fit_type : str
            The type of fit to perform on the data. The default is
            None, which will not perform any fit. The options are:
                'crystal_ball' - Perform a Crystal Ball fit on the data.
                'gaussian'     - Perform a Gaussian fit on the data.
        logx : bool
            A flag to indicate if the x-axis should be logarithmic.
            The default is False.
        logy : bool
            A flag to indicate if the y-axis should be logarithmic.
            The default is False.
        logz : bool
            A flag to indicate if the z-axis (colorbar) should be
            logarithmic. The default is False.
        draw_stat_error : bool
            A flag to indicate if the statistical error should be drawn
            on the plot. The default is False.
        
        Returns
        -------
        None.
        """
        ax.set_title(self._title)
        
        if show_option == '2d' and self._plotdata is not None:
            values = np.sum([v for v in self._plotdata.values()], axis=0)
            binedges = self._binedges[list(self._plotdata.keys())[0]]
            ax.imshow(values.T, extent=(binedges[0], binedges[-1], binedges[0], binedges[-1]), aspect='auto', origin='lower')
            ax.set_xlabel(self._variables[0]._xlabel if self._xtitle is None else self._xtitle)
            ax.set_ylabel(self._variables[1]._xlabel)
            ax.set_aspect('equal')
            
            # Draw the identity line. This must span the full range
            # of the plot, so we need to find the minimum and maximum
            # of the range for the plot.
            if draw_identity:
                min_range = min([binedges[0], binedges[-1]])
                max_range = max([binedges[0], binedges[-1]])
                ax.plot([min_range, max_range], [min_range, max_range], 'k--')

        if show_option == 'projection' and self._plotdata_diagonal is not None:
            labels, data = zip(*self._plotdata_diagonal.items())
            colors = [self._colors[label] for label in labels]
            bincenters = [self._binedges_diagonal[l][:-1] + np.diff(self._binedges_diagonal[l]) / 2 for l in labels]

            ax.hist(bincenters, weights=data, bins=self._variables[0]._nbins,
                    range=(-1,1) if self._xrange is None else self._xrange,
                    histtype='barstacked', label=labels, color=colors, stacked=True)
            ax.set_xlabel('(Y-X)/X' if self._xtitle is None else self._xtitle)
            ax.set_ylabel('Entries')
            ax.set_xlim(-1, 1) if self._xrange is None else ax.set_xlim(self._xrange)

            if fit_type is not None:
                super().fit_with_function(ax, bincenters[0], np.sum(data, axis=0), self._binedges_diagonal[labels[0]], fit_type, range=(-1,1) if self._xrange is None else self._xrange)

            if draw_stat_error:
                x = bincenters[0]
                y = np.sum(data, axis=0)
                xerr = np.diff(self._binedges_diagonal[labels[0]]) / 2
                yerr = np.sqrt(y)
                draw_error_boxes(ax, x, y, xerr, yerr, facecolor='gray', edgecolor='none', alpha=0.5, hatch='///')

            if invert_stack_order:
                h, l = ax.get_legend_handles_labels()
                if draw_stat_error:
                    h.append(plt.Rectangle((0, 0), 1, 1, fc='gray', alpha=0.5, hatch='///'))
                    l.append('MC Statistical Uncertainty')
                ax.legend(h[-2::-1]+h[-1:], l[-2::-1]+l[-1:])
            else:
                h, l = ax.get_legend_handles_labels()
                if draw_stat_error:
                    h.append(plt.Rectangle((0, 0), 1, 1, fc='gray', alpha=0.5, hatch='///'))
                    l.append('MC Statistical Uncertainty')
                ax.legend(h, l)
        
        if style.scilimits and not logy:
            ax.ticklabel_format(axis='y', scilimits=style.scilimits)
        if style.mark_pot:
            mark_pot(ax, self._exposure, style.mark_pot_horizontal)
        if style.mark_preliminary is not None:
            mark_preliminary(ax, style.mark_preliminary)
        # Set the axis to be logarithmic if requested.
        if logx:
            # Modify the x-axis limits to ensure that the lower limit
            # is greater than zero. The lower edge needs to be at least
            # 3 orders of magnitude less than the maximum value in the
            # plot.
            xr = self._variables[0]._range if self._xrange is None else self._xrange
            if xr[0] == 0:    
                xhigh_exporder = np.floor(np.log10(xr[1]))
                xlow = xhigh_exporder - 3
                ax.set_xlim(10**xlow, xr[1])
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')
