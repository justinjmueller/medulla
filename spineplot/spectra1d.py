import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from spectra import SpineSpectra
from style import Style
from variable import Variable

class SpineSpectra1D(SpineSpectra):
    """
    A class designed to encapsulate a single variable's spectrum for an
    ensemble of samples. This is a specialization of the SpineSpectra
    class that is intended to be used for 1D spectra. At its core, this
    is a simple histogram plot of a single variable.

    Attributes
    ----------
    _variable : Variable
        The Variable object for the spectrum.
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
    def __init__(self, variable, categories, colors, category_types) -> None:
        """
        Initializes the SpineSpectra1D.

        Parameters
        ----------
        variable : Variable
            The Variable object for the spectrum.
        categories : dict
            A dictionary of the categories for the spectrum. This
            serves as a map between the category label in the input
            TTree and the category label for the spectrum (and
            therefore what is shown in a single legend entry).
        colors : dict
            A dictionary of the colors for the categories in the
            spectrum. This serves as a map between the category label
            for the spectrum (value in the `_categories` dictionary)
            and the color to use for the histogram. The color can be
            any valid matplotlib color string or a cycle indicator
            (e.g. 'C0', 'C1', etc.).
        category_types : dict
            A dictionary of the types for the categories in the
            spectrum. This serves as a map between the category label
            for the spectrum (value in the `_categories` dictionary)
            and the type of plot to use for the histogram. The type
            should be either 'histogram' or 'scatter' to correspond to
            a stacked histogram or scatter plot, respectively.

        Returns
        -------
        None.
        """
        super().__init__([variable,], categories, colors)
        self._variable = self._variables[0]
        self._category_types = category_types

    def add_sample(self, sample, is_ordinate) -> None:
        """
        Adds a sample to the SpineSpectra1D object. The sample's data
        is extracted per category and stored for later plotting.
        Multiple samples may have overlapping categories, so the data
        is stored in a dictionary with the category as the key.

        Parameters
        ----------
        sample : Sample
            The sample to add to the SpineSpectra1D object.
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
            self._onebincount = {}
        data, weights = sample.get_data([self._variable._key,], with_mask=self._variable.mask)
        for category, values in data.items():
            values = values[0]
            if category not in self._categories.keys():
                continue
            if self._categories[category] not in self._plotdata:
                self._plotdata[self._categories[category]] = np.zeros(self._variable._nbins)
                self._onebincount[self._categories[category]] = 0
            h = np.histogram(values, bins=self._variable._nbins, range=self._variable._range, weights=weights[category])
            self._onebincount[self._categories[category]] += np.sum(weights[category])
            self._plotdata[self._categories[category]] += h[0]
            self._binedges[self._categories[category]] = h[1]

    def draw(self, ax, style, override_xlabel=None, show_component_number=False,
             show_component_percentage=False, invert_stack_order=False,
             fit_type=None, logx=False, logy=False) -> None:
        """
        Plots the data for the SpineSpectra1D object.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to draw the artist on.
        style : Style
            The style to use when drawing the artist. The default is
            None. This is intended to be used in cases where the artist
            has some configurable style options.
        override_xlabel : str
            An optional override for the x-axis label. The default is
            None, which will use the label from the Variable object.
        show_component_number : bool
            A flag to indicate if the component number should be shown
            in the legend. The default is False.
        show_component_percentage : bool
            A flag to indicate if the component percentage should be
            shown in the legend. The default is False.
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

        Returns
        -------
        None.
        """
        ax.set_xlabel(self._variable._xlabel if override_xlabel is None else override_xlabel)
        ax.set_ylabel('Candidates')
        ax.set_xlim(*self._variable._range)

        if self._plotdata is not None:
            labels, data = zip(*self._plotdata.items())
            colors = [self._colors[label] for label in labels]
            bincenters = [self._binedges[l][:-1] + np.diff(self._binedges[l]) / 2 for l in labels]

            histogram_mask = [li for li, label in enumerate(labels) if self._category_types[label] == 'histogram']
            scatter_mask = [li for li, label in enumerate(labels) if self._category_types[label] == 'scatter']

            denominator = np.sum([self._onebincount[labels[i]] for i in histogram_mask])
            counts = [x for x in self._onebincount.values()]

            if fit_type is not None:
                super().fit_with_function(ax, bincenters[0], np.sum(data, axis=0), self._binedges[labels[0]], fit_type, range=self._variable._range)

            if show_component_number and show_component_percentage:
                hlabel = lambda x : f'{np.sum(x):.1f}, {np.sum(x)/denominator:.2%}'
                slabel = lambda x : f'{np.sum(x):.1f}'
                labels = [f'{label} ({hlabel(d) if li in histogram_mask else slabel(d)})' for li, (label, d) in enumerate(zip(labels, counts))]
            elif show_component_number:
                labels = [f'{label} ({np.sum(d):.1f})' for label, d in zip(labels, counts)]
            elif show_component_percentage:
                labels = [f'{label} ({np.sum(d)/denominator:.2%})' if li in histogram_mask else label for li, (label, d) in enumerate(zip(labels, counts))]

            if invert_stack_order:
                reduce = lambda x : [x[i] for i in histogram_mask[::-1]]
            else:
                reduce = lambda x : [x[i] for i in histogram_mask]
            
            ax.hist(reduce(bincenters), weights=reduce(data), bins=self._variable._nbins, range=self._variable._range, label=reduce(labels), color=reduce(colors), **style.plot_kwargs)

            reduce = lambda x : [x[i] for i in scatter_mask]
            for i, label in enumerate(reduce(labels)):
                ax.errorbar(bincenters[scatter_mask[i]], data[scatter_mask[i]], yerr=np.sqrt(data[scatter_mask[i]]), fmt='o', label=label, color=colors[scatter_mask[i]])
        
        if invert_stack_order:
            h, l = ax.get_legend_handles_labels()
            ax.legend(h[::-1], l[::-1])
        else:
            ax.legend()
        if style.mark_pot:
            self.mark_pot(ax)
        if style.mark_preliminary is not None:
            self.mark_preliminary(ax, style.mark_preliminary)

        # Set the axis to be logarithmic if requested.
        if logx:
            # Modify the x-axis limits to ensure that the lower limit
            # is greater than zero. The lower edge needs to be at least
            # 3 orders of magnitude less than the maximum value in the
            # plot.
            if self._variable._range[0] == 0:    
                xhigh_exporder = np.floor(np.log10(self._variable._range[1]))
                xlow = xhigh_exporder - 3
                ax.set_xlim(10**xlow, self._variable._range[1])
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')
        
        if style.get_title() is not None:
            ax.set_title(style.get_title())