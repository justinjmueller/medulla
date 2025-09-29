import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from spectra import SpineSpectra
from style import Style
from variable import Variable
from systematic import Systematic
from utilities import mark_pot, mark_preliminary, draw_error_boxes

class SpineSpectra1D(SpineSpectra):
    """
    A class designed to encapsulate a single variable's spectrum for an
    ensemble of samples. This is a specialization of the SpineSpectra
    class that is intended to be used for 1D spectra. At its core, this
    is a simple histogram plot of a single variable.

    Attributes
    ----------
    _title : str
        The title of the spectrum. This will be placed at the top of
        the axis assigned to the artist.
    _xrange : tuple
        The range of the x-axis for the spectrum. This is a tuple of
        the lower and upper limits of the x-axis. If None, the range
        will be determined by the range set in the Variable object.
    _xtitle : str
        The label for the x-axis. If None, the label will be taken
        from the Variable object.
    _yrange : tuple, or float, optional
        If this is a tuple, it is the range of the y-axis for the
        spectrum. If this is a float, it will scale the maximum value
        of the histogram by this factor. If None, the range will be
        determined by the range of the histogram.
    _ytitle : str
        The label for the y-axis. If None, the label will be 'Candidates'.
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
    def __init__(self, variable, categories, colors, category_types,
                 title=None, xrange=None, xtitle=None,
                 yrange=None, ytitle=None) -> None:
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
        title : str, optional
            The title of the spectrum. This will be placed at the top
            of the axis assigned to the artist. The default is None.
        xrange : tuple, optional
            The range of the x-axis for the spectrum. This is a tuple
            of the lower and upper limits of the x-axis. If None, the
            range will be determined by the range set in the Variable
            object. The default is None.
        xtitle : str, optional
            The label for the x-axis. If None, the label will be taken
            from the Variable object. The default is None.
        yrange : tuple, or float, optional
            If this is a tuple, it is the range of the y-axis for the
            spectrum. If this is a float, it will scale the maximum
            value of the histogram by this factor. If None, the range
            will be determined by the range of the histogram. The
            default is None.
        ytitle : str, optional
            The label for the y-axis. If None, the label will be
            'Candidates'. The default is None.

        Returns
        -------
        None.
        """
        super().__init__([variable,], categories, colors, title,
                         xrange, xtitle, yrange, ytitle)
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
            xr = self._variable._range if self._xrange is None else self._xrange
            h = np.histogram(values, bins=self._variable._nbins, range=xr, weights=weights[category])
            self._onebincount[self._categories[category]] += np.sum(weights[category])
            self._plotdata[self._categories[category]] += h[0]
            self._binedges[self._categories[category]] = h[1]

    def draw(self, ax, style, show_component_number=False,
             show_component_percentage=False, invert_stack_order=False,
             fit_type=None, logx=False, logy=False, normalize=False,
             draw_error=None) -> None:
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
        normalize : bool
            A flag to indicate if the histogram should be normalized
            to unity. The default is False.
        draw_error : str, optional
            Indicates the name of the Systematic object to use for
            drawing the error boxes. The default is None.

        Returns
        -------
        None.
        """
        ax.set_xlabel(self._variable._xlabel if self._xtitle is None else self._xtitle)
        ax.set_ylabel('Candidates')
        ax.set_xlim(*self._variable._range if self._xrange is None else self._xrange)
        ax.set_title(self._title)

        if self._plotdata is not None:
            labels, data = zip(*self._plotdata.items())
            colors = [self._colors[label] for label in labels]
            bincenters = [self._binedges[l][:-1] + np.diff(self._binedges[l]) / 2 for l in labels]
            binwidths = [np.diff(self._binedges[l]) for l in labels]
            xr = self._variable._range if self._xrange is None else self._xrange

            histogram_mask = [li for li, label in enumerate(labels) if self._category_types[label] == 'histogram']
            scatter_mask = [li for li, label in enumerate(labels) if self._category_types[label] == 'scatter']

            denominator = np.sum([self._onebincount[labels[i]] for i in histogram_mask])
            counts = [x for x in self._onebincount.values()]

            if fit_type is not None:
                super().fit_with_function(ax, bincenters[0], np.sum(data, axis=0), self._binedges[labels[0]], fit_type, range=xr)

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
            
            scale = 1.0 if not normalize else 1.0 / np.sum(reduce(data))
            ax.hist(reduce(bincenters), weights=[scale*x for x in reduce(data)], bins=self._variable._nbins, range=xr, label=reduce(labels), color=reduce(colors), **style.plot_kwargs)
            if draw_error:
                systs = [s[draw_error] for s in self._systematics.values() if draw_error in s]
                cov = np.sum(s.get_covariance(self._variable._key) for s in systs)
                x = reduce(bincenters)[0]
                y = scale * np.sum(reduce(data), axis=0)
                xerr = [x / 2 for x in binwidths[0]]
                scov = Systematic.transform_as(cov, scale if not normalize else np.sum(reduce(data), axis=0))
                yerr = np.sqrt(np.diag(scov))
                draw_error_boxes(ax, x, y, xerr, yerr, facecolor='gray', edgecolor='none', alpha=0.5, hatch='///')

            reduce = lambda x : [x[i] for i in scatter_mask]
            for i, label in enumerate(reduce(labels)):
                scale = 1.0 if not normalize else 1.0 / np.sum(data[scatter_mask[i]])
                ax.errorbar(bincenters[scatter_mask[i]], scale*data[scatter_mask[i]], yerr=scale*np.sqrt(data[scatter_mask[i]]), fmt='o', label=label, color=colors[scatter_mask[i]])
        
        if invert_stack_order:
            h, l = ax.get_legend_handles_labels()
            if draw_error:
                h.append(plt.Rectangle((0, 0), 1, 1, fc='gray', alpha=0.5, hatch='///'))
                l.append(systs[0].label)
                ax.legend(h[-2::-1]+h[-1:], l[-2::-1]+l[-1:])
            else:
                ax.legend(h[::-1], l[::-1])
        else:
            h, l = ax.get_legend_handles_labels()
            if draw_error:
                h.append(plt.Rectangle((0, 0), 1, 1, fc='gray', alpha=0.5, hatch='///'))
                l.append(systs[0].label)
            ax.legend(h, l)

        if isinstance(self._yrange, (tuple, list)):
            ax.set_ylim(*self._yrange)
        elif isinstance(self._yrange, (int, float)):
            yl = ax.get_ylim()[1]
            ax.set_ylim(None, yl * self._yrange)

        # Set the axis to be logarithmic if requested.
        if logx:
            # Modify the x-axis limits to ensure that the lower limit
            # is greater than zero. The lower edge needs to be at least
            # 3 orders of magnitude less than the maximum value in the
            # plot.
            xr = self._variable._range if self._xrange is None else self._xrange
            if xr[0] == 0:    
                xhigh_exporder = np.floor(np.log10(xr[1]))
                xlow = xhigh_exporder - 3
                ax.set_xlim(10**xlow, xr[1])
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')

        # hadj and vadj are used to adjust the position of the POT and
        # preliminary labels horizontally and vertically, respectively.
        # This is necessary to ensure that the labels do not overlap
        # with plot elements. The following logic is meant to capture
        # all cases where the labels might overlap with the plot.
        hadj = 0
        vadj = 0

        # The scilimits option cannot be used with a logarithmic y-axis.
        # The hadj value is adjusted to ensure that the POT label does
        # not overlap with the scientific notation placed above the
        # y-axis.
        if style.scilimits and not logy:
            ax.ticklabel_format(axis='y', scilimits=style.scilimits)
            hadj = 0.035

        # The vadj value is adjusted to ensure that the POT label does
        # not overlap with the top axis of the plot when the y-axis is
        # logarithmic.
        if logy:
            vadj = 0.1
        
        # Add the POT and preliminary labels to the plot.
        if style.mark_pot:
            mark_pot(ax, self._exposure, style.mark_pot_horizontal, vadj=vadj)
        if style.mark_preliminary is not None:
            mark_preliminary(ax, style.mark_preliminary, hadj=hadj, vadj=vadj)