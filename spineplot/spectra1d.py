import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as _Rect
import warnings

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

        # Per-plot (per-variable) area normalization support.
        # We track each sample's total weighted yield for this artist's variable
        # so we can optionally rescale the final histogram stack to match the
        # ordinate sample area at draw time.
        if not hasattr(self, "_sample_areas") or self._sample_areas is None:
            self._sample_areas = {}  # {sample_name: total_weighted_yield}
            self._sample_plotdata = {}  # {sample_name: {category_label: hist}}
            self._ordinate_sample_name = None
            self._ordinate_area = None
        if is_ordinate:
            self._ordinate_sample_name = sample._name

        if self._plotdata is None:
            self._plotdata = {}
            self._binedges = {}
            self._onebincount = {}
        data, weights = sample.get_data([self._variable._key,], with_mask=self._variable.mask)
        sample_area = 0.0
        if sample._name not in self._sample_plotdata:
            self._sample_plotdata[sample._name] = {}
        for category, values in data.items():
            values = values[0]
            if category not in self._categories.keys():
                continue
            if self._categories[category] not in self._plotdata:
                self._plotdata[self._categories[category]] = np.zeros(self._variable._nbins)
                self._onebincount[self._categories[category]] = 0
            xr = self._variable._range if self._xrange is None else self._xrange
            # Use either uniform bins or customized bins
            h = np.histogram(values, 
                bins=self._variable._nbins if self._variable._custom_bins is None else self._variable._custom_bins ,
                range=xr, weights=weights[category])
            # Total (in-range) weighted yield for this sample/variable.
            # This is the “plot area” used by exposure_type='area' normalization.
            # Note: this corresponds to the integral/sum of bin contents, not an
            # integral over x with bin-width factors.
            sample_area += float(np.sum(h[0]))
            self._onebincount[self._categories[category]] = self._onebincount.get(self._categories[category], 0.0) + float(np.sum(weights[category]))
            self._plotdata[self._categories[category]] += h[0]
            # Also keep per-sample histograms so exposure_type='area' can rescale
            # non-ordinate samples to match the ordinate integral (per plot).
            if self._categories[category] not in self._sample_plotdata[sample._name]:
                self._sample_plotdata[sample._name][self._categories[category]] = np.zeros_like(h[0])
            self._sample_plotdata[sample._name][self._categories[category]] += h[0]
            self._binedges[self._categories[category]] = h[1]


        # Record this sample's total weighted yield for this variable.
        self._sample_areas[sample._name] = self._sample_areas.get(sample._name, 0.0) + sample_area
        if is_ordinate:
            self._ordinate_area = self._sample_areas[sample._name]
    def draw(self, ax, style, show_component_number=False,
             show_component_percentage=False, invert_stack_order=False,
             fit_type=None, logx=False, logy=False, normalize=False,
             draw_error=None, make_ratio_plot=False, ratio_scatter_label='Data', ratio_ylim=(0.5, 1.5), ratio_ylabel='Data/MC') -> None:
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

        # Optional ratio subplot (Data/MC) below the main spectrum.
        # We create an additional axes within the same figure using the existing
        # axis position. This avoids needing figure-level plumbing changes.
        ratio_ax = None
        if make_ratio_plot:
            fig = ax.figure
            pos = ax.get_position()
            # Split the original axis vertically: top for spectrum, bottom for ratio.
            gap = 0.02
            ratio_frac = 0.28
            ratio_h = pos.height * ratio_frac
            main_h = pos.height - ratio_h - gap
            # Resize the main axis and create the ratio axis below it.
            ax.set_position([pos.x0, pos.y0 + ratio_h + gap, pos.width, main_h])
            ratio_ax = fig.add_axes([pos.x0, pos.y0, pos.width, ratio_h], sharex=ax)
            # Clean up tick labels on the top axis; ratio axis owns x-labels.
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.set_xlabel('')
            ratio_ax.set_ylabel(ratio_ylabel)
            ratio_ax.set_ylim(*ratio_ylim)
            ratio_ax.grid(True, axis='y', alpha=0.3)
            ratio_ax.set_xlabel(self._variable._xlabel if self._xtitle is None else self._xtitle)

        # Track the total stacked histogram and its uncertainty envelope so we can
        # expand y-limits and avoid clipping error boxes.
        total_y = None
        total_yerr = None

        if self._plotdata is not None:
            def _legend_count(label, payload, label_type):
                """Compute a sensible legend count for a component.

                Priority:
                  1) `_onebincount[label]` when it exists.
                  2) Fallback to the plotted payload sum.

                Notes:
                  - Histogram payloads are binned arrays, so sum is yield.
                  - Scatter payloads in this code are also binned arrays.
                """
                """Compute a sensible legend count for a component.

                Priority:
                  1) `_onebincount[label]` when it exists.
                  2) Fallback to the plotted payload sum.

                Notes:
                  - Histogram payloads are binned arrays, so sum is yield.
                  - Scatter payloads in this code are also binned arrays.
                """
                if hasattr(self, '_onebincount') and self._onebincount is not None:
                    if label in self._onebincount:
                        return float(self._onebincount.get(label, 0.0))
                try:
                    arr = np.asarray(payload, dtype=float)
                    return float(_np.nansum(arr))
                except Exception:
                    return 0.0

            def _legend_count_scaled(payload):
                """Legend count computed from the *currently plotted* payload.

                This is used for exposure_type='area' so the legend reflects post-scaling yields.
                """
                try:
                    arr = np.asarray(payload, dtype=float)
                    return float(np.nansum(arr))
                except Exception:
                    return 0.0

            labels, data = zip(*self._plotdata.items())
            colors = [self._colors[label] for label in labels]
            bincenters = [self._binedges[l][:-1] + np.diff(self._binedges[l]) / 2 for l in labels]
            binwidths = [np.diff(self._binedges[l]) for l in labels]
            xr = self._variable._range if self._xrange is None else self._xrange

            histogram_mask = [li for li, label in enumerate(labels) if self._category_types[label] == 'histogram']
            scatter_mask = [li for li, label in enumerate(labels) if self._category_types[label] == 'scatter']
            
            denominator = np.sum([self._onebincount[labels[i]] for i in histogram_mask])
            counts = [_legend_count(label, payload, self._category_types[label]) for (label, payload) in zip(labels, data)]

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

            hist_bincenters = reduce(bincenters)
            hist_data = [scale*x for x in reduce(data)]
            hist_labels = reduce(labels)
            hist_colors = reduce(colors)

            # Total stacked histogram (used for uncertainty band + ylim envelope)
            total_y = scale * np.sum(reduce(data), axis=0)
            # Optional per-plot area normalization to the ordinate sample.
            # Triggered when the ordinate sample uses exposure_type == 'area'.
            #
            # Behavior: keep the ordinate sample unchanged, and rescale all other samples
            # so that their *total stack area* matches the ordinate integral for this plot.
            #
            # This requires per-sample histogram bookkeeping (filled in add_sample).
            area_scales = {}
            if getattr(self, '_exposure_type', None) == 'area' and getattr(self, '_ordinate_sample_name', None) is not None:
                ord_name = self._ordinate_sample_name
                ord_area = float(self._sample_areas.get(ord_name, 0.0))
                for sname, sarea in getattr(self, '_sample_areas', {}).items():
                    sarea = float(sarea)
                    if sname == ord_name:
                        area_scales[sname] = 1.0
                    elif sarea > 0 and ord_area > 0:
                        area_scales[sname] = ord_area / sarea
                    else:
                        area_scales[sname] = 1.0

            if area_scales:
                # Rebuild the merged category histograms from the per-sample ones after scaling.
                # Keep only histogram categories (scatter are drawn separately below).
                scaled_plotdata = {k: np.zeros(self._variable._nbins) for k in self._plotdata.keys()}
                for sname, cat_hists in getattr(self, '_sample_plotdata', {}).items():
                    sf = area_scales.get(sname, 1.0)
                    for cat, hvals in cat_hists.items():
                        if cat in scaled_plotdata:
                            scaled_plotdata[cat] += sf * hvals

                # Replace the merged histograms and recompute derived arrays.
                self._plotdata.update(scaled_plotdata)
                labels, data = zip(*self._plotdata.items())
                # After area scaling, keep undecorated category labels for lookups
                # (colors, bin edges, etc.), and build decorated labels only for the legend.
                base_labels = list(labels)

                histogram_mask = [li for li, label in enumerate(base_labels) if self._category_types[label] == 'histogram']
                scatter_mask = [li for li, label in enumerate(base_labels) if self._category_types[label] == 'scatter']
                denominator = np.sum([self._onebincount[base_labels[i]] for i in histogram_mask])
                counts = [_legend_count_scaled(payload) for payload in data]

                decorated_labels = base_labels
                if show_component_number and show_component_percentage:
                    hlabel = lambda x : f'{np.sum(x):.1f}, {np.sum(x)/denominator:.2%}'
                    slabel = lambda x : f'{np.sum(x):.1f}'
                    decorated_labels = [f'{label} ({hlabel(d) if li in histogram_mask else slabel(d)})'
                                      for li, (label, d) in enumerate(zip(base_labels, counts))]
                elif show_component_number:
                    decorated_labels = [f'{label} ({np.sum(d):.1f})' for label, d in zip(base_labels, counts)]
                elif show_component_percentage:
                    decorated_labels = [f'{label} ({np.sum(d)/denominator:.2%})' if li in histogram_mask else label
                                      for li, (label, d) in enumerate(zip(base_labels, counts))]

                # Use base labels for style/data lookups
                labels = base_labels
                colors = [self._colors[label] for label in base_labels]
                bincenters = [self._binedges[l][:-1] + np.diff(self._binedges[l]) / 2 for l in base_labels]
                binwidths = [np.diff(self._binedges[l]) for l in base_labels]

                if invert_stack_order:
                    reduce = lambda x : [x[i] for i in histogram_mask[::-1]]
                else:
                    reduce = lambda x : [x[i] for i in histogram_mask]

                scale = 1.0 if not normalize else 1.0 / np.sum(reduce(data))
                hist_bincenters = reduce(bincenters)
                hist_data = [scale*x for x in reduce(data)]
                hist_labels = reduce(decorated_labels)
                hist_colors = reduce(colors)
                total_y = scale * np.sum(reduce(data), axis=0)

                # Ensure scatter plots also get decorated legend labels in area mode.
                scatter_labels = [decorated_labels[i] for i in scatter_mask]

            # Use either uniform bins or customized bins
            ax.hist(hist_bincenters, weights=hist_data, 
                    bins=self._variable._nbins if self._variable._custom_bins is None else self._variable._custom_bins,
                    range=xr, label=hist_labels, color=hist_colors, **style.plot_kwargs)
            if draw_error:
                systs = [s[draw_error] for s in self._systematics.values() if draw_error in s]
                if systs:
                    # Sum covariances across samples for this systematic
                    cov = np.sum(s.get_covariance(self._variable._key) for s in systs)
                    # If exposure_type='area' is active, also scale the uncertainty to match
                    # the scaled stack: cov_total = sum_i (f_i^2 * cov_i).
                    if area_scales:
                        cov_scaled = np.zeros_like(cov)
                        for sname, sdict in self._systematics.items():
                            if draw_error not in sdict:
                                continue
                            sf = float(area_scales.get(sname, 1.0))
                            cov_scaled += (sf**2) * sdict[draw_error].get_covariance(self._variable._key)
                        cov = cov_scaled
                    x = reduce(bincenters)[0]
                    y = total_y
                    xerr = [bw / 2 for bw in binwidths[0]]
                    scov = Systematic.transform_as(cov, scale if not normalize else np.sum(reduce(data), axis=0))
                    yerr = np.sqrt(np.diag(scov))
                    total_yerr = yerr
                    draw_error_boxes(ax, x, y, xerr, yerr, facecolor='gray', edgecolor='none', alpha=0.5, hatch='///')
                else:
                    available = []
                    if self._systematics:
                        try:
                            available = list(next(iter(self._systematics.values())).keys())
                        except Exception:
                            available = []
                    warnings.warn(
                        f"draw_error='{draw_error}' not found in any sample's systematics. "
                        f"Available keys: {available}"
                    )
                    draw_error = None   # suppress legend entry

            scatter_labels = [labels[i] for i in scatter_mask] if 'scatter_labels' not in locals() else scatter_labels
            reduce = lambda x : [x[i] for i in scatter_mask]
            for i, label in enumerate(scatter_labels):
                scale = 1.0 if not normalize else 1.0 / np.sum(data[scatter_mask[i]])
                ax.errorbar(bincenters[scatter_mask[i]], scale*data[scatter_mask[i]], yerr=scale*np.sqrt(data[scatter_mask[i]]), fmt='o', label=label, color=colors[scatter_mask[i]])
        

            # Draw ratio panel if requested (first scatter component / total stack).
            if ratio_ax is not None and scatter_mask and total_y is not None:
                # Choose which scatter series to ratio (defaults to label 'Data').
                si0 = None
                try:
                    target = ratio_scatter_label
                    # First try exact match on the base labels.
                    for _si in scatter_mask:
                        if labels[_si] == target:
                            si0 = _si
                            break
                    # If not found and labels are decorated (e.g. 'Data (123)'), allow prefix match.
                    if si0 is None:
                        for _si in scatter_mask:
                            if str(labels[_si]).startswith(str(target)):
                                si0 = _si
                                break
                except Exception:
                    si0 = None
                if si0 is None:
                    si0 = scatter_mask[0]
                # Scatter payload is stored binned (one y per bin center).
                y_data = np.asarray(data[si0], dtype=float)
                # Apply the same normalization used for scatter drawing above.
                sscale = 1.0 if not normalize else 1.0 / np.sum(data[si0])
                y_data = sscale * y_data
                yerr_data = sscale * np.sqrt(np.clip(y_data / max(sscale, 1e-12), 0.0, None))

                y_mc = np.asarray(total_y, dtype=float)
                with np.errstate(divide='ignore', invalid='ignore'):
                    ratio = np.where(y_mc > 0, y_data / y_mc, np.nan)
                    ratio_err = np.where(y_mc > 0, yerr_data / y_mc, np.nan)

                x = bincenters[si0]
                ratio_ax.errorbar(x, ratio, yerr=ratio_err, fmt='o', color=colors[si0], markersize=4)

                # Optional MC uncertainty band in ratio (if available from draw_error).
                if total_yerr is not None:
                    ymcerr = np.asarray(total_yerr, dtype=float)
                    with np.errstate(divide='ignore', invalid='ignore'):
                        band = np.where(y_mc > 0, ymcerr / y_mc, np.nan)
                    # Use bin edges so the uncertainty band spans each bin (centered at x)
                    edges = self._binedges[labels[si0]]
                    # Turn edges (N+1) and band (N) into step-style arrays for filled rectangles
                    x_step = np.repeat(edges, 2)[1:-1]          # length 2*N
                    band_step = np.repeat(band, 2)             # length 2*N
                    ratio_ax.fill_between(x_step, 1.0 - band_step, 1.0 + band_step,
                                         color='gray', alpha=0.35, linewidth=0)
                    # Use the same label as the main-axis error legend when available,
                    # otherwise fall back to a generic description.
                    label = 'MC Uncertainty'
                    proxy = _Rect((0, 0), 1, 1, fc='gray', alpha=0.35, linewidth=0)
                    ratio_ax.legend([proxy], [label])
                ratio_ax.axhline(1.0, color='k', linewidth=1, alpha=0.6)

                # If top axis is logy, ratio panel should remain linear.
                ratio_ax.set_yscale('linear')
        if invert_stack_order:
            h, l = ax.get_legend_handles_labels()
            if draw_error:
                h.append(plt.Rectangle((0, 0), 1, 1, fc='gray', alpha=0.5, hatch='///'))
                l.append(systs[0].label or systs[0].name)
                ax.legend(h[-2::-1]+h[-1:], l[-2::-1]+l[-1:])
            else:
                ax.legend(h[::-1], l[::-1])
        else:
            h, l = ax.get_legend_handles_labels()
            if draw_error:
                h.append(plt.Rectangle((0, 0), 1, 1, fc='gray', alpha=0.5, hatch='///'))
                l.append(systs[0].label or systs[0].name)
            ax.legend(h, l)

        if isinstance(self._yrange, (tuple, list)):
            ax.set_ylim(*self._yrange)
        elif isinstance(self._yrange, (int, float)):
            yl = ax.get_ylim()[1]
            ax.set_ylim(None, yl * self._yrange)
        else:
            # If the user didn't explicitly set a y-range, expand the upper limit
            # to include the uncertainty envelope (y + yerr) so error boxes won't
            # get clipped.
            if total_y is not None:
                current_low, current_high = ax.get_ylim()

                y_upper = (total_y + total_yerr) if total_yerr is not None else total_y

                # Ignore non-finite values.
                env_max = np.nanmax(np.where(np.isfinite(y_upper), y_upper, np.nan))
                if np.isfinite(env_max):
                    if logy:
                        # For log scale, ylim must be positive. Also add multiplicative headroom.
                        low = current_low
                        if low <= 0:
                            pos = total_y[total_y > 0] if total_y is not None else np.array([])
                            low = np.nanmin(pos) if pos.size else 1e-6
                        high = max(current_high, env_max) * 1.3
                        ax.set_ylim(low, high)
                    else:
                        # Add a bit of linear headroom.
                        high = max(current_high, env_max) * 1.15
                        ax.set_ylim(current_low, high)

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

class SpineSpectraCutFlow(SpineSpectra):
    """
    An artist that shows the distribution of a single variable after
    each of an ordered sequence of cuts, drawn as overlaid step
    histograms (one per cut) on a shared axis. Each histogram is the
    *total* event count (all categories summed), coloured by cut
    position so the effect of each selection step is visible at a
    glance.

    Attributes
    ----------
    _cuts : dict
        Ordered mapping of cut branch name → human-readable cut label.
        The branch name must exist as a boolean (0/1) column in the
        sample data. Cuts are applied cumulatively in order: the first
        entry is the distribution with only that cut applied, the
        second with the first AND second cut applied, etc.
        A special key ``'no_cut'`` (or any key whose branch evaluates
        to all-True) can be used as the pre-selection baseline.
    _plotdata : dict
        {cut_label: np.ndarray of bin counts} accumulated across
        all registered samples.
    _binedges : np.ndarray
        Shared bin edges (same for every cut).
    """

    def __init__(self, variable, categories, colors, category_types,
                 cuts, title=None, xrange=None, xtitle=None,
                 yrange=None, ytitle=None) -> None:
        """
        Parameters
        ----------
        variable : Variable
            The variable whose distribution is shown.
        categories : dict
            Category-key → category-label mapping (passed through to
            the base class; used only for exposure tracking here).
        colors : dict
            Category-label → colour mapping (not used for line colours,
            which are taken from the style cycle, but kept for API
            compatibility).
        category_types : dict
            Unused by this artist but kept for API compatibility.
        cuts : dict
            Ordered mapping of cut branch name → cut label. Cuts are
            applied cumulatively in the listed order.
        title : str, optional
        xrange : tuple, optional
        xtitle : str, optional
        yrange : tuple or float, optional
        ytitle : str, optional
        """
        super().__init__([variable], categories, colors, title,
                         xrange, xtitle, yrange, ytitle)
        self._variable = self._variables[0]
        self._category_types = category_types
        # Store cuts as an ordered list of (branch, label) pairs so
        # iteration order is deterministic on all Python versions.
        self._cuts = list(cuts.items())   # [(branch, label), ...]
        self._plotdata = None             # {label: np.ndarray}
        self._binedges = None

    # ------------------------------------------------------------------
    def add_sample(self, sample, is_ordinate) -> None:
        """Accumulate per-cut histograms from *sample*."""
        super().add_sample(sample, is_ordinate)

        if self._plotdata is None:
            self._plotdata = {label: np.zeros(self._variable._nbins)
                              for _, label in self._cuts}
            self._binedges = None

        xr = self._variable._range if self._xrange is None else self._xrange
        bins = (self._variable._nbins
                if self._variable._custom_bins is None
                else self._variable._custom_bins)

        # Build the cumulative mask expression progressively.
        cumulative_branches = []
        for branch, label in self._cuts:
            cumulative_branches.append(branch)

            # Combine masks: each cut column is 1.0/0.0 in the data.
            # Build a pandas eval expression that ANDs all cuts so far.
            if len(cumulative_branches) == 1:
                mask_expr = f'{branch} == 1'
            else:
                mask_expr = ' & '.join(
                    f'({b} == 1)' for b in cumulative_branches)

            # get_data applies the mask and returns per-category arrays.
            data, weights = sample.get_data(
                [self._variable._key], with_mask=mask_expr)

            total_counts = np.zeros(self._variable._nbins)
            total_weights = []
            all_values = []

            for category, values in data.items():
                if category not in self._categories:
                    continue
                vals = values[0]
                w = weights[category]
                all_values.append(vals)
                total_weights.append(w)

            if all_values:
                combined_vals = np.concatenate(all_values)
                combined_w = np.concatenate(total_weights)
                h, edges = np.histogram(combined_vals, bins=bins,
                                        range=xr, weights=combined_w)
                total_counts += h
                if self._binedges is None:
                    self._binedges = edges

            self._plotdata[label] += total_counts

    # ------------------------------------------------------------------
    def draw(self, ax, style, logx=False, logy=False,
             normalize=False, show_counts=True) -> None:
        """
        Draw overlaid step histograms — one per cut — on *ax*.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
        style : Style
        logx : bool
        logy : bool
        normalize : bool
            Normalise each histogram to unit area independently.
        show_counts : bool
            Append the total event count to each legend label.
        """
        ax.set_xlabel(
            self._variable._xlabel if self._xtitle is None else self._xtitle)
        ax.set_ylabel('Candidates' if not normalize else 'A.U.')
        xr = self._variable._range if self._xrange is None else self._xrange
        ax.set_xlim(*xr)
        ax.set_title(self._title)

        # Optional ratio subplot (Data/MC) below the main spectrum.
        # We create an additional axes within the same figure using the existing
        # axis position. This avoids needing figure-level plumbing changes.
        ratio_ax = None
        if make_ratio_plot:
            fig = ax.figure
            pos = ax.get_position()
            # Split the original axis vertically: top for spectrum, bottom for ratio.
            gap = 0.02
            ratio_frac = 0.28
            ratio_h = pos.height * ratio_frac
            main_h = pos.height - ratio_h - gap
            # Resize the main axis and create the ratio axis below it.
            ax.set_position([pos.x0, pos.y0 + ratio_h + gap, pos.width, main_h])
            ratio_ax = fig.add_axes([pos.x0, pos.y0, pos.width, ratio_h], sharex=ax)
            # Clean up tick labels on the top axis; ratio axis owns x-labels.
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.set_xlabel('')
            ratio_ax.set_ylabel(ratio_ylabel)
            ratio_ax.set_ylim(*ratio_ylim)
            ratio_ax.grid(True, axis='y', alpha=0.3)

        if self._plotdata is None or self._binedges is None:
            return

        bins = (self._variable._nbins
                if self._variable._custom_bins is None
                else self._variable._custom_bins)
        bin_edges = self._binedges
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2

        # Use the style colour cycle for cut lines.
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = [c['color'] for c in prop_cycle]

        for ci, (_, label) in enumerate(self._cuts):
            counts = self._plotdata[label].copy()
            scale = 1.0 / np.sum(counts) if (normalize and np.sum(counts) > 0) else 1.0
            counts = counts * scale

            legend_label = label
            if show_counts:
                legend_label = f'{label} ({np.sum(self._plotdata[label]):.1f})'

            color = colors[ci % len(colors)]
            ax.stairs(counts, bin_edges, label=legend_label,
                      color=color, linewidth=1.8)

        ax.legend()

        # y-range
        if isinstance(self._yrange, (tuple, list)):
            ax.set_ylim(*self._yrange)
        elif isinstance(self._yrange, (int, float)):
            ax.set_ylim(None, ax.get_ylim()[1] * self._yrange)

        if logx:
            if xr[0] == 0:
                xhigh = np.floor(np.log10(xr[1]))
                ax.set_xlim(10 ** (xhigh - 3), xr[1])
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')

        hadj, vadj = 0, 0
        if style.scilimits and not logy:
            ax.ticklabel_format(axis='y', scilimits=style.scilimits)
            hadj = 0.035
        if logy:
            vadj = 0.1

        if style.mark_pot:
            mark_pot(ax, self._exposure, style.mark_pot_horizontal, vadj=vadj)
        if style.mark_preliminary is not None:
            mark_preliminary(ax, style.mark_preliminary, hadj=hadj, vadj=vadj)


class SpineSystematics(SpineSpectra):
    """
    An artist that plots the fractional (or absolute) systematic
    uncertainty as a function of a single variable, one band per
    systematic recipe. Each recipe's uncertainty is the per-bin
    fractional error derived from the diagonal of the covariance
    matrix built from the multisim/multisigma universes stored in
    the ordinate sample.

    The plot can be drawn in two modes (controlled by *mode* in
    ``draw_kwargs``):

    ``'fractional'`` (default)
        y-axis shows ``sigma_i / N_i`` per bin, where ``N_i`` is the
        central-value bin count.  Bins with zero counts are masked.

    ``'absolute'``
        y-axis shows ``sigma_i`` directly.

    TOML usage example::

        [[figure]]
        name = "sys_eleE"
        type = "SimpleFigure"
        style = 'basic_icarus'
        [[figure.artists]]
        type        = "SpineSystematics"
        variable    = "reco_leading_electron_energy"
        # Optional list of recipe names to include (defaults to all recipes).
        # systematics = ["xsec_ccqe", "xsec_mec"]
        [figure.artists.draw_kwargs]
        mode        = "fractional"   # or "absolute"
        show_total  = true           # overlay quadrature sum of shown recipes
        show_stat   = false          # include the statistical uncertainty
        logy        = false

    Attributes
    ----------
    _variable : Variable
        The variable whose systematic uncertainty is shown.
    _recipe_names : list[str] or None
        Names of recipes to include.  If None all available recipes are used.
    _covdata : dict
        {recipe_name: np.ndarray covariance matrix} collected from the ordinate sample.
    _cv : np.ndarray or None
        Central-value histogram (all categories summed).
    _binedges : np.ndarray or None
        Bin edges for the variable.
    _recipe_labels : dict
        {recipe_name: human-readable label} for legend entries.
    """

    def __init__(self, variable, categories, colors, category_types,
                 recipe_names=None,
                 all_recipe_names=None,
                 title=None, xrange=None, xtitle=None,
                 yrange=None, ytitle=None) -> None:
        super().__init__([variable], categories, colors, title,
                         xrange, xtitle, yrange, ytitle)
        self._variable       = self._variables[0]
        self._category_types = category_types
        # _allowed_names: the set of sysnames to draw.
        # Priority: explicit 'systematics' list > all recipe names > everything.
        if recipe_names is not None:
            self._allowed_names = set(recipe_names)
        elif all_recipe_names is not None:
            # Default: show only the recipe-combined entries, not raw branches
            self._allowed_names = set(all_recipe_names)
        else:
            self._allowed_names = None   # truly no filter — show everything
        self._covdata        = {}             # {recipe_name: cov_matrix}
        self._recipe_labels  = {}             # {recipe_name: label}
        self._cv             = None
        self._binedges       = None

    # ------------------------------------------------------------------
    def add_sample(self, sample, is_ordinate) -> None:
        """Collect covariance matrices and central-value histogram from the ordinate sample."""
        super().add_sample(sample, is_ordinate)

        # Only the ordinate sample is used for systematics
        if not is_ordinate:
            return

        vkey = self._variable._key
        bins = (self._variable._nbins
                if self._variable._custom_bins is None
                else self._variable._custom_bins)
        xr   = self._variable._range if self._xrange is None else self._xrange

        # Build central-value histogram (all categories summed)
        data, weights = sample.get_data([vkey], with_mask=self._variable.mask)
        cv        = np.zeros(self._variable._nbins)
        bin_edges = None
        for category, values in data.items():
            if category not in self._categories:
                continue
            vals = values[0]
            w    = weights[category]
            h, edges = np.histogram(vals, bins=bins, range=xr, weights=w)
            cv += h
            if bin_edges is None:
                bin_edges = edges
        self._cv       = cv
        self._binedges = bin_edges

        # Collect covariance matrices from the sample's processed systematics
        for sysname, syst in sample._systematics.items():
            # Filter: only include names in the allowed set (recipes or explicit list)
            if self._allowed_names is not None and sysname not in self._allowed_names:
                continue
            cov_key = f'{sysname}_{vkey}'
            if cov_key not in syst._covariances:
                continue
            cov = syst._covariances[cov_key].copy()
            if sysname in self._covdata:
                self._covdata[sysname] += cov
            else:
                self._covdata[sysname] = cov
            # Store human-readable label for legend (use recipe label if set)
            self._recipe_labels[sysname] = syst._label if syst._label else sysname

    # ------------------------------------------------------------------
    def draw(self, ax, style,
             mode='fractional',
             show_total=True,
             show_stat=False,
             logy=False) -> None:
        """
        Draw one step-line per systematic recipe.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
        style : Style
        mode : str
            ``'fractional'`` (default) or ``'absolute'``.
        show_total : bool
            Overlay the quadrature sum of all shown recipes as a dashed black line.
        show_stat : bool
            Include the statistical uncertainty in the plot and total.
        logy : bool
            Use a logarithmic y-axis.
        """
        if self._cv is None or self._binedges is None:
            return

        ax.set_xlabel(
            self._variable._xlabel if self._xtitle is None else self._xtitle)
        ax.set_ylabel('Fractional uncertainty' if mode == 'fractional'
                      else 'Absolute uncertainty')
        xr = self._variable._range if self._xrange is None else self._xrange
        ax.set_xlim(*xr)
        ax.set_title(self._title)

        # Optional ratio subplot (Data/MC) below the main spectrum.
        # We create an additional axes within the same figure using the existing
        # axis position. This avoids needing figure-level plumbing changes.
        ratio_ax = None
        if make_ratio_plot:
            fig = ax.figure
            pos = ax.get_position()
            # Split the original axis vertically: top for spectrum, bottom for ratio.
            gap = 0.02
            ratio_frac = 0.28
            ratio_h = pos.height * ratio_frac
            main_h = pos.height - ratio_h - gap
            # Resize the main axis and create the ratio axis below it.
            ax.set_position([pos.x0, pos.y0 + ratio_h + gap, pos.width, main_h])
            ratio_ax = fig.add_axes([pos.x0, pos.y0, pos.width, ratio_h], sharex=ax)
            # Clean up tick labels on the top axis; ratio axis owns x-labels.
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.set_xlabel('')
            ratio_ax.set_ylabel(ratio_ylabel)
            ratio_ax.set_ylim(*ratio_ylim)
            ratio_ax.grid(True, axis='y', alpha=0.3)

        bin_edges = self._binedges
        cv        = self._cv
        safe_cv   = np.where(cv > 0, cv, np.nan)   # avoid division by zero

        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors     = [c['color'] for c in prop_cycle]

        total_var = np.zeros(len(cv))
        plotted   = 0

        for sysname, cov in self._covdata.items():
            # Skip statistical uncertainty unless explicitly requested
            if not show_stat and 'statistical' in sysname.lower():
                continue
            sigma = np.sqrt(np.diag(cov))
            y     = sigma / safe_cv if mode == 'fractional' else sigma
            label = self._recipe_labels.get(sysname, sysname)
            color = colors[plotted % len(colors)]
            ax.stairs(y, bin_edges, label=label, color=color, linewidth=1.8)
            total_var += np.diag(cov)
            plotted   += 1

        if show_total and plotted > 1:
            sigma_total = np.sqrt(total_var)
            y_total     = sigma_total / safe_cv if mode == 'fractional' else sigma_total
            ax.stairs(y_total, bin_edges, label='Total', color='black',
                      linewidth=2.2, linestyle='--')

        ax.legend()

        if isinstance(self._yrange, (tuple, list)):
            ax.set_ylim(*self._yrange)
        elif isinstance(self._yrange, (int, float)):
            ax.set_ylim(None, ax.get_ylim()[1] * self._yrange)

        if logy:
            ax.set_yscale('log')

        hadj, vadj = 0, 0
        if style.scilimits and not logy:
            ax.ticklabel_format(axis='y', scilimits=style.scilimits)
            hadj = 0.035

        if style.mark_pot:
            mark_pot(ax, self._exposure, style.mark_pot_horizontal, vadj=vadj)
        if style.mark_preliminary is not None:
            mark_preliminary(ax, style.mark_preliminary, hadj=hadj, vadj=vadj)
