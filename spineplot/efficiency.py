import numpy as np
from scipy.stats import binom, beta
import pandas as pd
import matplotlib.pyplot as plt

from artists import SpineArtist
from style import Style
from variable import Variable
from utilities import mark_pot, mark_preliminary

class SpineEfficiency(SpineArtist):
    """
    A class designed to encapsulate the calculation of the efficiency
    of the selection as a function of the specified variable. The method
    employed is a Bayesian approach, where the true efficiency is assumed
    to be a binomial random variable with a uniform prior. The posterior
    distribution is then calculated using the binomial likelihood and the
    prior. 

    Attributes
    ----------
    _title : str
        The title of the artist. This will be placed at the top of the
        axis assigned to the artist.
    _xrange : tuple
        The range of the x-axis. If None, the range is taken from the
        Variable object.
    _xtitle : str
        The title of the x-axis. If None, the title is taken from the
        Variable object.
    _variable : Variable
        The variable to calculate the efficiency with respect to.
    _samples : list
        A list of samples to use in the calculation of the efficiency.
    _categories : dict
        A dictionary mapping the category key to the category name.
    _cuts : dict
        A dictionary mapping the cut key to the cut label.
    _show_option : str
        The option to use when showing the artist.
    _npts : int
        The number of points to use when calculating the efficiency.
    _posteriors : dict
        A dictionary containing the posterior distributions for the
        efficiency calculation.
    _totals : dict
        A dictionary containing the total number of events in each bin
        of the variable.
    _successes : dict
        A dictionary containing the number of successful events in each
        bin of the variable.
    """
    def __init__(self, variable, categories, cuts, title,
                 xrange=None, xtitle=None, show_option='table',
                 npts=100, show_counts=True):
        """
        Parameters
        ----------
        variable : Variable
            The variable to calculate the efficiency with respect to.
        categories : dict
            A dictionary mapping the category key to the category name.
        cuts : dict
            A dictionary mapping the cut key to the cut label.
        title : str
            The title of the artist.
        xrange : tuple, optional
            The range of the x-axis. If None, the range is taken from
            the Variable object. The default is None.
        xtitle : str, optional
            The title of the x-axis. If None, the title is taken from
            the Variable object. The default is None.
        show_option : str, optional
            The option to use when showing the artist. The default is
            'table.'
        npts : int, optional
            The number of points to use when calculating the efficiency.
            The default is 1e6.
        show_counts : bool, optional
            Toggle whether numerator/denominator counts are appended to table entries. The default is True.
        """
        super().__init__(title)
        self._xrange = xrange
        self._xtitle = xtitle
        self._variable = variable
        self._samples = list()
        self._categories = categories
        self._cuts = cuts
        self._show_option = show_option
        self._show_counts = show_counts
        self._npts = int(npts)
        self._posteriors = dict()
        self._totals = dict()
        self._successes = dict()
        self._selected_counts = dict()
        self._nuniv = 1000
        self._seed = 12345

    def draw(self, ax, show_option, percentage=True, show_seqeff=True,
             show_unseqeff=True, show_purity=False, yrange=None, npts=100, style=None,
             logx=False, logy=False, colors=None,
             draw_error=None, show_syst_band=True, syst_alpha=0.25,
             syst_color='gray', syst_label=None, show_legend=True):
        """
        Draw the artist on the given axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to draw the artist on.
        show_option : str, optional
            The option to use when showing the artist. The default is
            None. This is intended to be used in cases where the artist
            can be shown in different ways (e.g. 2D vs projection of 2D
            down to 1D).
        percentage : bool, optional
            A flag to indicate if the efficiency should be displayed
            as a percentage. The default is True.
        show_seqeff : bool, optional
            A flag to indicate if the sequential (cumulative)
            efficiency should be shown. The default is True.
        show_unseqeff : bool, optional
            A flag to indicate if the unsequential (single-cut)
            efficiency should be shown. The default is True.
        show_purity : bool, optional
            A flag to indicate if the purity should be shown in the
            efficiency table. The default is False.
        yrange : tuple, optional
            The range of the y-axis. The default is None.
        npts : int, optional
            The number of points to use when calculating the efficiency.
            The default is 1e6.
        style : Style, optional
            The style to use when drawing the artist. The default is
            None. This is intended to be used in cases where the artist
            has some configurable style options.
        logx : bool, optional
            A flag to indicate if the x-axis should be logarithmic.
            The default is False.
        logy : bool, optional
            A flag to indicate if the y-axis should be logarithmic.
            The default is False.

        Returns
        -------
        None.
        """
        ax.set_title(self._title)
        if style is None:
            style = Style()
        self._draw_error = draw_error

        # Get a list of the groups, while respecting the order of the
        # groups as they are configured in the analysis block. It is
        # possible that duplicates exist (multiple categories within
        # the same group), so care must be taken to ensure that the
        # groups are unique AND in the correct order.
        groups = list()
        for category in self._categories.values():
            if category not in groups:
                groups.append(category)

        if show_option == 'table':
            # Lambda formatter to round the values to two decimal
            # places, display as percentage (if toggled), and add super
            # and subscripts for the error values.
            if percentage:
                formatter = lambda x,y,z: rf'${100*x:.2f}^{{\ +{100*y:.2f}}}_{{\ -{100*z:.2f}}}$'
                diff_key = 'Differential\nEfficiency [%]'
                cumu_key = 'Cumulative\nEfficiency [%]'
                purity_key = 'Purity [%]'
            else:
                formatter = lambda x,y,z: rf'${x:.2f}^{{\ +{y:.2e}}}_{{\ -{z:.2e}}}$'
                diff_key = 'Differential\nEfficiency'
                cumu_key = 'Cumulative\nEfficiency'
                purity_key = 'Purity'

            # Clear up the axis because we are going to draw a table
            # on it (no need for any other plot elements).
            ax.axis('off')

            # Create the table data.
            results = pd.DataFrame({r'   ': list(), r'Cut': list(), r'Efficiency [%]': list(), r'Cumulative [%]': list()})
            group_endpoint = dict()

            # Loop over the groups requested in the plot.
            for group in groups:
                # Extract the cv, msigma, and psigma values for the
                # group.
                _, cv, msigma, psigma = self.reduce(group, significance=0.6827)
                _, pur_cv, pur_msigma, pur_psigma = self.calculate_purity(group, significance=0.6827)

                # Build per-cut formatted strings including numerator/denominator
                cuts_order = list(self._cuts.keys())
                diff_list = []
                cumu_list = []
                purity_list = []

                # precompute total selected across all groups for purity denominators
                total_selected_all = {f'seq_{c}': 0 for c in self._cuts.keys()}
                for g, cats in self._selected_counts.items():
                    for cat, cuts_dict in cats.items():
                        for k, v in cuts_dict.items():
                            if k.startswith('seq_'):
                                total_selected_all[k] = total_selected_all.get(k, 0) + int(v)

                # gather per-category counts for this group
                group_cats = self._selected_counts.get(group, {})
                total_events = sum(cat.get('total', 0) for cat in group_cats.values())

                for cut in cuts_order:
                    unseq_key = f'unbinned_unseq_{cut}'
                    seq_key = f'unbinned_seq_{cut}'

                    u_cv = cv.get(unseq_key, 0.0)
                    u_m = msigma.get(unseq_key, 0.0)
                    u_p = psigma.get(unseq_key, 0.0)
                    s_cv = cv.get(seq_key, 0.0)
                    s_m = msigma.get(seq_key, 0.0)
                    s_p = psigma.get(seq_key, 0.0)

                    n_unseq = sum(cat.get(f'unseq_{cut}', 0) for cat in group_cats.values())
                    n_seq = sum(cat.get(f'seq_{cut}', 0) for cat in group_cats.values())

                    if self._show_counts:
                        if total_events > 0:
                            unseq_str = f"{formatter(u_cv, u_m, u_p)} ({int(n_unseq)}/{int(total_events)})"
                            seq_str = f"{formatter(s_cv, s_m, s_p)} ({int(n_seq)}/{int(total_events)})"
                        else:
                            unseq_str = f"{formatter(u_cv, u_m, u_p)} (0/0)"
                            seq_str = f"{formatter(s_cv, s_m, s_p)} (0/0)"
                    else:
                        unseq_str = f"{formatter(u_cv, u_m, u_p)}"
                        seq_str = f"{formatter(s_cv, s_m, s_p)}"

                    diff_list.append(unseq_str)
                    cumu_list.append(seq_str)

                    if show_purity:
                        pur_key = f'purity_seq_{cut}'
                        p_cv = pur_cv.get(pur_key, 0.0)
                        p_m = pur_msigma.get(pur_key, 0.0)
                        p_p = pur_psigma.get(pur_key, 0.0)
                        n_group = sum(cat.get(f'seq_{cut}', 0) for cat in group_cats.values())
                        n_total = total_selected_all.get(f'seq_{cut}', 0)
                        if self._show_counts:
                            if n_total > 0:
                                purity_str = f"{formatter(p_cv, p_m, p_p)} ({int(n_group)}/{int(n_total)})"
                            else:
                                purity_str = f"{formatter(p_cv, p_m, p_p)} (0/0)"
                        else:
                            purity_str = f"{formatter(p_cv, p_m, p_p)}"
                        purity_list.append(purity_str)



                entry = {   r'   ': [group,] + [r'' for _ in range(1, len(self._cuts))],
                            r'Cut': self._cuts.values(),
                            diff_key: diff_list,
                            cumu_key: cumu_list }
                if show_purity:
                    entry[purity_key] = purity_list
                results = pd.concat([results, pd.DataFrame(entry)])
                group_endpoint[group] = len(results)

            # Trim the table to only show the columns that are
            # requested by the user.
            cols = [r'   ', r'Cut',] if len(groups) > 1 else [r'Cut',]
            if show_unseqeff:
                cols.append(diff_key)
            if show_seqeff:
                cols.append(cumu_key)
            if show_purity:
                cols.append(purity_key)
            results = results[cols]

            # Rename "Cumulative" to "Efficiency" if it is the only
            # efficiency to show.
            if show_seqeff and not show_unseqeff:
                results.rename(columns={cumu_key : 'Efficiency [%]' if percentage else 'Efficiency'}, inplace=True)
            
            table_data = [results.columns.to_list()] + results.values.tolist()
            table = ax.table(cellText=table_data, colLabels=None, loc='center', cellLoc='center', edges='T')
            table.scale(1, 2.75)
            for i in range(2, len(table_data)):
                if i == len(table_data) - 1:
                    for j in range(len(table_data[i])):
                        table[i, j].visible_edges = 'B'
                else:
                    for j in range(len(table_data[i])):
                        table[i, j].visible_edges = 'open'
                    if i in group_endpoint.values():
                        table[i, 0].visible_edges = 'B'

            def calc_bbox_yext(obj):
                figure = plt.gcf()
                bbox = obj.get_window_extent(renderer=figure.canvas.get_renderer())
                p0, p1 = figure.transFigure.inverted().transform(bbox)
                return p1[1] - p0[1]

            scale = 2.75
            while calc_bbox_yext(table) > 0.92:
                table.scale(1, 1/scale)
                scale -= 0.05
                table.scale(1, scale)
            if scale < 2.0:
                print(f'Warning: Table with title `{self._title}` is too large to fit on the figure (scale = {scale:.2f}). Consider extending the figure vertically.')

            # Mark the POT and preliminary information on the plot.
            if style.mark_pot:
                mark_pot(ax, self._exposure, style.mark_pot_horizontal, vadj=0.1)
            if style.mark_preliminary is not None:
                mark_preliminary(ax, style.mark_preliminary, vadj=0.1)

        elif show_option in ('differential', 'final_only'):
            # Lambda formatter to round the values to two decimal
            # places, and display as percentage if toggled.
            if percentage:
                fmt = lambda x : 100*x
            else:
                fmt = lambda x : x

            # For clarity, this artist will only draw either the
            # sequential or the non-sequential efficiency. The
            # configuration file technically has a toggle for both,
            # but if both are requested, the artist will only draw the
            # sequential efficiency. The default is to draw the
            # sequential efficiency.
            if show_seqeff:
                key_base = 'binned_seq_'
            elif show_unseqeff:
                key_base = 'binned_unseq_'
            else:
                key_base = 'binned_seq_'

            # Decide which cuts to draw:
            # - 'differential': all cuts
            # - 'final_only'  : only the last cut
            all_cuts = list(self._cuts.items())
            if show_option == 'final_only':
                cuts_to_draw = [all_cuts[-1]]
            else:
                cuts_to_draw = all_cuts

            # Note: if the user requests a differential plot, the
            # clarity of the plot is highly dependent on the number
            # of groups requested. If the number of groups is one,
            # it makes sense to denote the difference between cuts
            # with different colors. However, if the number of groups
            # is greater than one, it is better to color the groups
            # differently and denote the difference between cuts with
            # different marker styles. This can get a little messy,
            # but this gives the best "default" setting while leaving
            # the user the flexibility (and the responsibility) to 
            # ensure that the plot is clear and understandable.
            if len(groups) == 1:
                for ci, (cut, cutname) in enumerate(cuts_to_draw):
                    _, cv, msigma, psigma = self.reduce(groups[0], significance=0.6827)

                    # Optional systematic band (draw behind points).
                    if show_syst_band and draw_error is not None:
                        eff_raw = cv[f'{key_base}{cut}']
                        # print("Computing efficiency for cut:", cutname, cut)
                        sigma_eff = self.compute_syst_sigma_eff(groups[0], cut, sequential=(key_base == 'binned_seq_'))
                        if sigma_eff is not None:
                            _lab = syst_label if syst_label is not None else f'{draw_error} syst.'
                            ax.bar(
                                self._variable._bin_centers[groups[0]],
                                2 * fmt(sigma_eff),
                                width=self._variable._bin_widths[groups[0]],
                                bottom=fmt(eff_raw - sigma_eff),
                                color=colors.get(groups[0], style.get_color(0)) if colors else style.get_color(0),
                                alpha=syst_alpha,
                                zorder=0,
                                edgecolor='none',
                            )
                            # print(f"Debug: Adding systematic band for group {groups[0]}, cut {cut}, with sigma_eff = {sigma_eff}, raw_eff = {eff_raw}, cv = {cv[f'{key_base}{cut}']}")

                    ax.errorbar(
                        self._variable._bin_centers[groups[0]],
                        fmt(cv[f'{key_base}{cut}']),
                        xerr=self._variable._bin_widths[groups[0]] / 2.0,
                        yerr=[
                            fmt(msigma[f'{key_base}{cut}']),
                            fmt(psigma[f'{key_base}{cut}']),
                        ],
                        fmt='o',
                        color=colors.get(groups[0], style.get_color(ci)) if (colors and show_option == 'final_only') else style.get_color(ci),
                        label=cutname if show_option != 'final_only' else None,
                    )
            else:
                for gi, group in enumerate(groups):
                    for ci, (cut, cutname) in enumerate(cuts_to_draw):
                        _, cv, msigma, psigma = self.reduce(group, significance=0.6827)
                        # print("Computing efficiency for cut:", cutname, cut)

                        # Optional systematic band (draw behind points).
                        if show_syst_band and draw_error is not None:
                            eff_raw = cv[f'{key_base}{cut}']
                            sigma_eff = self.compute_syst_sigma_eff(group, cut, sequential=(key_base == 'binned_seq_'))
                            if sigma_eff is not None:
                                _lab = syst_label if syst_label is not None else f'{draw_error} syst.'
                                ax.bar(
                                    self._variable._bin_centers[group],
                                    2 * fmt(sigma_eff),
                                    width=self._variable._bin_widths[group],
                                    bottom=fmt(eff_raw - sigma_eff),
                                    color=colors.get(group, style.get_color(gi)) if colors else style.get_color(gi),
                                    alpha=syst_alpha,
                                    zorder=0,
                                    edgecolor='none',
                                )
                                # print(f"Debug: Adding systematic band for group {groups[gi]}, cut {cut}, with sigma_eff = {sigma_eff}, raw_eff = {eff_raw}, cv = {cv[f'{key_base}{cut}']}")

                        ax.errorbar(
                            self._variable._bin_centers[group],
                            fmt(cv[f'{key_base}{cut}']),
                            xerr=self._variable._bin_widths[group] / 2.0,
                            yerr=[
                                fmt(msigma[f'{key_base}{cut}']),
                                fmt(psigma[f'{key_base}{cut}']),
                            ],
                            fmt=style.get_marker(ci),
                            color=colors.get(group, style.get_color(gi)) if colors else style.get_color(gi),
                            label=group if show_option == 'final_only' else f'{group} : {cutname}',
                        )

            ax.set_xlabel(self._variable._xlabel if self._xtitle is None else self._xtitle)
            ax.set_ylabel('Efficiency [%]' if percentage else 'Efficiency')
            ax.set_xlim(self._variable._range if self._xrange is None else self._xrange)
            if yrange is not None:
                ax.set_ylim(yrange)

            # Set the axis to be logarithmic if requested.
            if logx:
                # Modify the x-axis limits to ensure that the lower limit
                # is greater than zero. The lower edge needs to be at least
                # 3 orders of magnitude less than the maximum value in the
                # plot.
                xr = self._variable._range if self._xrange is None else self._xrange
                if xr == 0:    
                    xhigh_exporder = np.floor(np.log10(xr[1]))
                    xlow = xhigh_exporder - 3
                    ax.set_xlim(10**xlow, xr[1])
                ax.set_xscale('log')
            if logy:
                ax.set_yscale('log')

            if show_legend:
                ax.legend()

            # Mark the POT and preliminary information on the plot.
            if style.scilimits:
                ax.ticklabel_format(axis='y', scilimits=style.scilimits)
            if style.mark_pot:
                mark_pot(ax, self._exposure, style.mark_pot_horizontal)
            if style.mark_preliminary is not None:
                mark_preliminary(ax, style.mark_preliminary, hadj=0.035 if style.scilimits is not None else 0)

    def add_sample(self, sample, is_ordinate):
        """
        Add a sample to the efficiency calculation. The calculation of 
        the efficiency uses a Bayesian approach, where the true
        efficiency is assumed to be a binomial random variable with a
        uniform prior. The posterior distribution is then calculated
        using the binomial likelihood and the prior. The nature of this
        calculation means that it can be performed on a sample-by-sample
        (and category-by-category) basis, which is why the sample is
        processed individually.

        Parameters
        ----------
        sample : Sample
            The sample to add to the efficiency calculation.
        is_ordinate : bool
            A flag to indicate if the sample is the ordinate sample.

        Returns
        -------
        None.
        """
        super().add_sample(sample, is_ordinate)
        self._samples.append(sample)

        # Calculate the efficiency for the sample
        self.calculate(sample)

    @staticmethod
    def multiply_posteriors(pos0, pos1):
        """
        Update the first posterior with the product of the two
        posterior distributions. This is a simple process, but it is
        important to normalize the posterior distribution after the
        multiplication. Additionally, the method is designed to handle
        both 1D and 2D posterior distributions (i.e. the method can
        handle both unbinned and binned efficiency calculations).

        Parameters
        ----------
        pos0 : np.ndarray
            The first posterior distribution.
        pos1 : np.ndarray
            The second posterior distribution.
        group0 : str
            The name of the first group, which is used to identify the
            posterior distribution if an exception is raised.
        group1 : str
            The name of the second group, which is used to identify the
            posterior distribution if an exception is raised.
        
        Returns
        -------
        pos : np.ndarray
            The product of the two posterior distributions.
        """
        pos = pos0*pos1
        if len(pos.shape) == 1:
            pos /= np.sum(pos)
        else:
            pos /= np.sum(pos, axis=-1)[:, np.newaxis]
        # Clip underflow noise to keep the distribution sparse
        pos[pos < 1e-10] = 0.0
        return pos

    def calculate(self, sample, significance=0.6827):
        """
        Calculate the efficiency of the selection on the sample as a
        function of the variable of interest. The method employed is a
        Bayesian approach, where the true efficiency is assumed to be a
        binomial random variable with a uniform prior. The posterior
        distribution is then calculated using the binomial likelihood
        and the prior. The nature of this calculation means that it can
        be performed on a sample-by-sample (and category-by-category)
        basis, which is why each category within the sample is processed
        individually.

        Parameters
        ----------
        sample : Sample
            The sample to calculate the efficiency for.
        significance : float, optional
            The significance level to use when calculating the
            efficiency. The default is 0.6827.

        Returns
        -------
        None.
        """
        # Create a linear space of efficiencies to calculate the
        # posterior distribution.
        efficiencies = np.linspace(0.0, 1, self._npts, dtype=np.float32)

        # Get the data for the binning variable and the configured
        # cuts. The data is returned as a dictionary with the key
        # being the category and the value consisting of the variable
        # and the cuts.
        data, _ = sample.get_data([self._variable._key, *self._cuts.keys()], with_mask=self._variable.mask)
        for category, values in data.items():
            if category not in self._categories:
                continue

            # If the group does not already have an entry in the
            # posteriors dictionary, create one.
            if self._categories[category] not in self._posteriors:
                self._posteriors[self._categories[category]] = {f'unbinned_seq_{c}' : np.ones(efficiencies.shape) for c in self._cuts.keys()}
                self._posteriors[self._categories[category]].update({f'unbinned_unseq_{c}' : np.ones(efficiencies.shape) for c in self._cuts.keys()})

                self._totals[self._categories[category]] = np.zeros(self._variable._nbins)

                # If the show option is set to 'differential', create
                # the necessary dictionaries to store the binned
                # efficiencies.
                if self._show_option in ('differential', 'final_only'):
                    
                    self._posteriors[self._categories[category]].update({f'binned_seq_{c}' : np.ones((self._variable._nbins, efficiencies.shape[0])) for c in self._cuts.keys()})
                    self._posteriors[self._categories[category]].update({f'binned_unseq_{c}' : np.ones((self._variable._nbins, efficiencies.shape[0])) for c in self._cuts.keys()})
                    
                    self._successes[self._categories[category]] = {f'binned_seq_{c}' : np.zeros(self._variable._nbins) for c in self._cuts.keys()}
                    self._successes[self._categories[category]].update({f'binned_unseq_{c}' : np.zeros(self._variable._nbins) for c in self._cuts.keys()})
            
            # Prepare bookkeeping for totals and per-cut selections used in tables/purity
            group_name = self._categories[category]
            if group_name not in self._selected_counts:
                self._selected_counts[group_name] = {}
            if category not in self._selected_counts[group_name]:
                init = {'total': 0}
                for _c in self._cuts.keys():
                    init[f'seq_{_c}'] = 0
                    init[f'unseq_{_c}'] = 0
                self._selected_counts[group_name][category] = init

            total_events = len(values[0])
            self._selected_counts[group_name][category]['total'] += int(total_events)

            # The calculation of efficiency as a function of some the
            # variable of interest requires the binning of the variable
            # into a number of bins.
            bin_edges = self._variable._bin_edges[self._categories[category]]
            nbins = len(bin_edges) - 1

            # Calculate the success and total for each bin
            indices = np.digitize(values[0], bin_edges, right=False)
            self._totals[self._categories[category]] += np.bincount(indices, minlength=len(bin_edges)+1)[1:-1]
            for ci, (cut, cutname) in enumerate(self._cuts.items()):
                # Sequential cuts (unbinned)
                total = total_events
                success = np.sum(np.all(values[1:ci+2], axis=0))
                self._posteriors[self._categories[category]][f'unbinned_seq_{cut}'] = SpineEfficiency.multiply_posteriors(self._posteriors[self._categories[category]][f'unbinned_seq_{cut}'], binom.pmf(success, total, efficiencies))
                self._selected_counts[group_name][category][f'seq_{cut}'] += int(success)
                # Non-sequential cuts (unbinned)
                success = np.sum(values[ci+1].to_numpy(bool))
                self._posteriors[self._categories[category]][f'unbinned_unseq_{cut}'] = SpineEfficiency.multiply_posteriors(self._posteriors[self._categories[category]][f'unbinned_unseq_{cut}'], binom.pmf(success, total, efficiencies))
                self._selected_counts[group_name][category][f'unseq_{cut}'] += int(success)

                # If the show option is set to 'differential', calculate
                # the binned efficiencies.
                if self._show_option in ('differential', 'final_only'):
                    # Sequential cuts (binned)
                    indices = np.digitize(values[0][np.all(values[1:ci+2], axis=0)], bin_edges, right=False)
                    success = np.bincount(indices, minlength=len(bin_edges)+1)[1:-1]
                    self._successes[self._categories[category]][f'binned_seq_{cut}'] += success
                    _s = self._successes[self._categories[category]][f'binned_seq_{cut}'][:, np.newaxis]  # (nbins, 1)
                    _t = self._totals[self._categories[category]][:, np.newaxis]                          # (nbins, 1)
                    self._posteriors[self._categories[category]][f'binned_seq_{cut}'] = binom.pmf(_s, _t, efficiencies[np.newaxis, :]).astype(np.float32)

                    # Non-sequential cuts (binned)
                    indices = np.digitize(values[0][values[ci+1] == 1], bin_edges, right=False)
                    success = np.bincount(indices, minlength=len(bin_edges)+1)[1:-1]
                    self._successes[self._categories[category]][f'binned_unseq_{cut}'] += success
                    _s = self._successes[self._categories[category]][f'binned_unseq_{cut}'][:, np.newaxis]  # (nbins, 1)
                    _t = self._totals[self._categories[category]][:, np.newaxis]                             # (nbins, 1)
                    self._posteriors[self._categories[category]][f'binned_unseq_{cut}'] = binom.pmf(_s, _t, efficiencies[np.newaxis, :]).astype(np.float32)

    def reduce(self, group, significance=0.6827):
        """
        Reduce the efficiency calculation to the final result (i.e.
        the central value and errors). This calculation sensibly
        treats the unbinned and binned efficiencies in the correct
        way: namely, the calculation preserves the shape of the
        unbinned and binned efficiencies. This is important for the
        usages downstream in the draw function.

        Parameters
        ----------
        group : str
            The key corresponding to the group of efficiencies to
            reduce.
        significance : float, optional
            The significance level to use when calculating the
            efficiency. The default is 0.6827.

        Returns
        -------
        final_posteriors : dict
            A dictionary containing the final efficiency posteriors
            for the requested group.
        cv : dict
            A dictionary containing the central value of the
            efficiency for the requested group.
        msigma : dict
            A dictionary containing the negative error on the
            efficiency for the requested group.
        psigma : dict
            A dictionary containing the positive error on the
            efficiency for the requested group.

        Raises
        ------
        ValueError
            If the group is not in the list of groups configured
            in the analysis block.
        """
        efficiencies = np.linspace(0.0, 1, self._npts, dtype=np.float32)
        
        # If the group is not in the list of groups configured in
        # the analysis block, raise an exception. Otherwise,
        # extract the final posteriors for the group.
        if group not in self._posteriors.keys():
            raise ValueError(f"Group '{group}' not in the list of groups configured in the analysis block!")
        final_posteriors = self._posteriors[group]
        
        # Initialize the dictionaries to store the central value
        # and errors on the efficiency.
        cv = {}
        msigma = {}
        psigma = {}
        
        # The significance level is used to calculate the error on
        # the efficiency. The error is calculated by finding the
        # efficiency value that corresponds to the significance
        # level in the cumulative distribution of the posterior
        # distribution.
        sig = [0.5-significance/2.0, 0.5+significance/2.0]

        for key in final_posteriors.keys():
            if len(final_posteriors[key].shape) == 1:
                final_posteriors[key] /= np.sum(final_posteriors[key])
                cv[key] = efficiencies[int(np.argmax(final_posteriors[key]))]
                cumulative = np.cumsum(final_posteriors[key])
                msigma[key] = cv[key]-efficiencies[int(np.argmax(cumulative > sig[0]))]
                psigma[key] = efficiencies[int(np.argmax(cumulative > sig[1]))]-cv[key]

                # There is a case where posterior distribution
                # peaks at 0.0 or 1.0 due to 100% success or
                # failure. In this case, the error is set to 0.0,
                # which intuitively means that the error bar covers
                # only allowed values of the efficiency.
                if msigma[key] < 0:
                    msigma[key] = 0
                if psigma[key] < 0:
                    psigma[key] = 0

            else:
                final_posteriors[key] /= np.sum(final_posteriors[key], axis=1)[:, np.newaxis]
                cv[key] = efficiencies[np.argmax(final_posteriors[key], axis=1).astype(int)]
                cumulative = np.cumsum(final_posteriors[key], axis=1)
                msigma[key] = cv[key]-efficiencies[np.argmax(cumulative > sig[0], axis=1)]
                psigma[key] = efficiencies[np.argmax(cumulative > sig[1], axis=1)]-cv[key]

                # There is a case where posterior distribution
                # peaks at 0.0 or 1.0 due to 100% success or
                # failure. In this case, the error is set to 0.0,
                # which intuitively means that the error bar covers
                # only allowed values of the efficiency.    
                msigma[key][msigma[key] < 0] = 0
                psigma[key][psigma[key] < 0] = 0

        return final_posteriors, cv, msigma, psigma

    def calculate_purity(self, group, significance=0.6827):
        """
        Calculate purity for a given group.
        
        Purity = (events selected in group) / (total events selected in all groups)
        
        Parameters
        ----------
        group : str
            The name of the group to calculate purity for.
        significance : float, optional
            The significance level for uncertainty calculation. Default is 0.6827 (1 sigma).
        
        Returns
        -------
        tuple : (group, cv, msigma, psigma)
            Dictionaries with cut names as keys and purity values/uncertainties as values.
        """

        if group not in self._selected_counts:
            return {}, {}, {}, {}

        cv = {}
        msigma = {}
        psigma = {}

        sig = [0.5 - significance/2.0, 0.5 + significance/2.0]

        # Precompute totals across all groups for each cut
        total_selected_all = {f'seq_{c}': 0 for c in self._cuts.keys()}
        for g, cats in self._selected_counts.items():
            for cat, cuts_dict in cats.items():
                for k, v in cuts_dict.items():
                    total_selected_all[k] = total_selected_all.get(k, 0) + int(v)

        # Sum selected in this group (group may contain multiple categories)
        group_selected = {f'seq_{c}': 0 for c in self._cuts.keys()}
        if group in self._selected_counts:
            for cat, cuts_dict in self._selected_counts[group].items():
                for k, v in cuts_dict.items():
                    group_selected[k] = group_selected.get(k, 0) + int(v)

        for cut in self._cuts.keys():
            key = f'seq_{cut}'
            n_group = group_selected.get(key, 0)
            n_total = total_selected_all.get(key, 0)

            if n_total == 0:
                cv[f'purity_{key}'] = 0.0
                msigma[f'purity_{key}'] = 0.0
                psigma[f'purity_{key}'] = 0.0
                continue

            purity = n_group / n_total
            cv[f'purity_{key}'] = purity

            # Bayesian uncertainty using Beta distribution
            n_background = n_total - n_group
            lower = beta.ppf(sig[0], n_group + 1, n_background + 1)
            upper = beta.ppf(sig[1], n_group + 1, n_background + 1)

            msigma[f'purity_{key}'] = purity - lower
            psigma[f'purity_{key}'] = upper - purity

        return group, cv, msigma, psigma


    def statistical_efficiency_covariance(
        self,
        x,
        mask_den,
        mask_num,
        bin_edges,
        empty_value=0.0,
    ):
        """Diagonal (independent-bins) statistical covariance of efficiency.

        For each bin b, with Den_b trials and Num_b successes:
            eff_b = Num_b / Den_b
            Var(eff_b) ~= eff_b * (1 - eff_b) / Den_b

        Notes
        -----
        - This is the usual binomial approximation for the uncertainty on an
          estimated efficiency in a bin.
        - Bins are treated as independent => covariance is diagonal.
        """
        mask_den = np.asarray(mask_den, dtype=bool)
        mask_num = np.asarray(mask_num, dtype=bool)
        if mask_den.shape != mask_num.shape:
            raise ValueError("mask_den and mask_num must have same shape")
        if not np.all(~mask_num | mask_den):
            raise ValueError("mask_num must be a subset of mask_den for efficiency.")

        bin_edges = np.asarray(bin_edges, dtype=float)
        B = len(bin_edges) - 1
        if B <= 0:
            return np.zeros((0, 0), dtype=float)

        x = np.asarray(x)
        b_idx = np.digitize(x, bin_edges) - 1
        in_range = (b_idx >= 0) & (b_idx < B)

        den_counts = np.bincount(b_idx[mask_den & in_range], minlength=B).astype(float)
        num_counts = np.bincount(b_idx[mask_num & in_range], minlength=B).astype(float)

        eff = np.full(B, float(empty_value), dtype=float)
        good = den_counts > 0
        eff[good] = num_counts[good] / den_counts[good]

        var = np.zeros(B, dtype=float)
        var[good] = eff[good] * (1.0 - eff[good]) / den_counts[good]

        return np.diag(var)

    def efficiency_covariance_for_systematic(
        self,
        sample,
        systematic,
        x_key,
        mask_den,
        mask_num,
        bin_edges,
        nuniv,
        seed,
        empty_value,
        return_yields=False,
    ):
        """Return an efficiency covariance matrix for a Systematic-like object.

        Handles:
        - Combined systematics (`_components`) by summing component covariances.
        - Handle-less statistical systematics (`_name == 'statistical'`).
        - Universe-based systematics via `efficiency_covariance()`.

        Returns None if the systematic cannot provide efficiency information.
        """
        if systematic is None:
            return None

        # Combined systematic: sum component covariances.
        components = getattr(systematic, '_components', None)
        if components:
            if return_yields:
                raise ValueError("Cannot use return_yields=True for combined systematics; request yields from components instead.")

            cov_sum = None
            for comp in components:
                cov_i = self.efficiency_covariance_for_systematic(
                    sample=sample,
                    systematic=comp,
                    x_key=x_key,
                    mask_den=mask_den,
                    mask_num=mask_num,
                    bin_edges=bin_edges,
                    nuniv=nuniv,
                    seed=seed,
                    empty_value=empty_value,
                    return_yields=return_yields,
                )
                # print(f"Debug: Component systematic '{getattr(comp, '_name', 'unknown')}' covariance:\n{cov_i}")
                if cov_i is None:
                    continue
                cov_sum = cov_i if cov_sum is None else (cov_sum + cov_i)
            return cov_sum

        # Handle-less leaf: interpret only if it is the statistical systematic.
        if getattr(systematic, '_handle', None) is None:
            if str(getattr(systematic, '_name', '')).lower() == 'statistical':
                x_all = sample._data[x_key].to_numpy()
                return self.statistical_efficiency_covariance(
                    x=x_all,
                    mask_den=mask_den,
                    mask_num=mask_num,
                    bin_edges=bin_edges,
                    empty_value=empty_value,
                )
            return None

        # Universe-based leaf.
        return systematic.efficiency_covariance(
            sample=sample,
            x_key=x_key,
            mask_den=mask_den,
            mask_num=mask_num,
            bin_edges=bin_edges,
            nuniv=nuniv,
            seed=seed,
            empty_value=empty_value,
            return_yields=return_yields,
        )

    def compute_efficiency_with_systematics(
        self,
        sample,
        x_key,
        mask_den,
        mask_num,
        bin_edges,
        draw_error,
        systematic=None,
        nuniv=None,
        seed=None,
        empty_value=0.0,
        return_yields=False,
    ):
        """Return efficiency covariance for the requested systematic key/name."""
        if nuniv is None:
            nuniv = self._nuniv
        if seed is None:
            seed = self._seed

        if not hasattr(sample, '_systematics'):
            return None
        # Resolve which systematic to use.
        # - If `systematic` is provided, use it directly.
        # - Otherwise, resolve from `draw_error` key/name in the sample.
        syst = systematic
        if syst is None:
            # Try to find the systematic by key first (typical for recipe names like 'total_stat').
            syst = sample._systematics.get(draw_error)

            # If that fails, fall back to matching by the underlying systematic name
            # (e.g. draw_error='statistical' but the stored key is 'NuMIFull_statistical').
            if syst is None:
                want = str(draw_error).lower()
                matches = [
                    s for s in sample._systematics.values()
                    if str(getattr(s, '_name', '')).lower() == want
                ]
                if len(matches) == 1:
                    syst = matches[0]

        # print(f"Debug: For sample '{sample._name}', draw_error '{draw_error}', found systematic: {syst}", flush=True)
        if syst is None:
            print(f"Debug: No systematic matching '{draw_error}' for sample '{sample._name}'; skipping.", flush=True)
            return None

        mask_den = np.asarray(mask_den, dtype=bool)
        mask_num = np.asarray(mask_num, dtype=bool)

        cov = self.efficiency_covariance_for_systematic(
            sample=sample,
            systematic=syst,
            x_key=x_key,
            mask_den=mask_den,
            mask_num=mask_num,
            bin_edges=bin_edges,
            nuniv=nuniv,
            seed=seed,
            empty_value=empty_value,
            return_yields=return_yields,
        )

        # print(f"Debug: Efficiency covariance for sample '{sample._name}', systematic '{draw_error}':\n{cov}")
        return cov


    def compute_syst_sigma_eff(self, group, cut, sequential=True):
        """Return per-bin systematic sigma on efficiency (absolute units, 0–1).

        Important: efficiency is a ratio. When combining multiple samples,
        combine per-universe (Num, Den) yields across samples first, then build
        per-universe efficiencies and compute the covariance.

        For recipe/combined systematics, recursively flatten to leaf systematics
        and de-duplicate leaves so the same branch isn't counted twice.
        """
        variable = self._variable
        x_key = variable._key
        bin_edges = np.asarray(variable._bin_edges[group], dtype=float)
        B = len(bin_edges) - 1
        if B <= 0:
            return None

        group_categories = [
            cat_id for cat_id, g_label in self._categories.items() if g_label == group
        ]
        if not group_categories:
            return None

        def _resolve_systematic(sample):
            if not hasattr(sample, '_systematics'):
                return None

            syst = sample._systematics.get(self._draw_error)
            if syst is not None:
                return syst

            want = str(self._draw_error).lower()
            matches = [
                s for s in sample._systematics.values()
                if str(getattr(s, '_name', '')).lower() == want
            ]
            if len(matches) == 1:
                return matches[0]
            return None

        def _iter_leaves(syst):
            components = getattr(syst, '_components', None)
            if components:
                for comp in components:
                    yield from _iter_leaves(comp)
                return
            yield syst

        def _leaf_key(leaf):
            name = str(getattr(leaf, '_name', ''))
            handle = getattr(leaf, '_handle', None)
            if handle is None:
                return name
            hname = getattr(handle, 'name', None)
            if hname is None:
                hname = getattr(handle, 'fName', None)
            if hname is None:
                hname = f'handle_{id(handle)}'
            return f"{name}::{hname}"

        def _cov_from_num_den(Num_tot, Den_tot, E_cv, empty_value=0.0):
            U, Bb = Num_tot.shape
            if Bb != B:
                raise ValueError("Bin-count mismatch while building efficiency covariance")
            if U < 1:
                return np.zeros((B, B), dtype=float)

            E_u = np.full((U, B), float(empty_value), dtype=float)
            good = Den_tot > 0.0
            E_u[good] = Num_tot[good] / Den_tot[good]

            diff = E_u - E_cv[np.newaxis, :]
            Cov = (diff.T @ diff) / float(U)
            return Cov

        NumDen_by_leaf = {}  # leaf_key -> (Num_tot, Den_tot)

        stat_den_counts = np.zeros(B, dtype=float)
        stat_num_counts = np.zeros(B, dtype=float)

        cv_den_counts = np.zeros(B, dtype=float)
        cv_num_counts = np.zeros(B, dtype=float)

        contributed = False

        for sample in self._samples:
            n_events = len(sample._data)

            var_mask_expr = getattr(variable, "mask", None)
            if var_mask_expr:
                mask_var = sample._data.eval(var_mask_expr).to_numpy(dtype=bool)
            else:
                mask_var = np.ones(n_events, dtype=bool)

            cat_series = sample._data[sample._category_branch]
            mask_group = cat_series.isin(group_categories).to_numpy(dtype=bool)

            base_mask = mask_var & mask_group
            mask_den = base_mask

            if cut:
                cuts_order = list(self._cuts.keys())
                if cut not in cuts_order:
                    raise ValueError(f"Cut '{cut}' not found in configured cuts: {cuts_order}")

                if sequential:
                    cut_index = cuts_order.index(cut)
                    required_cuts = cuts_order[: cut_index + 1]
                else:
                    required_cuts = [cut]

                missing = [c for c in required_cuts if c not in sample._data.columns]
                if missing:
                    continue

                mask_all = np.ones(n_events, dtype=bool)
                for c in required_cuts:
                    mask_all &= sample._data[c].to_numpy(dtype=bool)
                mask_num = base_mask & mask_all
            else:
                mask_num = base_mask

            if not mask_den.any():
                continue

            syst = _resolve_systematic(sample)
            if syst is None:
                continue

            # CV counts for this sample/group/cut (used for CV-centered covariance)
            x_all = sample._data[x_key].to_numpy()
            b_idx = np.digitize(x_all, bin_edges) - 1
            in_range = (b_idx >= 0) & (b_idx < B)
            cv_den_counts += np.bincount(b_idx[mask_den & in_range], minlength=B).astype(float)
            cv_num_counts += np.bincount(b_idx[mask_num & in_range], minlength=B).astype(float)

            # De-dup leaves within this sample to avoid counting the same branch twice
            # when recipes contain repeated entries.
            seen = set()

            for leaf in _iter_leaves(syst):
                leaf_name = str(getattr(leaf, '_name', '')).lower()
                leaf_handle = getattr(leaf, '_handle', None)

                sig = _leaf_key(leaf)
                if sig in seen:
                    continue
                seen.add(sig)

                if leaf_handle is None:
                    if leaf_name == 'statistical':
                        x_all = sample._data[x_key].to_numpy()
                        b_idx = np.digitize(x_all, bin_edges) - 1
                        in_range = (b_idx >= 0) & (b_idx < B)
                        stat_den_counts += np.bincount(b_idx[mask_den & in_range], minlength=B).astype(float)
                        stat_num_counts += np.bincount(b_idx[mask_num & in_range], minlength=B).astype(float)
                        contributed = True
                    continue

                out = self.compute_efficiency_with_systematics(
                    sample=sample,
                    x_key=x_key,
                    mask_den=mask_den,
                    mask_num=mask_num,
                    bin_edges=bin_edges,
                    draw_error=self._draw_error,
                    systematic=leaf,
                    nuniv=self._nuniv,
                    seed=self._seed,
                    empty_value=0.0,
                    return_yields=True,
                )
                if out is None:
                    continue

                _cov_i, Num, Den = out

                key = _leaf_key(leaf)
                if key not in NumDen_by_leaf:
                    NumDen_by_leaf[key] = (Num.copy(), Den.copy())
                else:
                    Num_tot, Den_tot = NumDen_by_leaf[key]
                    if Num_tot.shape != Num.shape:
                        raise ValueError("Universe-count mismatch while combining leaf yields")
                    NumDen_by_leaf[key] = (Num_tot + Num, Den_tot + Den)

                contributed = True

        if not contributed:
            return None

        cov_total = None

        E_cv = np.full(B, 0.0, dtype=float)
        good_cv = cv_den_counts > 0.0
        E_cv[good_cv] = cv_num_counts[good_cv] / cv_den_counts[good_cv]

        for _, (Num_tot, Den_tot) in NumDen_by_leaf.items():
            cov_i = _cov_from_num_den(Num_tot, Den_tot, E_cv, empty_value=0.0)
            cov_total = cov_i if cov_total is None else (cov_total + cov_i)

        if stat_den_counts.sum() > 0:
            eff = np.zeros(B, dtype=float)
            good = stat_den_counts > 0
            eff[good] = stat_num_counts[good] / stat_den_counts[good]
            var = np.zeros(B, dtype=float)
            var[good] = eff[good] * (1.0 - eff[good]) / stat_den_counts[good]
            cov_stat = np.diag(var)
            cov_total = cov_stat if cov_total is None else (cov_total + cov_stat)

        if cov_total is None:
            return None

        # print(f"Debug: CV denominator histogram for group '{group}', cut '{cut}': {cv_den_counts}")

        # Track denominator covariance for all contributing leaves.
        # denom_cov_total = None
        # for leaf_key, (_Num_leaf, Den_leaf) in NumDen_by_leaf.items():
        #     U_den = Den_leaf.shape[0]
        #     denom_diff = Den_leaf - cv_den_counts[np.newaxis, :]
        #     denom_cov_leaf = (denom_diff.T @ denom_diff) / float(U_den)
        #     print(
        #         f"Debug: Denominator covariance for leaf '{leaf_key}' (group '{group}', cut '{cut}'):\n{denom_cov_leaf}"
        #     )
        #     denom_cov_total = denom_cov_leaf if denom_cov_total is None else (denom_cov_total + denom_cov_leaf)

        # if denom_cov_total is not None:
        #     print(
        #         f"Debug: Total denominator covariance across leaves for group '{group}', cut '{cut}':\n{denom_cov_total}"
        #     )

        # print(f"Debug: Total efficiency covariance for group '{group}', cut '{cut}':\n{cov_total}", flush=True)


        sigma_eff = np.sqrt(np.clip(np.diag(cov_total), 0.0, np.inf))
        # print(
        #     f"Debug: Per-bin systematic sigma on efficiency for group '{group}', cut '{cut}': {sigma_eff}",
        #     flush=True,
        # )
        return sigma_eff
