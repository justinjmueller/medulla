import numpy as np
import pandas as pd
from typing import Optional

class ConfigException(Exception):
    pass

class Variable:
    """
    A class designed to encapsulate the configuration of a single
    variable for the analysis.

    Attributes
    ----------
    _name : str
        The name of the variable.
    _key : str
        The key/name of the branch in the input TTree containing the
        variable data.
    _range : tuple
        The range of the variable.
    _nbins : int
        The number of bins for the variable.
    _xlabel : str
        The x-axis label for the variable.
    _mask : string
        A mask formula to apply to the variable.
    _bin_edges : numpy.ndarray
        The bin edges for the variable.
    _bin_centers : numpy.ndarray
        The bin centers for the variable.
    _bin_widths : numpy.ndarray
        The bin widths for the variable.
    _valid : dict
        A dictionary containing the validity check for the variable
        in each sample. The key is the sample name and the value is
        a boolean indicating whether the variable is present in the
        sample. This is intended to be checked before using the
        variable in the analysis.
    """
    def __init__(
        self,
        name : str,
        key : str,
        range : Optional[tuple] = None,
        nbins : Optional[int] = None,
        binning_scheme : str = 'equal_width',
        xlabel : Optional[str] = None,
        mask : Optional[str] = None,
        custom_bins : Optional[np.ndarray] = None
    ) -> None:
        """
        Initializes the Variable object with the given kwargs.

        Parameters
        ----------
        name : str
            The name of the variable.
        key : str
            The key/name of the branch in the input TTree containing
            the variable data.
        range : tuple
            The range of the variable.
        nbins : int
            The number of bins for the variable.
        binning_scheme : str
            The binning scheme for the variable; that is, how the bin
            edges are determined. Several options are available:            
            - 'equal_width': Creates a fixed number of bins with equal
            width over the specified range.
            - 'equal_population': Creates bins with an equal number of
            entries (as close as possible) in each bin.
            - 'custom_bins': Uses user-defined bin edges specified in
            the custom_bins parameter.
        xlabel : str
            The x-axis label for the variable.
        mask : string, optional
            A mask formula to apply to the variable. The default is
            None.
        custom_bins : np.ndarray, optional
            The customized bin edges (required if binning_scheme is
            'custom_bins').
        #TODO: consider refactoring the binning scheme handling.

        Returns
        -------
        None
        """
        self._name = name
        self._key = key
        self._range = range
        self._nbins = nbins
        self._binning_scheme = binning_scheme
        self._xlabel = xlabel
        self._mask = mask
        self._custom_bins = custom_bins

        # Initialize the binning structures as empty dictionaries.
        # These will be populated in the finalize() method.
        self._bin_edges = dict()
        self._bin_centers = dict()
        self._bin_widths = dict()

        # Each binning scheme requires a different set of parameters to
        # be defined. Check the validity of the provided parameters
        # based on the selected binning scheme. Raise a ConfigException
        # if any required parameters are missing or invalid.

        # First, check that the binning scheme option is a supported
        # one.
        if self._binning_scheme not in ['equal_width',
                                        'equal_population',
                                        'custom_bins']:
            raise ConfigException(
                "Invalid binning scheme. Available options are: "
                "'equal_width', 'equal_population', 'custom_bins'."
                f" Variable: '{self._name}'."
            )

        # If the user requests equal-population binning, ensure that
        # the number of bins is at least 1 and the range is valid.
        if self._binning_scheme == 'equal_population':
            if self._nbins < 1:
                raise ConfigException(
                    "Number of bins must be at least 1 when using"
                    " 'equal_population' binning scheme. Variable: "
                    f"'{self._name}'."
                )
            if self._range[0] >= self._range[1]:
                raise ConfigException(
                    "Invalid range specified for variable when using"
                    " 'equal_population' binning scheme. The lower "
                    "bound must be less than the upper bound. "
                    f"Variable: '{self._name}'."
                )

        # If the user requests custom bins, ensure that the custom
        # bins are actually provided, contain at least two edges, and
        # are in strictly increasing order.
        elif self._binning_scheme == 'custom_bins':
            self._range = self._custom_bins[0], self._custom_bins[-1]
            if self._custom_bins is None:
                raise ConfigException(
                    "Custom bins must be provided when using"
                    " 'custom_bins' binning scheme. Variable: "
                    f"'{self._name}'."
                )
            if len(self._custom_bins) < 2:
                raise ConfigException(
                    "At least two bin edges must be provided for"
                    " 'custom_bins' binning scheme. Variable: "
                    f"'{self._name}'."
                )
            if not np.all(np.diff(self._custom_bins) > 0):
                raise ConfigException(
                    "Custom bin edges must be in strictly increasing"
                    " order. Variable: '{self._name}'."
                )
        
        # If the user requests equal-width binning, ensure that the
        # number of bins is at least 1 and the range is valid.
        elif self._binning_scheme == 'equal_width':
            if self._nbins < 1:
                raise ConfigException(
                    "Number of bins must be at least 1 when using"
                    " 'equal_width' binning scheme. Variable: "
                    f"'{self._name}'."
                )
            if self._range[0] >= self._range[1]:
                raise ConfigException(
                    "Invalid range specified for variable when using"
                    " 'equal_width' binning scheme. The lower bound"
                    " must be less than the upper bound. Variable: "
                    f"'{self._name}'."
                )

    def finalize(
        self,
        samples : list[Sample],
        categories : dict
    ) -> None:
        """
        Finalizes the Variable creation by performing the necessary
        checks on the provided samples to ensure that the variable
        definition is valid for all samples.

        This function also prepares the binning structures for the
        variable based on the user-specified binning scheme. This is
        done once per group of categories in the analysis, which can be
        relevant for the equal-population binning scheme.

        If the variable is not found in a sample, a flag is set for
        the variable instance to indicate that the variable is missing
        and should not be used. This missing variable flag is intended
        to be checked whenever an artist or other analysis component
        attempts to use the variable.

        Parameters
        ----------
        samples : list[Sample]
            A list of Sample objects to check for the variable's
            presence.
        categories : dict
            A dictionary mapping category names to their respective
            group names.

        Returns
        -------
        None.
        """
        # Perform the validity check for each sample. For now, this is
        # a simple check to see if the variable key is present in the
        # Sample object.
        self._valid = {s.name : self._key in s._data for s in samples}

        # If the variable is valid in all samples, proceed to prepare
        # the binning structures.
        if not all(self._valid.values()):
            return

        # Prepare the binning structures for each group of categories.
        # This is relevant for the equal-population binning scheme,
        # where the bin edges may differ for each group.
        for s in samples:
            # Prepare the data for each group of categories, respecting
            # the mask if provided.
            groups = {v: [] for v in categories.values()}
            for k, v in s.get_data([self._key], self._mask)[0].items():
                if k in categories.keys():
                    groups[categories[k]].append(v[0])

            # Now compute the bin edges for each group based on the
            # specified binning scheme.
            data = pd.concat(v)
            mask = (
                (data >= self._range[0]) &
                (data <= self._range[1])
            )
            for g, v in groups.items():
                # Equal-population binning is intended to create bins
                # of approximately equal population based on the data
                # distribution.
                if v and self._binning_scheme == 'equal_population':
                    self._bin_edges[g] = np.percentile(
                        data[mask],
                        np.linspace(0, 100, self._nbins+1)
                    )

                # Custom bins use user-defined bin edges and provide
                # the finest level of control over the binning.
                elif self._binning_scheme == 'custom_bins':
                    self._bin_edges[g] = np.array(self._custom_bins)
                    self._nbins = len(self._custom_bins) - 1
                    self._range = self._bin_edges[g][[0, -1]]
                    
                # Default to equal-width binning if no valid scheme is
                # specified. Creates bins of equal width (uniformly)
                # over the specified range.
                else:
                    self._bin_edges[g] = np.linspace(
                        self._range[0],
                        self._range[1],
                        self._nbins+1
                    )

                # Compute the bin centers and widths from the bin
                # edges. Fortunately, this is shared across all binning
                # schemes once the edges are determined.
                stack = np.stack([
                    self._bin_edges[g][1:],
                    self._bin_edges[g][:-1]
                ], axis=1)
                self._bin_centers[g] = 0.5*np.sum(
                    stack,
                    axis=1
                )
                self._bin_widths[g] = np.subtract(
                    stack[:, 0],
                    stack[:, 1]
                )
    @property
    def mask(self):
        """
        Getter method for the mask attribute.

        Parameters
        ----------
        None.

        Returns
        -------
        str
            The mask formula for the variable.
        """
        return self._mask