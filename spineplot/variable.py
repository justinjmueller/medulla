import numpy as np
import pandas as pd

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
    _validity_check : dict
        A dictionary containing the validity check for the variable
        in each sample. The key is the sample name and the value is
        a boolean indicating whether the variable is present in the
        sample. This is intended to be checked before using the
        variable in the analysis.
    """
    def __init__(self, name, key, range, nbins,
                 binning_scheme='equal_width', xlabel=None,
                 mask=None) -> None:
        """
        Initializes the Variable object with the given kwargs.

        Parameters
        ----------
        name : str
            The name of the variable.
        key : str
            The key/name of the branch in the input TTree containing the
            variable data.
        range : tuple
            The range of the variable.
        nbins : int
            The number of bins for the variable.
        binning_scheme : str
            The binning scheme for the variable. This can be either
            'equal_width' or 'equal_population'. The default is
            'equal_width,' which creates bins of equal width
            irrespective of the number of entries in each bin.
        xlabel : str
            The x-axis label for the variable.
        mask : string, optional
            A mask formula to apply to the variable. The default is None.

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
        self._validity_check = {}
        self._bin_edges = {}
        self._bin_centers = {}
        self._bin_widths = {}

    def check_data(self, categories, sample_name, sample):
        """
        Provides some functionality to check the data for the variable.
        Specifically, this function checks that each sample does indeed
        have the key for the variable. Additionally, it is also useful
        to parse the categories and provide a separate binning for each
        category group in the case of equal-population binning.

        If the variable is not found in a sample, a flag is set for
        the variable instance to indicate that the variable is missing
        and should not be used. If the variable is used later in the 
        analysis, the missing variable flag will be checked and an
        exception will be raised.

        Parameters
        ----------
        categories : dict
            A dictionary containing the categories for the analysis.
            The key is the category enumeration and the value is the
            name of the group that the enumerated category belongs to.
        sample_name : str
            The name of the Sample to be checked for validity.
        sample : Sample
            The Sample to be checked for validity/availability of the
            variable.

        Returns
        -------
        None.
        """
        self._validity_check[sample_name] = (self._key in sample._data.keys())

        if all(self._validity_check.values()):
            groups = {v: [] for v in categories.values()}
            for k, v in sample.get_data([self._key], self._mask)[0].items():
                if k in categories.keys():
                    groups[categories[k]].append(v[0])

            for g, v in groups.items():
                if v and self._binning_scheme == 'equal_population':
                    all_entries = pd.concat(v)
                    range_mask = ((all_entries >= self._range[0]) & (all_entries <= self._range[1]))
                    self._bin_edges[g] = np.percentile(all_entries[range_mask], np.linspace(0, 100, self._nbins+1))
                    self._bin_centers[g] = 0.5*(self._bin_edges[g][1:] + self._bin_edges[g][:-1])
                    self._bin_widths[g] = self._bin_edges[g][1:] - self._bin_edges[g][:-1]
                else:
                    self._bin_edges[g] = np.linspace(self._range[0], self._range[1], self._nbins+1)
                    self._bin_centers[g] = 0.5*(self._bin_edges[g][1:] + self._bin_edges[g][:-1])
                    self._bin_widths[g] = self._bin_edges[g][1:] - self._bin_edges[g][:-1]

    @property
    def mask(self):
        """
        Getter method for the mask attribute.

        Parameters
        ----------
        None.

        Returns
        -------
        string
            The mask formula for the variable.
        """
        return self._mask