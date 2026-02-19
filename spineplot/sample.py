import numpy as np
import pandas as pd
import re
import uproot

from systematic import Systematic

class Sample:
    """
    A class designed to encapsulate the data for a single sample and
    all associated functionality and metadata. Data is loaded in
    batches using uproot.iterate to reduce peak memory usage; the
    full dataset is never held in memory as a single DataFrame.

    Attributes
    ----------
    _name : str
        The name of the sample.
    _exposure_type : str
        The exposure type for the sample. This can be either 'pot' or
        'livetime'.
    _file_handle : uproot.reading.ReadOnlyDirectory
        The file handle for the input ROOT file.
    _exposure_pot : float
        The exposure of the sample in POT.
    _exposure_livetime : float
        The exposure of the sample in livetime.
    _category_branch : str
        The name of the branch in the TTree containing the category
        labels.
    _columns : set
        The set of column names available in the sample (including
        any precomputed columns).
    _trees : list
        The list of TTree names in the ROOT file to load for the sample.
    _precompute : dict
        A dictionary of new branches to compute from existing branches.
    _presel : str
        A pre-selection expression to apply to the sample data.
    _override_category : int
        The category value to override the category branch with.
    _fillna : int or float
        The value used to fill NaN entries in the category branch.
    _weight : float
        The exposure-scaling weight for the sample. Set by set_weight().
    _systematics : dict
        A dictionary of Systematic objects for the Sample.
    _print_sys : bool
        A boolean flag that toggles the printing of integrated
        systematic uncertainties for the sample.
    _presel_mask : np.ndarray
        A boolean array aligned with all events in the file (before
        presel) indicating which events pass the pre-selection. Stored
        for systematic weight alignment.
    """
    def __init__(self, name, rf, category_branch, key, exposure_type, trees,
                 fillna=None, systematics=None, override_exposure=None, precompute=None,
                 presel=None, override_category=None, print_sys=False) -> None:
        """
        Initializes the Sample object with the given name and key.

        Parameters
        ----------
        name : str
            The name of the sample.
        rf : uproot.reading.ReadOnlyDirectory
            The file handle for the input ROOT file.
        category_branch : str
            The name of the branch in the TTree containing the category
            labels. This categorical information is referenced in the
            configuration file to designate specific components of the
            sample and apply different styles to them.
        key : str
            The key/name of the TDirectory in the ROOT file input
            containing the sample data.
        exposure_type : str
            The exposure type for the sample. This can be either 'pot'
            or 'livetime'. This is used for matching the exposure of
            the sample to the target sample.
        trees : list
            The list of TTree names in the ROOT file to load for the
            sample.
        override_exposure : float, optional
            The exposure value to override the exposure in the ROOT file
            with. The default is None.
        precompute : dict, optional
            A dictionary of new branches to compute from the existing
            branches in the sample. The keys are the names of the new
            branches and the values are the expressions to compute the
            new branches. The default is None.
        presel : str, optional
            A pre-selection string to apply to the sample data. The
            default is None.
        override_category : int
            The category to override the category branch with if it is
            configured. Else, the category branch is left as is.

        Returns
        -------
        None.
        """
        self._name = name
        self._exposure_type = exposure_type
        self._file_handle = rf[f'events/{key}']
        self._exposure_pot = self._file_handle['POT'].to_numpy()[0][0]
        self._exposure_livetime = self._file_handle['Livetime'].to_numpy()[0][0]
        self._category_branch = category_branch
        self._print_sys = print_sys
        self._trees = trees
        self._precompute = precompute
        self._presel = presel
        self._override_category = override_category
        self._fillna = fillna
        self._weight = 1.0

        if override_exposure is not None:
            self.override_exposure(override_exposure, exposure_type)

        # Collect column names from tree keys without loading all data.
        self._columns = set()
        for tree in trees:
            self._columns.update(self._file_handle[tree].keys())
        if precompute is not None:
            self._columns.update(precompute.keys())

        if self._category_branch not in self._columns:
            raise ValueError(f'Category branch `{self._category_branch}` not found in sample `{self._name}`.')

        # Compute and store the presel mask by iterating in batches.
        # The mask is aligned with all events in the file (before
        # presel) and is needed for systematic weight alignment.
        presel_masks = []
        nan_count = 0
        for chunk in self._iter_raw_batches():
            if presel is not None:
                mask = chunk.eval(presel).to_numpy(dtype=bool)
            else:
                mask = np.ones(len(chunk), dtype=bool)
            presel_masks.append(mask)

            # Count NaN categories in the pre-selection subset.
            presel_chunk = chunk[mask]
            nan_count += int(presel_chunk[self._category_branch].isna().sum())

        self._presel_mask = (
            np.concatenate(presel_masks) if presel_masks else np.array([], dtype=bool)
        )

        # Warn about NaN categories found after presel.
        if nan_count > 0:
            if fillna is not None:
                print(
                    f'Filling NaN category in Sample `{self._name}`'
                    f' with `{fillna}`. You asked for this behavior!'
                )
            else:
                print(
                    f'Found NaN category in Sample `{self._name}` with'
                    f' {nan_count} occurrence(s). Masking NaNs...\n'
                    f' [!!!] Please note that this automatic behavior'
                    f' may lead to issues with systematics!'
                )

        # Initialize the systematics dictionary for the sample. Note:
        # the sample will always have a statistical uncertainty.
        self._systematics = dict()

        # Load the systematic uncertainties for the sample. If no
        # systematics are provided, the sample is assumed to have no
        # systematic uncertainties.
        if systematics is not None:    
            for sys in systematics:
                systs = [k for k in self._file_handle[sys].keys() if k not in ['Run', 'Subrun', 'Evt']]
                self._systematics.update({syst: Systematic(syst, self._file_handle[sys][syst]) for syst in systs})
        
        # Add statistical uncertainty. This can always be added to the
        # sample, because it is not dependent on some external source
        # of weights.
        self._systematics.update({f'{self._name}_statistical': Systematic('statistical', None)})

    def _iter_raw_batches(self, step_size='100 MB'):
        """
        Internal generator that yields raw DataFrame chunks from all
        configured trees, applying precompute expressions and the
        category override to each chunk.

        Parameters
        ----------
        step_size : str or int, optional
            The size of each batch passed to uproot.iterate.  Accepts
            the same values as uproot's step_size parameter (e.g.
            '100 MB', 10000). The default is '100 MB'.

        Yields
        ------
        chunk : pd.DataFrame
            A batch of events with precompute columns added.
        """
        for tree in self._trees:
            for chunk in self._file_handle[tree].iterate(
                library='pd', step_size=step_size
            ):
                if self._precompute is not None:
                    for k, v in self._precompute.items():
                        chunk[k] = chunk.eval(v)
                if self._override_category is not None:
                    chunk[self._category_branch] = self._override_category
                yield chunk

    def iterate_batches(self, variables, with_mask=None, step_size='100 MB'):
        """
        Iterate over the sample data in batches, yielding one
        (data_dict, weights_dict) pair per batch and per category.
        The presel, NaN-category handling, and any additional mask are
        applied to each chunk before yielding.

        This is the primary method for memory-efficient data access.
        Artists should prefer this method over get_data() when they
        accumulate results incrementally (e.g. histograms).

        Parameters
        ----------
        variables : list[str]
            The column names to include in the yielded data.
        with_mask : str, optional
            A pandas eval expression used as an additional row filter.
            The default is None.
        step_size : str or int, optional
            The size of each batch. The default is '100 MB'.

        Yields
        ------
        data : dict
            Maps int category label -> list of pd.Series (one per
            variable in *variables*).
        weights : dict
            Maps int category label -> pd.Series of per-event weights.
        """
        for chunk in self._iter_raw_batches(step_size=step_size):
            # Apply presel filter.
            if self._presel is not None:
                presel_mask = chunk.eval(self._presel).to_numpy(dtype=bool)
                chunk = chunk[presel_mask]

            if len(chunk) == 0:
                continue

            # Handle NaN categories per batch.
            nan_mask = chunk[self._category_branch].isna()
            if nan_mask.any():
                if self._fillna is not None:
                    chunk = chunk.copy()
                    chunk.loc[nan_mask, self._category_branch] = self._fillna
                else:
                    chunk = chunk[~nan_mask]
                if len(chunk) == 0:
                    continue

            # Apply optional additional mask.
            if with_mask is not None:
                extra_mask = chunk.eval(with_mask).to_numpy(dtype=bool)
            else:
                extra_mask = np.ones(len(chunk), dtype=bool)

            data = {}
            weights = {}
            for category in np.unique(chunk[self._category_branch]):
                if np.isnan(category):
                    continue
                cat_mask = (chunk[self._category_branch] == category) & extra_mask
                cat_data = chunk[cat_mask]
                data[int(category)] = [cat_data[v] for v in variables]
                weights[int(category)] = pd.Series(
                    np.full(len(cat_data), self._weight),
                    index=cat_data.index,
                )

            yield data, weights

    def load_column(self, column):
        """
        Load a single column from the sample, applying presel and NaN
        category handling, and return it as a numpy array.  This
        method is intended for use by Systematic.process() which
        requires per-event values aligned with the presel-filtered
        dataset.

        Parameters
        ----------
        column : str
            The name of the column to load.

        Returns
        -------
        np.ndarray
            The column values after applying presel and NaN filtering,
            in the same order as events are stored in the file.
        """
        chunks = []
        for chunk in self._iter_raw_batches():
            if self._presel is not None:
                mask = chunk.eval(self._presel).to_numpy(dtype=bool)
                chunk = chunk[mask]
            if len(chunk) == 0:
                continue
            nan_mask = chunk[self._category_branch].isna()
            if nan_mask.any():
                if self._fillna is not None:
                    chunk = chunk.copy()
                    chunk.loc[nan_mask, self._category_branch] = self._fillna
                else:
                    chunk = chunk[~nan_mask]
            if len(chunk) > 0:
                chunks.append(chunk[[column]])
        if chunks:
            return pd.concat(chunks)[column].to_numpy()
        return np.array([])

    def override_exposure(self, exposure, exposure_type='pot') -> None:
        """
        Overrides the exposure for the sample. This is useful for
        setting the exposure for samples for which the exposure is not
        valid. The exposure type can be either 'pot' or 'livetime'. It
        is not recommended to use this method unless the exposure is
        known to be incorrect.

        Parameters
        ----------
        exposure : float
            The exposure to set for the sample.
        exposure_type : str
            The type of exposure to set. This can be either 'pot' or
            'livetime'. The default is 'pot'.

        Returns
        -------
        None.
        """
        if exposure_type == 'pot':
            self._exposure_pot = exposure
        else:
            self._exposure_livetime = exposure

    def set_weight(self, target=None) -> None:
        """
        Sets the weight for the sample to the target value.

        Parameters
        ----------
        target : Sample
            The Sample object to use as the exposure normalization
            target. This is used to scale the weight of this sample to
            the target sample. If None, the weight is set to 1.
        
        Returns
        -------
        None.
        """
        if target is None:
            self._weight = 1.0
        elif self._exposure_type == 'pot':
            self._weight = target._exposure_pot / self._exposure_pot
        else:
            self._weight = target._exposure_livetime / self._exposure_livetime

        print(f"Setting weight for {self._name} to {self._weight:.2e}")
        for syst in self._systematics.values():
            syst.set_weight(self._weight)

    def get_data(self, variables, with_mask=None) -> dict:
        """
        Returns the data for the given variable(s) in the sample. The
        data is accumulated from all batches and returned as a
        dictionary with the category as the key and the data for the
        requested variable as the value.

        Note: for memory-critical code paths (e.g. large spectra) it
        is preferable to use iterate_batches() directly to avoid
        accumulating the full dataset in memory.

        Parameters
        ----------
        variables : list[str]
            The names of the variables to retrieve.
        with_mask : str, optional
            A mask formula to apply to the variable. The default is None.

        Returns
        -------
        data : dict
            The data for the requested variable in the sample. The data
            is stored as a dictionary with the category as the key and
            the data (a list of pandas Series, one per variable) as the
            value.
        weights : dict
            The weights for the requested variable in the sample. The
            weights are stored as a dictionary with the category as the
            key and the weights (a pandas Series) as the value.
        """
        data = {}
        weights = {}
        for batch_data, batch_weights in self.iterate_batches(variables, with_mask):
            for category, values in batch_data.items():
                if category not in data:
                    data[category] = [pd.Series(dtype=float) for _ in variables]
                    weights[category] = pd.Series(dtype=float)
                for i, v in enumerate(values):
                    data[category][i] = pd.concat([data[category][i], v])
                weights[category] = pd.concat([weights[category], batch_weights[category]])
        return data, weights

    def process_systematics(self, recipes) -> None:
        """
        Processes the systematic uncertainties for the sample.

        Parameters
        ----------
        recipes : list
            A list of dictionaries containing the recipes for
            combinations of systematic uncertainties.

        Returns
        -------
        None.
        """
        for syst in self._systematics.values():
            syst.process(self, self._presel_mask)
                
        # Each recipe has a name, which is used to identify the
        # combination of systematic uncertainties, and a pattern,
        # which is used to build a list of Systematic objects to
        # combine.
        for recipe in recipes:
            # Exclude the "nsigma" branches
            exclude = ['_nsigma', '_sigma']
            exclude_pat = '|'.join(re.escape(x) for x in exclude)

            # Compile the regex pattern for matching the systematic
            pattern = recipe['pattern']
            regxp = re.compile(rf'^(?!.*(?:{exclude_pat})).*{pattern}.*$')
            systematics = [syst for k, syst in self._systematics.items() if regxp.match(k)]

            # If there are no systematics to combine, skip the recipe.
            if len(systematics) == 0:
                continue

            # Combine the systematics and add the new Systematic object
            syst = Systematic.combine(systematics, recipe['name'], recipe.get('label', None))
            self._systematics[syst._name] = syst
        
        # Print the systematics for the sample (if requested).
        if self._print_sys:
            for sysname, syst in self._systematics.items():
                print(syst)

    def register_variable(self, variable, categories) -> None:
        """
        Registers a new Variable object with the Sample object. This
        allows the Sample object to call the Variable object's method
        to check the Variable's validity in the Sample. Additionally,
        this allows the Sample object to create or populate a
        Systematic object with a covariance matrix for the Variable.

        Parameters
        ----------
        variable : Variable
            The Variable object to register with the Sample object.
        categories : dict
            A dictionary containing the categories for the analysis.
            The key is the category enumeration and the value is the
            name of the group that the enumerated category belongs to.
        
        Returns
        -------
        None.
        """
        variable.check_data(categories, self._name, self)
        for syst in self._systematics.values():
            syst.register_variable(variable)

    def __str__(self) -> str:
        """
        Returns a string representation of the Sample object.
        
        Parameters
        ----------
        None.

        Returns
        -------
        res : str
            A string representation of the Sample object.
        """
        res = f'{"Sample:":<15}{self._name}'
        res += f'\n{"POT:":<15}{self._exposure_pot:.2e}'
        res += f'\n{"Livetime:":<15}{self._exposure_livetime:.2e}'
        return res

    def evaluate_formula(self, formula) -> pd.Series:
        """
        Evaluates the given formula across all batches of the sample
        data (post presel) and returns the concatenated result.

        Parameters
        ----------
        formula : str
            The formula to evaluate.

        Returns
        -------
        result : pd.Series
        """
        chunks = []
        for chunk in self._iter_raw_batches():
            if self._presel is not None:
                mask = chunk.eval(self._presel).to_numpy(dtype=bool)
                chunk = chunk[mask]
            if len(chunk) > 0:
                chunks.append(chunk.eval(formula))
        return pd.concat(chunks) if chunks else pd.Series(dtype=float)