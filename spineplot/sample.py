import numpy as np
import pandas as pd
import re
import uproot

from systematic import Systematic

class Sample:
    """
    A class designed to encapsulate the data for a single sample and
    all associated functionality and metadata.

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
    _data : pd.DataFrame
        The data comprising the sample.
    _systematics : dict
        A dictionary of Systematic objects for the Sample.
    _print_sys : bool
        A boolean flag that toggles the printing of integrated
        systematic uncertainties for the sample.
    _presel_mask : pd.Series
        A mask to apply to the sample data for pre-selection.
    """
    def __init__(self, name, rf, category_branch, key, exposure_type,
                 trees, systematics=None, override_exposure=None, precompute=None,
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

        if override_exposure is not None:
            self.override_exposure(override_exposure, exposure_type)

        self._data = pd.concat([self._file_handle[tree].arrays(library='pd') for tree in trees])
        if self._category_branch not in self._data.columns:
            raise ValueError(f'Category branch `{self._category_branch}` not found in sample `{self._name}`.')
        if override_category is not None:
            self._data[self._category_branch] = override_category
        
        if precompute is not None:
            for k, v in precompute.items():
                self._data[k] = self._data.eval(v)
        
        if presel is not None:
            self._presel_mask = self._data.eval(presel)
            self._data = self._data[self._presel_mask]
        else:
            self._presel_mask = np.ones(len(self._data), dtype=bool)

        # Check category branch for NaNs
        if np.isnan(self._data[self._category_branch]).any():
            occurrences = len(self._data[self._data[self._category_branch].isna()])
            print(f'Found NaN category in Sample `{self._name}` with {occurrences} occurrence(s). Masking NaNs...')
            self._data = self._data[~self._data[self._category_branch].isna()]

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
            self._data['weight'] = 1
        elif self._exposure_type == 'pot':
            scale = target._exposure_pot / self._exposure_pot
        else:
            scale = target._exposure_livetime / self._exposure_livetime

        print(f"Setting weight for {self._name} to {scale:.2e}")
        self._data['weight'] = scale
        for syst in self._systematics.values():
            syst.set_weight(scale)        

    def get_data(self, variables, with_mask=None) -> dict:
        """
        Returns the data for the given variable(s) in the sample. The
        data is returned as a dictionary with the category as the key
        and the data for the requested variable as the value.

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
            the data (a pandas Series) as the value.
        weights : dict
            The weights for the requested variable in the sample. The
            weights are stored as a dictionary with the category as the
            key and the weights (a pandas Series) as the value.
        """
        data = {}
        weights = {}
        if with_mask is not None:
            mask = self._data.eval(with_mask).to_numpy(dtype=bool)
        else:
            mask = np.ones(len(self._data), dtype=bool)
        for category in np.unique(self._data[self._category_branch]):
            if np.isnan(category):
                occurrences = len(self._data[self._data[self._category_branch].isna()])
                print(f'Found NaN category in {self._name} ({occurrences} occurrences). Masking NaNs...')
                continue
            data[int(category)] = list()
            for v in variables:
                data[int(category)].append(self._data[((self._data[self._category_branch] == category) & mask)][v])  
            if 'weight' in self._data.columns:  
                weights[int(category)] = self._data[((self._data[self._category_branch] == category) & mask)]['weight']
            else:
                weights[int(category)] = None
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
        Evaluates the given formula for the sample data.

        Parameters
        ----------
        formula : str
            The formula to evaluate.

        Returns
        -------
        result : pd.Series
        """
        return self._data.eval(formula)