import numpy as np
import pandas as pd
import uproot

from variable import Variable

class Systematic:
    """
    A class designed to encapsulate the systematic uncertainty on a Sample
    due to a single source of uncertainty. A few types of systematic
    uncertainties are supported:
    - Statistical uncertainty: The uncertainty due to the finite size of
        the sample. This is calculated using a simple sqrt(N)/N formula.
    - Multisim uncertainty: The uncertainty due to the variation of
        parameters in the simulation. This is calculated by generating
        a set of universes and calculating the covariance matrix for
        the systematic parameter. The weights for the universes are
        stored in a TTree by code external to this tool.
    - Multisigma uncertainty: The uncertainty due to the variation of
        parameters in the simulation. This is calculated by generating
        a set of universes and calculating the covariance matrix for
        the systematic parameter. The weights are stored in a TTree by
        code external to this tool as "n-sigma" weights. This class
        interpolates the weights at the desired sigma level to generate
        a set of universe weights.

    The Systematic object is designed to be used in conjunction with the
    Sample object. The Sample object contains the data for the analysis
    and critically the exposure information. All systematic
    uncertainties must be rescaled to the target exposure before being
    used in an analysis. This class provides a method to rescale the
    covariance matrices to ensure that the systematic uncertainties are
    properly treated.

    Attributes
    ----------
    _name : str
        The name of the systematic uncertainty.
    _label : str
        The label of the systematic uncertainty. Intended to be used in
        plots for identifying the systematic parameter.
    _handle : uproot.models.TBranch.Model_TBranchElement
        The handle to the branch containing the weights for the
        systematic parameter.
    _variables : dict
        The Variable objects to be used for the calculation of the
        impact of the systematic uncertainty. The keys are the names of
        the variables and the values are the Variable objects
        themselves.
    _universe_weights : numpy.ndarray
        An array of weights for the systematic uncertainty. This has
        shape (nevents, nuniv) where nevents is the number of events in
        the sample and nuniv is the number of universes generated. This
        is not meant to be kept around for long periods of time due to
        memory concerns and is cleared internally by the `clear`
        method.
    _covariances : dict
        A dictionary containing the covariance matrices for each of the
        variables. The keys are the names of the variables and the values
        are numpy.ndarray objects representing the covariance matrices.
    _std : float
        The one-bin fractional uncertainty for the systematic parameter.
        This intuitively can be thought of as the uncertainty on a simple
        count of selected candidates due to the influence of the
        systematic parameter.

    Methods
    -------
    register_variable(variable)
        Register a Variable object with the Systematic object.
    process(sample, nuniv=1000)
        Processes the systematic uncertainty for the given sample for
        all configured Variables.
    get_covariance(variable)
        Retrieve the covariance matrix for the given variable.
    set_weight(weight)
        Rescales the covariance matrices using the given weight.
    combine(systematics, name, label)
        Combine a list of Systematic objects into a single Systematic
        object.
    """
    def __init__(self, name, handle, label=None):
        """
        Initializes the Systematic object with the given name and Variable.

        Parameters
        ----------
        name : str
            The name of the systematic uncertainty.
        handle : uproot.models.TBranch.Model_TBranchElement
            The handle to the branch containing the weights for the
            systematic parameter.
        """
        self._name = name
        self._label = label
        self._handle = handle
        self._variables = dict()
        self._components = None

    def register_variable(self, variable):
        """
        Register a Variable object with the Systematic object.

        Parameters
        ----------
        variable : Variable
            The Variable object to register.
        """
        self._variables[variable._key] = variable

    def process(self, sample, mask, nuniv=1000) -> np.ndarray:
        """
        Processes the systematic uncertainty for the given sample for
        all configured Variables.

        Parameters
        ----------
        sample : Sample
            The parent Sample object containing the dataset.
        mask : pd.Series
            A mask to apply to the sample data. This reflects the
            `presel` mask applied to the sample.
        nuniv : int, optional
            The number of universes to generate. The default is 1000.

        Returns
        -------
        numpy.ndarray
            An array of weights for the systematic uncertainty.
        """
        # Check that the handle is valid. This is the "usual" case
        # where the systematic weights are stored in a TTree.
        if self._handle is not None:
            universe_weights = self.get_universe_weights(
                sample=sample,
                mask=mask,
                nuniv=nuniv,
            )
            self._universe_weights = universe_weights

            # The universe weights are used to characterize the
            # systematic uncertainty for each of the Variables by
            # method of covariance matrix. Each Variable has an
            # associated binning and field name, which are used to bin
            # the events for each set of universe weights.
            self._covariances = dict()
            for name, variable in self._variables.items():
                data = sample._data[name].to_numpy()
                bin_edges = list(variable._bin_edges.values())[0]
                bin_indices = np.digitize(data, bin_edges) - 1
                valid_indices = (bin_indices >= 0) & (bin_indices < len(bin_edges) - 1)
                bin_indices = bin_indices[valid_indices]
                
                # Universes
                histogram = np.zeros((len(bin_edges) - 1, self._universe_weights.shape[1]))
                filtered_weights = self._universe_weights[valid_indices, :]
                np.add.at(histogram, bin_indices, filtered_weights)

                # Central value
                cv_histogram = np.zeros(len(bin_edges) - 1)
                np.add.at(cv_histogram, bin_indices, 1)

                # Covariance matrix calculated with respect to the central
                # value.
                diff = histogram - cv_histogram[:, np.newaxis]
                self._covariances[f'{self._name}_{name}'] = (diff @ diff.T) / (self._universe_weights.shape[1])

                # print(f"Debug: systematic '{self._name}', variable '{name}': covariance matrix:\n{self._covariances[f'{self._name}_{name}']}")
                # print(f"Debug: cv histogram:\n {cv_histogram}, universe histogram:\n {histogram}, diff:\n {diff}, shape: {self._universe_weights.shape}, nuniv: {nuniv}")
                # print(f"Debug: numpy covariance:\n {np.cov(histogram, rowvar=True, ddof=1)}")
                # One-bin uncertainty
                diff = np.sum(diff, axis=0)
                self._std = np.sqrt((diff @ diff.T) / (self._universe_weights.shape[1]))
                self._std /= self._universe_weights.shape[0]

            self._universe_weights = None
        # The handle is None, which is taken to be the case where we
        # calculate the statistical uncertainty.
        else:
            self._covariances = dict()
            for name, variable in self._variables.items():
                data = sample._data[name].to_numpy()
                bin_edges = list(variable._bin_edges.values())[0]
                bin_indices = np.digitize(data, bin_edges) - 1
                valid_indices = (bin_indices >= 0) & (bin_indices < len(bin_edges) - 1)
                bin_indices = bin_indices[valid_indices]
                
                histogram = np.zeros(len(bin_edges) - 1)
                np.add.at(histogram, bin_indices, 1)

                self._covariances[f'{self._name}_{name}'] = np.diag(histogram)
                self._std = np.sqrt(histogram.sum())
                self._std /= histogram.sum()

    def get_universe_weights(self, sample, mask, nuniv=1000, seed=None) -> np.ndarray:
        """
        Return per-event universe weights after applying a boolean mask.

        Parameters
        ----------
        sample : Sample
            Parent Sample object (unused except for alignment checks).
        mask : array-like of bool
            Boolean mask over sample._data rows.
        nuniv : int, optional
            Number of universes to generate for multisigma inputs.
        seed : int, optional
            Random seed for reproducibility.

        Returns
        -------
        numpy.ndarray
            Array of shape (N_selected, nuniv) with per-event universe weights.
        """
        # print(f"Getting universe weights for systematic '{self._name}' from branch {self._handle} with mask of length {len(mask)} and nuniv={nuniv}.")
        if self._handle is None:
            raise ValueError(
                f"Systematic '{self._name}' has no handle; cannot build universe weights."
            )

        mask = np.asarray(mask, dtype=bool)

        # Read raw weights from the TTree branch.
        # IMPORTANT: `mask` is defined over rows of sample._data. Since sample._data
        # may have been filtered after loading, we first align the weights to the
        # current DataFrame rows using its (integer) index, then apply `mask`.
        weights_all = np.stack(self._handle.array(library='np'))

        if sample is not None and hasattr(sample, '_data'):
            try:
                row_idx = sample._data.index.to_numpy()
                # If indices are valid integer positions into the raw array, align.
                if row_idx.ndim == 1 and row_idx.size > 0 and np.issubdtype(row_idx.dtype, np.integer):
                    if row_idx.max() < weights_all.shape[0] and row_idx.min() >= 0:
                        weights_all = weights_all[row_idx, :]
            except Exception:
                pass

        weights_array = weights_all[mask, :]  # shape (N_sel, Uraw)

        # Multisigma interpolation if shape is 6 or 7
        if weights_array.shape[1] in (6, 7):
            if weights_array.shape[1] == 6:
                # GENIE multisigma: -1, +1, -2, +2, -3, +3
                sigma_levels_raw = np.array([-1, 1, -2, 2, -3, 3], dtype=float)
                order = np.argsort(sigma_levels_raw)
                sigma_levels = sigma_levels_raw[order]
                W = np.asarray(weights_array, dtype=float)[:, order]
            else:
                # Detector multisigma: -3, -2, -1, 0, +1, +2, +3
                sigma_levels = np.linspace(-3, 3, 7)
                W = np.asarray(weights_array, dtype=float)

            rng = np.random.default_rng(seed)
            # Draw per-event random sigma values for each universe
            random_sigmas = rng.normal(0.0, 1.0, size=(W.shape[0], nuniv))

            # Interpolate each event's 6/7-point curve at its own random_sigmas
            universe_weights = np.empty((W.shape[0], nuniv), dtype=float)
            for i in range(W.shape[0]):
                universe_weights[i] = np.interp(random_sigmas[i], sigma_levels, W[i])

            return universe_weights

        # Multisim or already-universe weights
        return np.asarray(weights_array, dtype=float)

    def efficiency_covariance(
        self,
        sample,
        x_key,
        mask_den,
        mask_num,
        bin_edges,
        nuniv=10,
        seed=12345,
        empty_value=0.0,
        return_yields=False,
    ) -> np.ndarray:
        """
        Compute per-bin efficiency covariance for this Systematic:

            e_b^(u) = N_b^(u) / D_b^(u)

        where N_b^(u) and D_b^(u) are numerator/denominator yields in bin b
        for universe u, built from the same per-event universe weights.

        Parameters
        ----------
        sample : Sample
            Parent Sample object containing the dataset.
        x_key : str
            Column name in sample._data used for binning (the efficiency variable).
        mask_den : array-like of bool
            Denominator mask over sample._data rows.
        mask_num : array-like of bool
            Numerator mask over sample._data rows (must be subset of mask_den).
        bin_edges : array-like
            Bin edges for the efficiency variable (length B+1).
        nuniv : int, optional
            Number of universes to generate for multisigma inputs.
        seed : int, optional
            Random seed for reproducibility.
        empty_value : float, optional
            Efficiency value to assign in bins with zero denominator.

        Returns
        -------
        numpy.ndarray
            Efficiency covariance matrix of shape (B, B).
        """
        # Combined systematics (created via Systematic.combine) have no branch handle.
        # For efficiency, we combine by summing component efficiency covariances.
        if self._handle is None and getattr(self, '_components', None):
            if return_yields:
                raise ValueError("Combined systematic cannot return Num/Den; call efficiency_covariance(return_yields=True) on its components.")
            cov_total = None
            for comp in self._components:
                cov_i = comp.efficiency_covariance(
                    sample=sample,
                    x_key=x_key,
                    mask_den=mask_den,
                    mask_num=mask_num,
                    bin_edges=bin_edges,
                    nuniv=nuniv,
                    seed=seed,
                    empty_value=empty_value,
                )
                cov_total = cov_i if cov_total is None else (cov_total + cov_i)

            if cov_total is None:
                raise ValueError(
                    f"Combined systematic '{self._name}' has no components; cannot compute efficiency covariance."
                )
            return cov_total

        if self._handle is None:
            raise ValueError(
                f"Systematic '{self._name}' has no handle; cannot compute efficiency covariance."
            )

        mask_den = np.asarray(mask_den, dtype=bool)
        mask_num = np.asarray(mask_num, dtype=bool)

        if mask_den.shape[0] != len(sample._data):
            raise ValueError("mask_den length does not match sample._data")
        if mask_num.shape[0] != len(sample._data):
            raise ValueError("mask_num length does not match sample._data")

        # Ensure numerator is subset of denominator
        if not np.all(~mask_num | mask_den):
            raise ValueError("mask_num must be a subset of mask_den for efficiency.")

        bin_edges = np.asarray(bin_edges, dtype=float)
        B = len(bin_edges) - 1
        if B <= 0:
            print("Warning: No bins defined for efficiency covariance calculation.")
            return np.zeros((0, 0), dtype=float)

        # Base mask: events that participate in the denominator
        base = mask_den.copy()

        # Variable values for base events
        x_all = sample._data[x_key].to_numpy()
        x = x_all[base]

        # Universe weights for base events only
        # print(f"Getting universe weights for efficiency covariance of systematic '{self._name}' with base mask of length {base.sum()}, nuniv={nuniv}, and seed {seed}.")
        W = self.get_universe_weights(sample, mask=base, nuniv=nuniv, seed=seed)
        # W has shape (N_base, U); we want universes as "rows"
        N_base, U = W.shape

        # Reduce masks to base
        den_base = mask_den[base]
        num_base = mask_num[base]

        # Bin indices for base events
        b_idx = np.digitize(x, bin_edges) - 1
        in_range = (b_idx >= 0) & (b_idx < B)

        # Per-universe numerator and denominator in each bin
        Num = np.zeros((U, B), dtype=float)
        Den = np.zeros((U, B), dtype=float)

        for bi in range(B):
            sel_bin = (b_idx == bi) & in_range

            sel_den = sel_bin & den_base
            sel_num = sel_bin & num_base

            if np.any(sel_den):
                # Sum over selected events (axis=0) -> shape (U,)
                Den[:, bi] = W[sel_den].sum(axis=0)
            if np.any(sel_num):
                Num[:, bi] = W[sel_num].sum(axis=0)
            # print(f"Debug: bin {bi}: {np.sum(sel_num)} numerator events, {np.sum(sel_den)} denominator events.")


        # Per-universe efficiency in each bin
        E_u = np.full((U, B), float(empty_value), dtype=float)
        good = Den > 0.0
        # print(f"Debug: efficiency numerator and denominator sums per universe and bin:\nNum:\n{Num[good]}\nDen:\n{Den[good]}")
        E_u[good] = Num[good] / Den[good]

        # Build CV efficiency histogram (same spirit as `process`, i.e. diff wrt CV)
        den_cv = np.bincount(b_idx[den_base & in_range], minlength=B).astype(float)
        num_cv = np.bincount(b_idx[num_base & in_range], minlength=B).astype(float)
        E_cv = np.full(B, float(empty_value), dtype=float)
        good_cv = den_cv > 0.0
        E_cv[good_cv] = num_cv[good_cv] / den_cv[good_cv]

        # Hand covariance around CV efficiency: Cov = (diff^T diff) / U
        diff = E_u - E_cv[np.newaxis, :]
        Cov = (diff.T @ diff) / float(U)

        # print(f"Debug: cv efficiency per bin:\n{E_cv}")
        # print("Systematic efficiency covariance matrix (CV-centered, hand-built):\n", Cov)

        if return_yields:
            return Cov, Num, Den
        return Cov

    def get_covariance(self, variable) -> np.ndarray:
        """
        Retrieve the covariance matrix for the given variable.

        Parameters
        ----------
        variable : str
            The name of the variable to retrieve the covariance matrix
            for.

        Returns
        -------
        numpy.ndarray
            The covariance matrix for the given variable.
        """
        return self._covariances[f'{self._name}_{variable}']

    def set_weight(self, weight):
        """
        Rescales the covariance matrices using the given weight. This
        is necessary for cases where the analysis is comprised of an
        ensemble of samples with different exposures. The proper
        treatment in this case is to "divide out" the exposure carried
        by the parent Sample, then rescale the covariance matrices by
        the exposure of the ordinate sample. This ratio is exactly the
        `weight` parameter. Because these are covariance matrices, the
        applied rescaling is the square of the weight.

        Parameters
        ----------
        weight : float
            The weight for the systematic uncertainty.

        Returns
        -------
        None.
        """
        for kvar, vvar in self._variables.items():
            self._covariances[f'{self._name}_{kvar}'] *= weight**2

    @staticmethod
    def combine(systematics, name, label) -> 'Systematic':
        """
        Combine a list of Systematic objects into a single Systematic
        object. This is done by adding the covariance matrices of the
        Systematic objects. The underlying assumption is that the
        systematic parameters encapsulated by the Systematic objects are
        totally uncorrelated.

        Parameters
        ----------
        systematics : list
            A list of Systematic objects to combine.
        name : str
            The name of the new Systematic object.
        label : str
            The label of the new Systematic object.

        Returns
        -------
        Systematic
            A new Systematic object with the covariance matrices added.
        """
        new_systematic = Systematic(name, None, label)
        # Retain the component systematics for quantities that require
        # per-event universe weights (e.g. efficiency covariances).
        new_systematic._components = list(systematics)
        new_systematic._variables = systematics[0]._variables
        new_systematic._covariances = dict()
        for kvar, vvar in new_systematic._variables.items():
            if any([kvar not in sys._variables.keys() for sys in systematics]):
                msg = f'Variable {kvar} not found in all Systematic objects.'
                raise ValueError(msg)
            new_systematic._covariances[f'{name}_{kvar}'] = np.sum(
                [sys._covariances[f'{sys._name}_{kvar}'] for sys in systematics],
                axis=0
            )
            new_systematic._std = np.sqrt(np.sum([sys._std**2 for sys in systematics]))
        
        return new_systematic

    @staticmethod
    def transform_as(cov, param):
        """
        Apply a scale correction to the covariance matrix. This is
        intended to be used when the covariance matrix is meant to be
        rescaled to a different exposure, or for the case where the
        spectrum is being area normalized. The former case is realized
        by a "simple" scale correction (i.e. the covariance matrix is
        multiplied by the square of the scale). The latter case is
        realized by a "normalized" scale correction (i.e. the
        covariance matrix is properly transformed to account for how
        the area normalization impacts the covariance matrix).

        Parameters
        ----------
        cov : numpy.ndarray
            The covariance matrix to transform.
        param : float or np.ndarray
            The parameter to use for the transformation. For a "simple"
            scale correction, this is a float. For a "normalized" scale
            correction, this is a numpy array with the same shape as
            `cov` and containing the bin contents of the histogram.
        
        Returns
        -------
        numpy.ndarray
            The transformed covariance matrix.
        """
        if np.isscalar(param):
            return cov * param**2

        else:
            A = np.sum(param)
            y = param.reshape(-1, 1)
            delta = np.eye(len(y))
            jacobian = (delta * A - y) / A**2
            return jacobian @ cov @ jacobian.T

    def __repr__(self):
        s = f'--Systematic({self._name}, {self._label})--'
        s += f'\n\tFractional uncertainty: {self._std:.2%}'
        return s

    @property
    def name(self):
        return self._name
    
    @property
    def label(self):
        return self._label

    @property
    def std(self):
        return self._std
