import os
import toml
import uproot
from matplotlib import pyplot as plt

from sample import Sample
from figure import SpineFigure, SimpleFigure
from spectra1d import SpineSpectra1D
from spectra2d import SpineSpectra2D
from efficiency import SpineEfficiency
from confusion import ConfusionMatrix
from roc import ROCCurve
from ternary import Ternary
from style import Style
from variable import Variable

class ConfigException(Exception):
    pass

class Analysis:
    """
    The Analysis class is used to containerize the plotting of the
    ensemble of samples for a variety of variables and plotting styles.

    Attributes
    ----------
    _toml_path : str
        The path to the TOML configuration file for the analysis.
    _config : dict
        The configuration dictionary loaded from the TOML file.
    _output_path : str
        The path to the output directory for the analysis.
    _ordinate_sample : str
        The name of the sample to use as the ordinate for the analysis.
    _category_tree : str
        The label of the branch corresponding to the category label in
        the input TTree.
    """
    def __init__(self, toml_path, rf_path) -> None:
        """
        Initializes the Analysis object with the given kwargs.

        Parameters
        ----------
        toml_path : str
            The path to the TOML configuration file for the analysis.
        rf_path : str
            The path to the ROOT file containing the data.

        Returns
        -------
        None
        """
        self._toml_path = toml_path
        self._config = toml.load(self._toml_path)
        for table in self._config.get('this_includes', []):
            Analysis.handle_include(self._config, table)
        rf = uproot.open(rf_path)

        # Load the output path
        if 'output' not in self._config.keys():
            raise ConfigException(f"No output path defined in the TOML file. Please check for a valid output configuration block in the TOML file ('{toml_path}').")
        self._output_path = self._config['output']['path']

        # Load the categories table
        self._categories = dict()
        for ci, cat in enumerate(self._config['analysis']['category_assignment']):
            self._categories.update({c : self._config['analysis']['category_labels'][ci] for c in cat})
        self._colors = {c : self._config['analysis']['category_colors'][ci] for ci, c in enumerate(self._config['analysis']['category_labels'])}
        self._category_types = {c : self._config['analysis']['category_types'][ci] for ci, c in enumerate(self._config['analysis']['category_labels'])}

        # Initialize the samples
        if 'samples' not in self._config.keys():
            raise ConfigException(f"No samples defined in the TOML file. Please check for a valid sample configuration block in the TOML file ('{toml_path}').")
        self._samples = {name: Sample(name, rf, self._config['analysis']['category_branch'], **self._config['samples'][name]) for name in self._config['samples']}

        # Load the plot styles table
        if 'style' not in self._config.keys():
            raise ConfigException(f"No plot styles defined in the TOML file. Please check for a valid plot style configuration block (table='styles') in the TOML file ('{toml_path}').")
        self._styles = {x['name'] : Style(**x) for x in self._config['style']}

        # Load the variable table
        if 'variables' not in self._config.keys():
            raise ConfigException(f"No variables defined in the TOML file. Please check for a valid variable configuration block (table='variables') in the TOML file ('{toml_path}').")
        self._variables = {name: Variable(name, **self._config['variables'][name]) for name in self._config['variables']}

        # Register variables with samples
        if 'systematic_recipe' in self._config.keys():
            recipes = self._config['systematic_recipe']
        else:
            recipes = list()
        for s in self._samples.values():
            for v in self._variables.values():
                s.register_variable(v, self._categories)
            s.process_systematics(recipes)

        # Load the artists table
        if 'figure' not in self._config.keys():
            raise ConfigException(f"No figures defined in the TOML file. Please check for a valid figure configuration block (table='figure') in the TOML file ('{toml_path}').")
        self._figures = dict()
        self._artists = list()
        for fig in self._config['figure']:
            if fig['type'] == 'SimpleFigure':
                with self._styles[fig['style']] as style:
                    self._figures[fig['name']] = SimpleFigure(fig.get('figsize', style.default_figsize), style, fig.get('title', style.default_title))
                    for x in fig['artists']:

                        # Check if the artist is restricted to certain
                        # groups. This allows for some additional
                        # amount of control over the plotting.
                        restrict_categories = {}
                        group_setting = x.get('groups', [])
                        if group_setting:
                            for g in group_setting:
                                restrict_categories.update({k : v for k,v in self._categories.items() if v == g})
                        else:
                            restrict_categories = self._categories.copy()
                    
                        if x['type'] == 'SpineSpectra1D':
                            # Check if the variable is present in all samples
                            if not all(self._variables[x['variable']]._validity_check.values()):
                                missing_samples = [k for k, v in self._variables[x['variable']]._validity_check.items() if not v]
                                raise ConfigException(f"Variable '{x['variable']}' not found in all samples ({' '.join(missing_samples)}).")
                            
                            # Create the artist
                            art = SpineSpectra1D(self._variables[x['variable']], restrict_categories,
                                                 self._colors, self._category_types, x.get('title', None),
                                                 x.get('xrange', None), x.get('xtitle', None),
                                                 x.get('yrange', None), x.get('ytitle', None))
                            draw_kwargs = x.get('draw_kwargs', {})
                            draw_kwargs['draw_error'] = draw_kwargs.get('draw_error', None)
                            self._figures[fig['name']].register_spine_artist(art, draw_kwargs=draw_kwargs)
                            self._artists.append(art)
                        elif x['type'] == 'SpineSpectra2D':
                            # Check if the variables are present in all samples
                            if not all(self._variables[x['xvariable']]._validity_check.values()) or not all(self._variables[x['yvariable']]._validity_check.values()):
                                missing_samples = [k for k, v in self._variables[x['xvariable']]._validity_check.items() if not v] + [k for k, v in self._variables[x['yvariable']]._validity_check.items() if not v]
                                raise ConfigException(f"Variable '{x['xvariable']}' or '{x['yvariable']}' not found in all samples ({' '.join(missing_samples)}).")
                            
                            # Create the artist
                            art = SpineSpectra2D([self._variables[x['xvariable']], self._variables[x['yvariable']]],
                                                  restrict_categories, self._colors, self._category_types,
                                                  x.get('title', None), x.get('xrange', None), x.get('xtitle', None),
                                                  x.get('yrange', None), x.get('ytitle', None))
                            self._figures[fig['name']].register_spine_artist(art, draw_kwargs=x.get('draw_kwargs', {}))
                            self._artists.append(art)
                        elif x['type'] == 'SpineEfficiency':
                            # Check if the variable is present in all samples
                            if not all(self._variables[x['variable']]._validity_check.values()):
                                missing_samples = [k for k, v in self._variables[x['variable']]._validity_check.items() if not v]
                                raise ConfigException(f"Variable '{x['variable']}' not found in all samples ({' '.join(missing_samples)}).")
                            
                            # Grab artist settings
                            show_option = x.get('draw_kwargs', {}).get('show_option', 'table')
                            npts = x.get('draw_kwargs', {}).get('npts', 1e6)
                            
                            # Create the artist
                            art = SpineEfficiency(self._variables[x['variable']], restrict_categories,
                                                  x['cuts'], x.get('title', None), x.get('xrange', None),
                                                  x.get('xtitle', None), show_option, npts)
                            self._figures[fig['name']].register_spine_artist(art, draw_kwargs=x.get('draw_kwargs', {}))
                            self._artists.append(art)

                        elif x['type'] == 'ConfusionMatrix':
                            # Check if the true and predicted labels are present in all samples
                            if not all([self._variables[x['true_labels']]._validity_check.values(),
                                        self._variables[x['predicted_labels']]._validity_check.values()]):
                                missing_samples = [k for k, v in self._variables[x['true_labels']]._validity_check.items() if not v] + \
                                                  [k for k, v in self._variables[x['predicted_labels']]._validity_check.items() if not v]
                                raise ConfigException(f"Variable '{x['true_labels']}' or '{x['predicted_labels']}' not found in all samples ({' '.join(missing_samples)}).")
                            
                            # Create the artist
                            kwargs = {k : v for k, v in x.items() if k not in ['type', 'true_labels', 'predicted_labels', 'groups', 'draw_kwargs']}
                            art = ConfusionMatrix(restrict_categories, self._variables[x['true_labels']],
                                                  self._variables[x['predicted_labels']], **kwargs)
                            self._figures[fig['name']].register_spine_artist(art, draw_kwargs=x.get('draw_kwargs', {}))
                            self._artists.append(art)

                        elif x['type'] == 'ROCCurve':
                            # Check if the discriminant scores are present in all samples
                            if not all([self._variables[disc]._validity_check.values() for disc in x['discriminant_scores']]):
                                missing_samples = [k for k, v in self._variables[disc]._validity_check.items() if not v for disc in x['discriminant_scores']]
                                raise ConfigException(f"Variable '{disc}' not found in all samples ({' '.join(missing_samples)}).")
                            
                            # Create the artist
                            disc = {g : self._variables[x['discriminant_scores'][i]] for i, g in enumerate(x['groups'])}
                            art = ROCCurve(restrict_categories, x['labels'], x['pos_label'], disc,
                                           x.get('background', None), x.get('title', None))
                            self._figures[fig['name']].register_spine_artist(art, draw_kwargs=x.get('draw_kwargs', {}))
                            self._artists.append(art)

                        elif x['type'] == 'Ternary':
                            # Check if the variables are present in all samples
                            if not all([self._variables[x[v]]._validity_check.values() for v in ['var0', 'var1', 'var2']]):
                                missing_samples = [k for k, v in self._variables[x[v]]._validity_check.items() for v in ['var0', 'var1', 'var2']]
                                raise ConfigException(f"Variable '{v}' not found in all samples ({' '.join(missing_samples)}).")

                            # Create the artist
                            art = Ternary(self._variables[x['var0']], self._variables[x['var1']], self._variables[x['var2']],
                                          restrict_categories, x.get('title', None))
                            self._figures[fig['name']].register_spine_artist(art, draw_kwargs=x.get('draw_kwargs', {}))
                            self._artists.append(art)
                            

    def override_exposure(self, sample_name, exposure, exposure_type='pot') -> None:
        """
        Overrides the exposure for the given sample. This is useful for
        setting the exposure for samples for which the exposure is not
        valid. The exposure type can be either 'pot' or 'livetime'. It
        is not recommended to use this method unless the exposure is
        known to be incorrect.
        
        Parameters
        ----------
        sample_name : str
            The name of the sample to override.
        exposure : float
            The exposure to set for the sample.
        exposure_type : str
            The type of exposure to set. This can be either 'pot' or
            'livetime'. The default is 'pot'.
        
        Returns
        -------
        None.
        """
        if sample_name not in self._samples.keys():
            raise ConfigException(f"Sample '{sample_name}' not found in sample list when attempting to override exposure. Please check the sample configuration block in the TOML file ('{self._toml_path}').")
        self._samples[sample_name].override_exposure(exposure, exposure_type)

    def run(self, close_figs=True) -> None:
        """
        Runs the analysis on the samples.

        Parameters
        ----------
        close_figs : bool
            Whether to close the figures after saving them. The default
            is True. This is a useful toggle between two use cases:
            interactive mode (False) and batch mode (True).

        Returns
        -------
        None.
        """
        if self._config['analysis']['ordinate_sample'] not in self._samples.keys():
            raise ConfigException(f"Ordinate sample '{self._config['analysis']['ordinate_sample']}' not found in sample list. Please check the sample configuration block (table='samples') in the TOML file ('{self._toml_path}').")
        ordinate = self._samples[self._config['analysis']['ordinate_sample']]
        for s in self._samples.values():
            s.set_weight(target=ordinate)

        for artist in self._artists:
            for sample in self._samples.values():
                artist.add_sample(sample, sample==ordinate)

        # Check if the output path exists. If not, create it.
        if not os.path.exists(self._output_path):
            os.makedirs(self._output_path)
        for figname, figure in self._figures.items():
            figure.create()
            figure.figure.savefig(f"{self._output_path}/{figname}.png")
            if close_figs:
                figure.close()

    def run_interactively(self, figure) -> SpineFigure:
        """
        Runs the analysis on the samples and creates the figure. This
        method is useful for interactive plotting in a Jupyter notebook.

        Parameters
        ----------
        figure : str
            The name of the figure to create.
        
        Returns
        -------
        SpineFigure
            The figure object.
        """
        if self._config['analysis']['ordinate_sample'] not in self._samples.keys():
            raise ConfigException(f"Ordinate sample '{self._config['analysis']['ordinate_sample']}' not found in sample list. Please check the sample configuration block (table='samples') in the TOML file ('{self._toml_path}').")
        ordinate = self._samples[self._config['analysis']['ordinate_sample']]
        for s in self._samples.values():
            s.set_weight(target=ordinate)

        for artist in self._artists:
            for sample in self._samples.values():
                artist.add_sample(sample, sample==ordinate)

        self._figures[figure].create()
        return self._figures[figure].figure

    @staticmethod
    def handle_include(config, table):
        """
        Handles the inclusion of other configuration files in the main
        configuration file. The include directive may also contain some
        other optional fields which correspond to specific actions:

        - choose: a dictionary of key-value pairs where the key is the
        name of a table and the value is a list of corresponding sub-
        tables to include.

        Parameters
        ----------
        config : dict
            The configuration dictionary to update.
        table : dict
            The block representing the include directive.
        
        Returns
        -------
        None.
        """
        with open(table['file'], 'r') as f:
            c = toml.load(f)
            if 'choose' in table.keys():
                for key, value in table['choose'].items():
                    if key in config.keys():
                        config[key].update({k: c[key][k] for k in value})
                    else:
                        config[key] = {v: c[key][v] for v in value}
            else:
                for key, value in c.items():
                    if key in config.keys():
                        config[key].update(value)
                    else:
                        config[key] = value