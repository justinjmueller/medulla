import os
import toml
import uproot
from matplotlib import pyplot as plt

# TODO: Make imports consistent (absolute vs relative)
from .core.sample import Sample
from .core.figure import SpineFigure, SimpleFigure
from .artists.spectra1d import SpineSpectra1D
from .artists.spectra2d import SpineSpectra2D
from .artists.efficiency import SpineEfficiency
from .artists.confusion import ConfusionMatrix
from .artists.roc import ROCCurve
from .artists.ternary import Ternary
from .core.style import Style
from .core.variable import Variable

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
        # TODO: Can all the `.keys()` calls be removed for config
        # checks?
        self._toml_path = toml_path
        self._config = toml.load(self._toml_path)
        for table in self._config.get('this_includes', []):
            Analysis.handle_include(self._config, table)
        rf = uproot.open(rf_path)

        # Load the output path
        if 'output' not in self._config.keys():
            raise ConfigException(
                f"No output path defined in the TOML file. Please"
                f" check for a valid output configuration block in the"
                f" TOML file ('{toml_path}')."
            )
        # TODO: Switch to pathlib?
        self._output_path = self._config['output']['path']

        # Load the categories table
        self._categories = dict()
        for ci, cat in enumerate(self._config['analysis']['category_assignment']):
            self._categories.update({c : self._config['analysis']['category_labels'][ci] for c in cat})
        self._colors = {c : self._config['analysis']['category_colors'][ci] for ci, c in enumerate(self._config['analysis']['category_labels'])}
        self._category_types = {c : self._config['analysis']['category_types'][ci] for ci, c in enumerate(self._config['analysis']['category_labels'])}

        # Initialize the samples for the analysis. The configuration
        # file must contain a valid sample configuration block, failing
        # which an exception is raised.
        self.initialize_samples(rf)

        # Initialize the styles for the analysis. A default style is
        # always provided, so no exception is raised if no style block
        # is present in the configuration file.
        self.initialize_styles()

        # Initialize the variables for the analysis. The configuration
        # file must contain a valid variable configuration block,
        # failing which an exception is raised.
        self.initialize_variables()

        # Register variables with samples
        if 'systematic_recipe' in self._config.keys():
            recipes = self._config['systematic_recipe']
        else:
            recipes = list()
        # TODO: Re-enable this functionality after refactoring the
        # variable registration and systematic processing.
        #for s in self._samples.values():
            #for v in self._variables.values():
            #    s.register_variable(v, self._categories)
            #s.process_systematics(recipes)

        # Initialize the figures and artists for the analysis. The
        # configuration file must contain a valid figure configuration
        # block, failing which an exception is raised.
        self.initialize_figures_and_artists()

    def initialize_samples(
        self,
        rf : uproot.reading.ReadOnlyDirectory
    ) -> None:
        """
        Initializes the samples for the analysis using the provided
        configuration file.

        Parameters
        ----------
        rf : uproot.reading.ReadOnlyDirectory
            The ROOT file containing the data.

        Returns
        -------
        None.

        Raises
        ------
        ConfigException
            No sample configuration block is present in the TOML
            configuration file.
        """
        # Check if the samples table exists
        if 'samples' not in self._config.keys():
            raise ConfigException(
                f"No samples defined in the TOML file. Please check"
                f" for a valid sample configuration block in the TOML"
                f" file ('{toml_path}')."
            )
        if isinstance(self._config['samples'], list):
            # Support "list-style" sample definitions (cleaner, in my
            # personal opinion)
            self._samples = {x['name']: Sample(
                x['name'],
                rf,
                self._config['analysis']['category_branch'],
                **x
            ) for x in self._config['samples']}
        else:
            # Support "dict-style" sample definitions (legacy support
            # and not as clean in my personal opinion)
            self._samples = {name: Sample(
                name,
                rf,
                self._config['analysis']['category_branch'],
                **self._config['samples'][name]
            ) for name in self._config['samples']}

    def initialize_styles(self) -> None:
        """
        Initializes the styles for the analysis using the provided
        configuration file. Unlike the other initialization methods,
        this one provides a default style if none are defined in the
        configuration file and does not raise an exception.

        Parameters
        ----------
        None.

        Returns
        -------
        None.
        """
        # Provide a default style to simplify the user experience.
        self._styles = {
            'default' : Style(
                'default',
                style_sheet='default',
                default_figsize=(8,6),
                title='',
                mark_preliminary='Work-in-Progress',
                plot_kwargs={'histtype': 'barstacked', 'stacked': True}
            )
        }

        # Load the plot styles table (if it exists)
        if 'style' in self._config.keys():
            self._styles.update(
                {x['name'] : Style(**x) for x in self._config['style']}
            )

    def initialize_variables(self) -> None:
        """
        Initializes the variables for the analysis using the provided
        configuration file.

        Parameters
        ----------
        None.

        Returns
        -------
        None.

        Raises
        ------
        ConfigException
            No variable configuration block is present in the TOML
            configuration file.
        """
        # Check if the variables table exists. This is a required
        # configuration block, so an exception is raised if it is not
        # present.
        if 'variables' not in self._config.keys():
            raise ConfigException(
                f"No variables defined in the TOML file. Please check"
                f" for a valid variable configuration block"
                f" (table='variables') in the TOML file"
                f" ('{toml_path}').")
        
        # Load the variables. Two styles of variable definitions
        # are supported for user convenience.
        if isinstance(self._config['variables'], list):
            # Support "list-style" variable definitions (cleaner, in
            # my personal opinion)
            self._variables = {x['name']: Variable(
                **x
            ) for x in self._config['variables']}
        else:
            # Support "dict-style" variable definitions (legacy support
            # and not as clean in my personal opinion)
            self._variables = {name: Variable(
                name,
                **self._config['variables'][name]
            ) for name in self._config['variables']}

        # Now "finalize" the Variable by checking its validity (well-
        # configured and present in all samples) and setting up its
        # binning scheme. It's important to note that Variables failing
        # the validity check do not raise an exception here; instead,
        # the exception is raised when the Variable is used in an
        # artist later on. Such exceptions would not actually protect
        # the user.
        for v in self._variables.values():
            v.finalize(
                list(self._samples.values()),
                self._categories
            )

    def initialize_figures_and_artists(self) -> None:
        """
        Initializes the figures and artists for the analysis using the
        provided configuration file.

        Parameters
        ----------
        None.

        Returns
        -------
        None.
        
        Raises
        ------
        ConfigException
            No figure configuration block is present in the TOML
            configuration file.
        """
        # Load the artists table
        if 'figure' not in self._config.keys():
            raise ConfigException(
                f"No figures defined in the TOML file. Please check"
                f" for a valid figure configuration block"
                f" (table='figure') in the TOML file ('{toml_path}').")
        self._figures = dict()
        self._artists = list()

        for figure in self._config['figure']:
            # Provide a default style if none is specified
            if 'style' not in figure.keys():
                figure['style'] = 'default'
            
            # Create the figure with the specified style using a
            # context manager to ensure proper style application
            with self._styles[figure['style']] as style:
                if figure['type'] == 'SimpleFigure':
                    self._figures[figure['name']] = SimpleFigure(
                        figure.get('figsize', style.default_figsize),
                        style,
                        figure.get('title', style.default_title)
                    )
                else:
                    raise ConfigException(
                        f"Figure type '{figure['type']}' does not"
                        f" correspond to a valid figure type."
                    )
                
                # Register the artists for the figure. Depending on the
                # figure type, multiple artist types may be supported.
                for artist in figure['artists']:
                    self.register_artist_to_figure(
                        artist,
                        self._figures[figure['name']]
                    )

    def register_artist_to_figure(
        self,
        artcfg : dict,
        figure : SpineFigure
    ) -> None:
        """
        Registers an artist to the given figure using the provided
        configuration. This method will raise an exception if the
        artist type requested is not supported.
        
        Parameters
        ----------
        artcfg : dict
            The artist configuration dictionary.
        figure : SpineFigure
            The figure to register the artist to.
        
        Returns
        -------
        None.

        Raises
        ------
        ConfigException
            Raised in a few cases:
            - The variable(s) required for the artist are not present
              in all samples.
            - The artist type does not correspond to a valid artist
              type.
        """
        # Do a quick check that the requested artist type is valid
        valid_artists = [
            'SpineSpectra1D',
            'SpineSpectra2D',
            'SpineEfficiency',
            'ConfusionMatrix',
            'ROCCurve',
            'Ternary'
        ]
        if artcfg['type'] not in valid_artists:
            raise ConfigException(
                f"Artist type '{artcfg['type']}' does not correspond"
                f" to a valid artist type."
            )

        # Check if the artist is restricted to certain groups. This
        # allows for some additional amount of control over the
        # plotting beyond what is provided by the category assignment
        # configuration.
        restrict = {}
        group_setting = artcfg.get('groups', [])
        if group_setting:
            for g in group_setting:
                restrict.update(
                    {k : v for k,v in self._categories.items() if v==g}
                )
        else:
            restrict = self._categories.copy()

        # TODO: Ensure draw_kwargs are handled consistently across
        # all artist types and directly mapped to kwargs of the 'draw'
        # method of each artist.

        # Now, create the artist and register it to the figure. We do
        # this in a series of if-else statements, thereby making it 
        # easy to extend the functionality in the future.
        if artcfg['type'] == 'SpineSpectra1D':
            # Do a full check on the variable required for the artist,
            # and take the opportunity to update the artist dictionary
            # with the validated Variable object.
            artcfg['variable'] = self.validate_variable(
                'variable',
                artcfg
            )

            # Create the artist
            art = SpineSpectra1D(
                restrict,
                self._colors,
                self._category_types,
                **artcfg
            )
            draw_kwargs = artcfg.get('draw_kwargs', {})
            draw_kwargs['draw_error'] = draw_kwargs.get('draw_error', None)
            figure.register_spine_artist(art, draw_kwargs=draw_kwargs)
            self._artists.append(art)
        elif artcfg['type'] == 'SpineSpectra2D':
            # Do a full check on the variables required for the artist,
            # and take the opportunity to update the artist dictionary
            # with the validated Variable objects.
            artcfg['xvariable'] = self.validate_variable(
                'xvariable',
                artcfg
            )
            artcfg['yvariable'] = self.validate_variable(
                'yvariable',
                artcfg
            )

            # Create the artist
            art = SpineSpectra2D(
                restrict,
                self._colors,
                self._category_types,
                **artcfg
            )
            figure.register_spine_artist(
                art,
                draw_kwargs=artcfg.get('draw_kwargs', {})
            )
            self._artists.append(art)

        elif artcfg['type'] == 'SpineEfficiency':
            # Do a full check on the variable required for the artist,
            # and take the opportunity to update the artist dictionary
            # with the validated Variable object.
            artcfg['variable'] = self.validate_variable(
                'variable',
                artcfg
            )

            # Create the artist
            art = SpineEfficiency(
                restrict,
                **artcfg
            )
            figure.register_spine_artist(
                art,
                draw_kwargs=artcfg.get('draw_kwargs', {})
            )
            self._artists.append(art)

        elif artcfg['type'] == 'ConfusionMatrix':
            # Do a full check on the variables required for the artist,
            # and take the opportunity to update the artist dictionary
            # with the validated Variable objects.
            artcfg['true_labels'] = self.validate_variable(
                'true_labels',
                artcfg
            )
            artcfg['predicted_labels'] = self.validate_variable(
                'predicted_labels',
                artcfg
            )

            # Create the artist
            art = ConfusionMatrix(
                restrict,
                **artcfg
            )
            figure.register_spine_artist(
                art,
                draw_kwargs=artcfg.get('draw_kwargs', {})
            )
            self._artists.append(art)
        
        elif artcfg['type'] == 'ROCCurve':
            # Do a full check on the variable required for the artist,
            # and take the opportunity to update the artist dictionary
            # with the validated Variable objects.
            for disc in artcfg['discriminant_scores']:
                artcfg[disc] = disc
            artcfg['discriminant_scores'] = {
                g : self.validate_variable(
                    artcfg['discriminant_scores'][i],
                    artcfg
                ) for i, g in enumerate(restrict.keys())
            }

            # Create the artist
            art = ROCCurve(
                restrict,
                **artcfg
            )
            figure.register_spine_artist(
                art,
                draw_kwargs=artcfg.get('draw_kwargs', {})
            )
            self._artists.append(art)

        elif artcfg['type'] == 'Ternary':
            # Do a full check on the variables required for the artist,
            # and take the opportunity to update the artist dictionary
            # with the validated Variable objects.
            artcfg['var0'] = self.validate_variable(
                'var0',
                artcfg
            )
            artcfg['var1'] = self.validate_variable(
                'var1',
                artcfg
            )
            artcfg['var2'] = self.validate_variable(
                'var2',
                artcfg
            )

            # Create the artist
            art = Ternary(
                restrict,
                **artcfg
            )
            figure.register_spine_artist(
                art,
                draw_kwargs=artcfg.get('draw_kwargs', {})
            )
            self._artists.append(art)
        
        else:
            raise ConfigException(
                f"Artist type '{artcfg['type']}' does not correspond"
                f" to a valid artist type."
            )

    def validate_variable(
        self,
        vkey : str,
        artcfg: dict
    ) -> Variable:
        """
        Validates that the variable required for the artist is present
        in all samples.

        Parameters
        ----------
        vkey : str
            The key of the variable in the artist configuration.
        artcfg : dict
            The artist configuration dictionary.
        
        Returns
        -------
        var : Variable
            The validated Variable object.

        Raises
        ------
        ConfigException
            Raised in a few cases:
            - The variable required for the artist is not present in
              the artist configuration.
            - The variable required for the artist is not present in
              the variable list.
            - The variable required for the artist is not present in
              all samples (not valid).
        """
        # Check if the variable key is present in the artist
        # configuration block.
        if vkey not in artcfg:
            raise ConfigException(
                f"Artist configuration for artist of type"
                f" '{artcfg['type']}' is missing required field"
                f" '{vkey}'."
            )
        
        # Check if the variable is present in the variable list within
        # the class.
        if artcfg[vkey] not in self._variables:
            raise ConfigException(
                f"Variable '{artcfg[vkey]}' not found in"
                f" variable list when attempting to create artist"
                f" of type '{artcfg['type']}'. Please check the"
                f" variable configuration block in the TOML file"
                f" ('{self._toml_path}')."
            )

        # Now it is safe to retrieve the Variable object.
        var = self._variables[artcfg[vkey]]

        # Finally, run a validity check to ensure the variable is
        # present in all samples.
        if not all(var._valid.values()):
            miss = [k for k, v in var._valid.items() if not v]
            raise ConfigException(
                f"Variable '{var._name}' not found in all samples:"
                f" ({' '.join(miss)})."
            )
        return var

    def override_exposure(
        self,
        sample_name : str,
        exposure : float,
        exposure_type : Optional[str] = 'pot'
    ) -> None:
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
        exposure_type : Optional[str]
            The type of exposure to set. This can be either 'pot' or
            'livetime'. The default is 'pot'.
        
        Returns
        -------
        None.
        """
        # Check if the sample exists
        if sample_name not in self._samples.keys():
            raise ConfigException(
                f"Sample '{sample_name}' not found in sample list"
                f" when attempting to override exposure. Please check"
                f" the sample configuration block in the TOML file"
                f" ('{self._toml_path}')."
            )
        
        # Override the exposure for the sample
        self._samples[sample_name].override_exposure(
            exposure,
            exposure_type
        )

    def run(
        self,
        close_figs : Optional[bool] = True
    ) -> None:
        """
        Runs the analysis on the samples.

        Parameters
        ----------
        close_figs : Optional[bool]
            Whether to close the figures after saving them. The default
            is True. This is a useful toggle between two use cases:
            interactive mode (False) and batch mode (True).

        Returns
        -------
        None.
        """
        # Check if the ordinate sample exists
        ordinate_sample = self._config['analysis']['ordinate_sample']
        if ordinate_sample not in self._samples.keys():
            raise ConfigException(
                f"Ordinate sample '{ordinate_sample}' not found in"
                f" sample list. Please check the sample configuration"
                f" block (table='samples') in the TOML file"
                f" ('{self._toml_path}').")

        # Set the weights for all samples based on the ordinate sample
        ordinate = self._samples[ordinate_sample]
        for s in self._samples.values():
            s.set_weight(target=ordinate)

        # Add all samples to all artists
        for artist in self._artists:
            for sample in self._samples.values():
                artist.add_sample(sample, sample==ordinate)

        # Check if the output path exists. If not, create it.
        if not os.path.exists(self._output_path):
            os.makedirs(self._output_path)

        # Create all figures and save them to the output path
        for figname, figure in self._figures.items():
            figure.create()
            # TODO: Allow for different output formats
            figure.figure.savefig(f"{self._output_path}/{figname}.png")
            if close_figs:
                figure.close()

    def run_interactively(
        self,
        figure : str
    ) -> SpineFigure:
        """
        Runs the analysis on the samples and creates the figure. This
        method is useful for interactive plotting in a Jupyter
        notebook.

        Parameters
        ----------
        figure : str
            The name of the figure to create.
        
        Returns
        -------
        SpineFigure
            The figure object.
        """
        # Check if the ordinate sample exists
        ordinate_sample = self._config['analysis']['ordinate_sample']
        if ordinate_sample not in self._samples.keys():
            raise ConfigException(
                f"Ordinate sample '{ordinate_sample}' not found in"
                f" sample list. Please check the sample configuration"
                f" block (table='samples') in the TOML file"
                f" ('{self._toml_path}').")
        
        # Set the weights for all samples based on the ordinate sample
        ordinate = self._samples[ordinate_sample]
        for s in self._samples.values():
            s.set_weight(target=ordinate)

        # Add all samples to all artists
        for artist in self._artists:
            for sample in self._samples.values():
                artist.add_sample(sample, sample==ordinate)

        # Check if the output path exists. If not, create it.
        if not os.path.exists(self._output_path):
            os.makedirs(self._output_path)

        # Create the figure and return it
        if figure not in self._figures.keys():
            raise ConfigException(
                f"Figure '{figure}' not found in figure list. Please"
                f" check the figure configuration block"
                f" (table='figure') in the TOML file "
                f" ('{self._toml_path}').")
        self._figures[figure].create()
        return self._figures[figure].figure

    @staticmethod
    def handle_include(
        config : dict,
        table : dict
    ) -> None:
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
                        config[key].update(
                            {k: c[key][k] for k in value}
                        )
                    else:
                        config[key] = {v: c[key][v] for v in value}
            else:
                for key, value in c.items():
                    if key in config.keys():
                        config[key].update(value)
                    else:
                        config[key] = value