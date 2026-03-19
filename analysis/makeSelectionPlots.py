import matplotlib.pyplot as plt
import os
import sys
import pickle
os.environ['MEDULLA_PLOT_DIR'] =  '/nashome/m/micarrig/icarus/nuESpine/medulla/spineplot/'
sys.path.append('/nashome/m/micarrig/icarus/nuESpine/medulla/spineplot/')
from efficiency import SpineEfficiency
from analysis import Analysis

# Make the plot based on the provided variable, and save in dir
def makePlot(var, dir):

    fig = ana.run_interactively(var)
    fig.set_size_inches(10,8)
    for ax in fig.axes:
        leg = ax.get_legend()
        if leg is not None:
            leg.get_frame().set_facecolor('white')
            leg.get_frame().set_alpha(0.2)
    plt.savefig(f'{dir}/plot_{var}.png', dpi=300, bbox_inches='tight')
    with open(f'{dir}/plot_{var}.pkl', 'wb') as f:
        pickle.dump(fig, f)
    plt.close()

if __name__ == "__main__":

    #ana = Analysis('/nashome/m/micarrig/icarus/nuESpine/medulla/analysis/nue_plots.toml', 'root://fndcadoor.fnal.gov:/icarus/scratch/users/micarrig/nueCCInclusiveSpine14/merged_output.root')
    ana = Analysis('/nashome/m/micarrig/icarus/nuESpine/medulla/analysis/nue_plots.toml', '/nashome/m/micarrig/icarus/nuESpine/gundam_nueCCInclusive.root')

    plot_category = 'neutrino_type'

    plotVars = [
        'true_neutrino_energy',
        'reco_leading_electron_energy',
        'reco_ele_beam_cosTheta',
        'reco_leading_electron_start_dedx',
        'reco_leading_electron_axial_spread',
        'reco_leading_electron_directional_spread',
        'reco_leading_electron_start_x',
        'reco_leading_electron_start_y',
        'reco_leading_electron_start_z',
        'reco_leading_electron_start_dir_x',
        'reco_leading_electron_start_dir_y',
        'reco_leading_electron_start_dir_z',
        'reco_leading_electron_mip_softmax',
        'reco_leading_electron_photon_softmax',
        'reco_leading_electron_muon_softmax',
        'reco_leading_electron_proton_softmax',
        'reco_leading_electron_pion_softmax',
        'reco_leading_electron_electron_softmax',
        'reco_leading_electron_primary_softmax',
        'reco_leading_electron_secondary_softmax',
        'reco_leading_electron_vertex_distance',
        'reco_muon_multiplicity',
        'reco_photon_multiplicity',
        'reco_pion_multiplicity',
        'reco_proton_multiplicity',
        'reco_electron_multiplicity',
        'true_leading_electron_energy',
        'true_ele_beam_cosTheta',
        'true_leading_electron_start_x',
        'true_leading_electron_start_y',
        'true_leading_electron_start_z',
        'true_leading_electron_start_dir_x',
        'true_leading_electron_start_dir_y',
        'true_leading_electron_start_dir_z',
        'true_muon_multiplicity',
        'true_photon_multiplicity',
        'true_pion_multiplicity',
        'true_proton_multiplicity',
        'true_electron_multiplicity',
        'true_momentum_transfer',
        'true_hadronic_invariant_mass',

        # 'reco_flash_time_mod_cutflow',
        # 'true_neutrino_energy_cutflow',
        # 'reco_leading_electron_energy_cutflow',
        # 'reco_ele_beam_cosTheta_cutflow',
    ]

    outputDir = f'/nashome/m/micarrig/icarus/nuESpine/plots/{plot_category}/'
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    for var in plotVars:
        makePlot(var, outputDir)
