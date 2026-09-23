import matplotlib.pyplot as plt
import os
import sys
import pickle
os.environ['MEDULLA_PLOT_DIR'] =  '/nashome/m/micarrig/icarus/nuESpine/medulla/spineplot/'
sys.path.append('/nashome/m/micarrig/icarus/nuESpine/medulla/spineplot/')
from analysis import Analysis

def makePlot(var, dir):
    fig = ana.run_interactively(var)
    fig.set_size_inches(10, 8)
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

    ana = Analysis('/nashome/m/micarrig/icarus/nuESpine/medulla/analysis/nminus1_plots.toml', '/exp/icarus/data/users/micarrig/merged_v27p2_fixed.root')

    plotVars = [
        'nminus1_fiducial_x',
        'nminus1_fiducial_y',
        'nminus1_fiducial_z',
        'nminus1_flash_cut',
        'nminus1_containment_x',
        'nminus1_containment_y',
        'nminus1_containment_z',
        'nminus1_particle_multiplicity',
        'nminus1_electron_softmax',
        'nminus1_primary_softmax',
        'nminus1_no_muons',
        'nminus1_vertex_distance',
        'nminus1_electron_energy',
    ]

    outputDir = '/nashome/m/micarrig/icarus/nuESpine/plots_sept15/nminus1/'
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    for var in plotVars:
        makePlot(var, outputDir)
