import matplotlib.pyplot as plt
import os
import sys
import pickle
os.environ['MEDULLA_PLOT_DIR'] =  '/Users/micarrig/Desktop/medulla/spineplot/'
sys.path.append('/Users/micarrig/Desktop/medulla/spineplot/')
from efficiency import SpineEfficiency
from analysis import Analysis

# Make the plot based on the provided variable, and save in dir
def makePlot(var, dir):

    fig = ana.run_interactively(var)
    fig.set_size_inches(10,8)
    plt.legend(
        facecolor="white",   # background color
        #edgecolor="black",   # optional: legend border
        framealpha=0.2       # 0 = fully transparent, 1 = fully opaque
    )
    plt.savefig(f'{dir}/plot_{var}.png', dpi=300, bbox_inches='tight')
    with open(f'{dir}/plot_{var}.pkl', 'wb') as f:
        pickle.dump(fig, f)
    plt.close()

if __name__ == "__main__":

    # ana = Analysis('/nashome/m/micarrig/icarus/nuESpine/nue_plots.toml', 'root://fndcadoor.fnal.gov:/icarus/scratch/users/micarrig/nueCCInclusiveSpine12/merged_output.root')
    ana = Analysis('/nashome/m/micarrig/icarus/nuESpine/medulla/analysis/systematics.toml', '/exp/icarus/data/users/micarrig/merged_v27p2_fixed.root')

    plot_category = 'systematics_total'


    # Each of these now has a signal-only and background-only variant
    # (see analysis/systematics.toml: SpineSystematics artists restricted
    # via `groups = [0]` / `groups = [1, 2, 3, 4, 5]`), producing separate
    # figures/plots for each rather than summing signal and background
    # into a single systematics band.
    plotVars = [
        'sys_eleE_signal',
        'sys_eleE_background',
        'sys_eleCosThetaNUMI_signal',
        'sys_eleCosThetaNUMI_background',
        'sys_nuE_signal',
        'sys_nuE_background',
        'sys_trueEleE_signal',
        'sys_trueEleE_background',
        'sys_trueEleCosThetaNUMI_signal',
        'sys_trueEleCosThetaNUMI_background',
    ]

    outputDir = f'/nashome/m/micarrig/icarus/nuESpine/plots_sept15/{plot_category}/'
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    for var in plotVars:
        makePlot(var, outputDir)
