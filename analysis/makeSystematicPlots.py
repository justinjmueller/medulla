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
    ana = Analysis('/nashome/m/micarrig/icarus/nuESpine/medulla/analysis/systematics.toml', '/nashome/m/micarrig/icarus/nuESpine/gundam_nueCCInclusive.root')

    plot_category = 'systematics_total'


    plotVars = [
        'sys_eleE',
        'sys_eleCosThetaNUMI',
        'sys_nuE',
        'sys_trueEleE',
        'sys_trueEleCosThetaNUMI'
    ]

    outputDir = f'/nashome/m/micarrig/icarus/nuESpine/plots/{plot_category}/'
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    for var in plotVars:
        makePlot(var, outputDir)