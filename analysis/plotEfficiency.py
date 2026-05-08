import sys, os
import uproot
import ROOT as r
import numpy as np
import array as arr
import awkward as ak
import matplotlib.pyplot as plt
from memory_profiler import profile
import pickle

sys.path.append('/nashome/m/micarrig/sbnd/nueCCNp/medulla/spineplot/')

from efficiency import SpineEfficiency
from analysis import Analysis

r.gStyle.SetOptStat(0)

def main():
    # Load the analysis
    #ana = Analysis('/nashome/m/micarrig/icarus/nuESpine/efficiency_analysis.toml', 'root://fndcadoor.fnal.gov:/icarus/scratch/users/micarrig/nueCCInclusiveSpine12/merged_output.root')
    #ana = Analysis('/nashome/m/micarrig/icarus/nuESpine/efficiency_recoInt.toml', 'root://fndcadoor.fnal.gov:/icarus/scratch/users/micarrig/nueCCInclusiveSpine12/merged_output.root')
    ana = Analysis('/nashome/m/micarrig/sbnd/nueCCNp/medulla/analysis/efficiency_nueCCNp.toml', 'root://fndcadoor.fnal.gov:/sbnd/scratch/users/micarrig/nueCCNp_v2/output/merged_output.root')

    outputDir = 'plots/selection_efficiency_trueTree'

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    # #Run interactively for specific figure
    # table = ana.run_interactively('selection_efficiency_table')
    # plt.savefig(f'{outputDir}/selection_efficiency_table.png', dpi=300)
    # with open(f'{outputDir}/selection_efficiency_table.pkl', 'wb') as f:
    #     pickle.dump(table, f)

    # eff_nuE = ana.run_interactively('efficiency_vs_energy')
    # plt.savefig(f'{outputDir}/efficiency_vs_energy.png', dpi=300)
    # with open(f'{outputDir}/efficiency_vs_energy.pkl', 'wb') as f:
    #     pickle.dump(eff_nuE, f)

    # eff_eleE = ana.run_interactively('efficiency_vs_electronEnergy')
    # plt.savefig(f'{outputDir}/efficiency_vs_electronEnergy.png', dpi=300)
    # with open (f'{outputDir}/efficiency_vs_electronEnergy.pkl', 'wb') as f:
    #     pickle.dump(eff_eleE, f)

    # eff_eleAngle = ana.run_interactively('efficiency_vs_electronAngle')
    # plt.savefig(f'{outputDir}/efficiency_vs_electronAngle.png', dpi=300)
    # with open(f'{outputDir}/efficiency_vs_electronAngle.pkl', 'wb') as f:
    #     pickle.dump(eff_eleAngle, f)

    # eff_pTrans = ana.run_interactively('efficiency_vs_trueMomentumTransfer')
    # plt.savefig(f'{outputDir}/efficiency_vs_trueMomentumTransfer.png', dpi=300)
    # with open(f'{outputDir}/efficiency_vs_trueMomentumTransfer.pkl', 'wb') as f:
    #     pickle.dump(eff_pTrans, f)

    # eff_invMass = ana.run_interactively('efficiency_vs_trueInvariantMass')
    # plt.savefig(f'{outputDir}/efficiency_vs_trueInvMass.png', dpi=300)
    # with open(f'{outputDir}/efficiency_vs_trueInvMass.pkl', 'wb') as f:
    #     pickle.dump(eff_invMass, f)

    # eff_dpT = ana.run_interactively('efficiency_vs_recoDPT')
    # plt.savefig(f'{outputDir}/efficiency_vs_recoDPT.png', dpi=300)
    # with open(f'{outputDir}/efficiency_vs_recoDPT.pkl', 'wb') as f:
    #     pickle.dump(eff_dpT, f)

    # eff_dalphaT = ana.run_interactively('efficiency_vs_recoDalphaT')
    # plt.savefig(f'{outputDir}/efficiency_vs_recoDalphaT.png', dpi=300)
    # with open(f'{outputDir}/efficiency_vs_recoDalphaT.pkl', 'wb') as f:
    #     pickle.dump(eff_dalphaT, f)

    eff_elePolar = ana.run_interactively('efficiency_vs_recoElectronPolarAngle')
    plt.savefig(f'{outputDir}/efficiency_vs_recoElectronPolarAngle.png', dpi=300)
    with open(f'{outputDir}/efficiency_vs_recoElectronPolarAngle.pkl', 'wb') as f:
        pickle.dump(eff_elePolar, f)

    eff_eleOpening = ana.run_interactively('efficiency_vs_recoElectronOpeningAngle')
    plt.savefig(f'{outputDir}/efficiency_vs_recoElectronOpeningAngle.png', dpi=300)
    with open(f'{outputDir}/efficiency_vs_recoElectronOpeningAngle.pkl', 'wb') as f:
        pickle.dump(eff_eleOpening, f)

    # eff_flashTime = ana.run_interactively('efficiency_vs_recoFlashTime')
    # plt.savefig(f'{outputDir}/efficiency_vs_recoFlashTime.png', dpi=300)
    # with open(f'{outputDir}/efficiency_vs_recoFlashTime.pkl', 'wb') as f:
    #     pickle.dump(eff_flashTime, f)


if __name__ == "__main__":

    main()