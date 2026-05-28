import sys
import uproot
import ROOT as r
import numpy as np
import array as arr
import awkward as ak
import matplotlib.pyplot as plt
from memory_profiler import profile
import pickle

sys.path.append('/nashome/m/micarrig/icarus/nuESpine/medulla/spineplot/')

from efficiency import SpineEfficiency
from analysis import Analysis

r.gStyle.SetOptStat(0)

def main():
    # Load the analysis
    #ana = Analysis('/nashome/m/micarrig/icarus/nuESpine/efficiency_analysis.toml', 'root://fndcadoor.fnal.gov:/icarus/scratch/users/micarrig/nueCCInclusiveSpine12/merged_output.root')
    #ana = Analysis('/nashome/m/micarrig/icarus/nuESpine/efficiency_recoInt.toml', 'root://fndcadoor.fnal.gov:/icarus/scratch/users/micarrig/nueCCInclusiveSpine12/merged_output.root')
    ana = Analysis('/nashome/m/micarrig/icarus/nuESpine/efficiency_analysis.toml', '/nashome/m/micarrig/icarus/nuESpine/output_nueCCInclusive2.root')

    # Run interactively for specific figure
    # table = ana.run_interactively('selection_efficiency_table')
    # plt.savefig('plots/selection_efficiency_table.png', dpi=300)
    # with open('plots/selection_efficiency_table.pkl', 'wb') as f:
    #     pickle.dump(table, f)

    # eff_nuE = ana.run_interactively('efficiency_vs_energy')
    # plt.savefig('plots/efficiency_vs_energy.png', dpi=300)
    # with open('plots/efficiency_vs_energy.pkl', 'wb') as f:
    #     pickle.dump(eff_nuE, f)

    # eff_eleE = ana.run_interactively('efficiency_vs_electronEnergy')
    # plt.savefig('plots/efficiency_vs_electronEnergy.png', dpi=300)
    # with open ('plots/efficiency_vs_electronEnergy.pkl', 'wb') as f:
    #     pickle.dump(eff_eleE, f)

    # eff_eleAngle = ana.run_interactively('efficiency_vs_electronAngle')
    # plt.savefig('plots/efficiency_vs_electronAngle.png', dpi=300)
    # with open('plots/efficiency_vs_electronAngle.pkl', 'wb') as f:
    #     pickle.dump(eff_eleAngle, f)

    # eff_pTrans = ana.run_interactively('efficiency_vs_trueMomentumTransfer')
    # plt.savefig('plots/efficiency_vs_trueMomentumTransfer.png', dpi=300)
    # with open('plots/efficiency_vs_trueMomentumTransfer.pkl', 'wb') as f:
    #     pickle.dump(eff_pTrans, f)

    # eff_invMass = ana.run_interactively('efficiency_vs_trueInvariantMass')
    # plt.savefig('plots/efficiency_vs_trueInvMass.png', dpi=300)
    # with open('plots/efficiency_vs_trueInvMass.pkl', 'wb') as f:
    #     pickle.dump(eff_invMass, f)

    eff_flashTime = ana.run_interactively('efficiency_vs_recoFlashTime')
    plt.savefig('plots/efficiency_vs_recoFlashTime.png', dpi=300)
    with open('plots/efficiency_vs_recoFlashTime.pkl', 'wb') as f:
        pickle.dump(eff_flashTime, f)


if __name__ == "__main__":

    main()