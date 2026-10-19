import sys, os
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
    # ana = Analysis('/exp/icarus/app/users/micarrig/nuESpine/medulla/analysis/efficiency_analysis.toml', '/exp/icarus/data/users/micarrig/merged_v30.root')
    # ana = Analysis('/nashome/m/micarrig/icarus/nuESpine/medulla/analysis/efficiency_withError.toml', '/exp/icarus/data/users/micarrig/merged_v30.root')
    ana = Analysis('/exp/icarus/app/users/micarrig/nuESpine/medulla/analysis/efficiency_recoInt.toml', '/exp/icarus/data/users/micarrig/merged_v30.root')

    output_dir = '/exp/icarus/app/users/micarrig/nuESpine/plots_sept15/efficiency_neutrino_type_trueTree/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    table = ana.run_interactively('selection_efficiency_table')
    plt.savefig(f'{output_dir}/selection_efficiency_table.png', dpi=300)
    with open(f'{output_dir}/selection_efficiency_table.pkl', 'wb') as f:
        pickle.dump(table, f)

    eff_nuE = ana.run_interactively('efficiency_vs_energy')
    plt.savefig(f'{output_dir}/efficiency_vs_energy.png', dpi=300)
    with open(f'{output_dir}/efficiency_vs_energy.pkl', 'wb') as f:
        pickle.dump(eff_nuE, f)

    eff_eleE = ana.run_interactively('efficiency_vs_electronEnergy')
    plt.savefig(f'{output_dir}/efficiency_vs_electronEnergy.png', dpi=300)
    with open (f'{output_dir}/efficiency_vs_electronEnergy.pkl', 'wb') as f:
        pickle.dump(eff_eleE, f)

    eff_eleAngle = ana.run_interactively('efficiency_vs_electronAngle')
    plt.savefig(f'{output_dir}/efficiency_vs_electronAngle.png', dpi=300)
    with open(f'{output_dir}/efficiency_vs_electronAngle.pkl', 'wb') as f:
        pickle.dump(eff_eleAngle, f)

    eff_pTrans = ana.run_interactively('efficiency_vs_trueMomentumTransfer')
    plt.savefig(f'{output_dir}/efficiency_vs_trueMomentumTransfer.png', dpi=300)
    with open(f'{output_dir}/efficiency_vs_trueMomentumTransfer.pkl', 'wb') as f:
        pickle.dump(eff_pTrans, f)

    eff_invMass = ana.run_interactively('efficiency_vs_trueInvariantMass')
    plt.savefig(f'{output_dir}/efficiency_vs_trueInvMass.png', dpi=300)
    with open(f'{output_dir}/efficiency_vs_trueInvMass.pkl', 'wb') as f:
        pickle.dump(eff_invMass, f)

    eff_tEleEle = ana.run_interactively('efficiency_vs_trueLeadingElectronEnergy')
    plt.savefig(f'{output_dir}/efficiency_vs_trueLeadingElectronEnergy.png', dpi=300)
    with open(f'{output_dir}/efficiency_vs_trueLeadingElectronEnergy.pkl', 'wb') as f:
        pickle.dump(eff_tEleEle, f) 

    eff_tEleAngle = ana.run_interactively('efficiency_vs_trueLeadingElectronCosTheta')
    plt.savefig(f'{output_dir}/efficiency_vs_trueLeadingElectronAngle.png', dpi=300)
    with open(f'{output_dir}/efficiency_vs_trueLeadingElectronAngle.pkl', 'wb') as f:
        pickle.dump(eff_tEleAngle, f)   

    # eff_flashTime = ana.run_interactively('efficiency_vs_recoFlashTime')
    # plt.savefig(f'{output_dir}/efficiency_vs_recoFlashTime.png', dpi=300)
    # with open(f'{output_dir}/efficiency_vs_recoFlashTime.pkl', 'wb') as f:
    #     pickle.dump(eff_flashTime, f)

    # eff_nu = ana.run_interactively('efficiency_vs_trueNu')
    # plt.savefig(f'{output_dir}/efficiency_vs_trueNu.png', dpi=300)
    # with open(f'{output_dir}/efficiency_vs_trueNu.pkl', 'wb') as f:
    #     pickle.dump(eff_nu, f)


    # eff_eleE_total = ana.run_interactively('efficiency_vs_trueLeadingElectronEnergy_total')
    # plt.savefig(f'{output_dir}/efficiency_vs_trueLeadingElectronEnergy_total.png', dpi=300)
    # with open(f'{output_dir}/efficiency_vs_trueLeadingElectronEnergy_total.pkl', 'wb') as f:
    #     pickle.dump(eff_eleE_total, f)

    # eff_eleE_flux = ana.run_interactively('efficiency_vs_trueLeadingElectronEnergy_flux')
    # plt.savefig(f'{output_dir}/efficiency_vs_trueLeadingElectronEnergy_flux.png', dpi=300)
    # with open(f'{output_dir}/efficiency_vs_trueLeadingElectronEnergy_flux.pkl', 'wb') as f:
    #     pickle.dump(eff_eleE_flux, f)

    # eff_eleE_detector = ana.run_interactively('efficiency_vs_trueLeadingElectronEnergy_detector')
    # plt.savefig(f'{output_dir}/efficiency_vs_trueLeadingElectronEnergy_detector.png', dpi=300)
    # with open(f'{output_dir}/efficiency_vs_trueLeadingElectronEnergy_detector.pkl', 'wb') as f:
    #     pickle.dump(eff_eleE_detector, f)

    # eff_eleE_xsec = ana.run_interactively('efficiency_vs_trueLeadingElectronEnergy_xsec')
    # plt.savefig(f'{output_dir}/efficiency_vs_trueLeadingElectronEnergy_xsec.png', dpi=300)
    # with open(f'{output_dir}/efficiency_vs_trueLeadingElectronEnergy_xsec.pkl', 'wb') as f:
    #     pickle.dump(eff_eleE_xsec, f)

    # eff_eleE_totalStat = ana.run_interactively('efficiency_vs_trueLeadingElectronEnergy_totalStat')
    # plt.savefig(f'{output_dir}/efficiency_vs_trueLeadingElectronEnergy_totalStat.png', dpi=300)
    # with open(f'{output_dir}/efficiency_vs_trueLeadingElectronEnergy_totalStat.pkl', 'wb') as f:
    #     pickle.dump(eff_eleE_totalStat, f)


if __name__ == "__main__":

    main()