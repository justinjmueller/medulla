import numpy as np
import pandas as pd
from tqdm import tqdm
import sys
import uproot
from ROOT import TFile, TEfficiency, TH1D, TGraphAsymmErrors, RDataFrame, TCanvas
import ROOT
from array import array
import ROOT, json
from math import nan


file_name = sys.argv[1]
if not file_name.endswith('.root'):
    print("Please provide a valid root file")
    sys.exit(1)

numi_dir = sys.argv[2] if len(sys.argv) > 2 else 'events/NuMIFull'
explicit_trees = sys.argv[3:] if len(sys.argv) > 3 else None

horn_current = 'fhc'
file_nu = uproot.open(file_name)
file_flux = uproot.open('/nashome/m/micarrig/icarus/nuESpine/flux.root')

flux_g4numi    = file_flux[f'g4numi_reweight_v03_01-->v03_02;1/{horn_current};1']
flux_beam_focus = file_flux[f'beam_focusing_uncertainties;1/{horn_current};1']
flux_pca        = file_flux['pca;1/principal_components;1']

ppfx_hweights_numu    = file_flux[f'ppfx_flux_weights/hweights_{horn_current}_numu;1']
ppfx_hweights_numubar = file_flux[f'ppfx_flux_weights/hweights_{horn_current}_numubar;1']
ppfx_hweights_nue     = file_flux[f'ppfx_flux_weights/hweights_{horn_current}_nue;1']
ppfx_hweights_nuebar  = file_flux[f'ppfx_flux_weights/hweights_{horn_current}_nuebar;1']

if explicit_trees:
    tree_names = explicit_trees
else:
    dir_obj = file_nu[numi_dir]
    tree_names = [k.split(';')[0] for k in dir_obj.keys(filter_classname='TTree')]

print(f"Processing directory '{numi_dir}' with trees: {tree_names}")

sigma = np.array([-3, -2, -1, 0, 1, 2, 3])
abs_sigma = np.abs(sigma)

N_PCA = 15


def compute_weights(nu_df):
    hysyst_beam_horn_2kA = []
    hysyst_beam_horn1_x_3mm = []
    hysyst_beam_horn1_y_3mm = []
    hysyst_beam_spot_1_3mm = []
    hysyst_beam_spot_1_7mm = []
    hysyst_beam_horn2_x_3mm = []
    hysyst_beam_horn2_y_3mm = []
    hysyst_beam_horns_0mm_water = []
    hysyst_beam_horns_2mm_water = []
    hysyst_beam_Beam_shift_x_1mm = []
    hysyst_beam_Beam_shift_y_1mm = []
    hysyst_beam_Target_z_7mm = []
    hpc = [[] for _ in range(N_PCA)]

    hysyst_beam_horn_2kA_sigma = []
    hysyst_beam_horn1_x_3mm_sigma = []
    hysyst_beam_horn1_y_3mm_sigma = []
    hysyst_beam_spot_1_3mm_sigma = []
    hysyst_beam_spot_1_7mm_sigma = []
    hysyst_beam_horn2_x_3mm_sigma = []
    hysyst_beam_horn2_y_3mm_sigma = []
    hysyst_beam_horns_0mm_water_sigma = []
    hysyst_beam_horns_2mm_water_sigma = []
    hysyst_beam_Beam_shift_x_1mm_sigma = []
    hysyst_beam_Beam_shift_y_1mm_sigma = []
    hysyst_beam_Target_z_7mm_sigma = []
    hpc_sigma = [[] for _ in range(N_PCA)]

    hnom_k0l_weights = []
    hnom_kpm_weights = []
    hnom_pipm_weights = []
    hnom_mu_weights = []
    run = []
    events = []
    subrun = []
    ppfx_cv_weight = []

    for e, event in tqdm(nu_df.iterrows()):
        run.append(event['Run'])
        subrun.append(event['Subrun'])
        events.append(event['Evt'])
        pdg = event['true_pdg']
        parent_pdg = event['true_parent_pdg']

        parent_pdg = np.nan_to_num(parent_pdg, nan=-999)
        pdg = np.nan_to_num(pdg, nan=-999)

        nu_e = event['true_neutrino_energy']

        if abs(int(parent_pdg)) == 311:  # K0
            if int(pdg) == 12:
                hnom_k0l_weights.append(flux_g4numi['hnom_nue_k0l_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_nue_k0l_weights;1'].axes[0].edges(), nu_e) - 1])
            elif int(pdg) == -12:
                hnom_k0l_weights.append(flux_g4numi['hnom_nuebar_k0l_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_nuebar_k0l_weights;1'].axes[0].edges(), nu_e) - 1])
            elif int(pdg) == 14:
                hnom_k0l_weights.append(flux_g4numi['hnom_numu_k0l_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_numu_k0l_weights;1'].axes[0].edges(), nu_e) - 1])
            elif int(pdg) == -14:
                hnom_k0l_weights.append(flux_g4numi['hnom_numubar_k0l_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_numubar_k0l_weights;1'].axes[0].edges(), nu_e) - 1])
            else:
                hnom_k0l_weights.append(flux_g4numi['hnom_k0l_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_k0l_weights;1'].axes[0].edges(), nu_e) - 1])
            hnom_mu_weights.append(float("nan"))
            hnom_pipm_weights.append(float("nan"))
            hnom_kpm_weights.append(float("nan"))
        elif abs(int(parent_pdg)) == 321:  # Kpm
            if int(pdg) == 12:
                hnom_kpm_weights.append(flux_g4numi['hnom_nue_kpm_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_nue_kpm_weights;1'].axes[0].edges(), nu_e) - 1])
            elif int(pdg) == -12:
                hnom_kpm_weights.append(flux_g4numi['hnom_nuebar_kpm_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_nuebar_kpm_weights;1'].axes[0].edges(), nu_e) - 1])
            elif int(pdg) == 14:
                hnom_kpm_weights.append(flux_g4numi['hnom_numu_kpm_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_numu_kpm_weights;1'].axes[0].edges(), nu_e) - 1])
            elif int(pdg) == -14:
                hnom_kpm_weights.append(flux_g4numi['hnom_numubar_kpm_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_numubar_kpm_weights;1'].axes[0].edges(), nu_e) - 1])
            else:
                hnom_kpm_weights.append(flux_g4numi['hnom_kpm_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_kpm_weights;1'].axes[0].edges(), nu_e) - 1])
            hnom_mu_weights.append(float("nan"))
            hnom_pipm_weights.append(float("nan"))
            hnom_k0l_weights.append(float("nan"))
        elif abs(int(parent_pdg)) == 211:  # pipm
            if int(pdg) == 12:
                hnom_pipm_weights.append(flux_g4numi['hnom_nue_pipm_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_nue_pipm_weights;1'].axes[0].edges(), nu_e) - 1])
            elif int(pdg) == -12:
                hnom_pipm_weights.append(flux_g4numi['hnom_nuebar_pipm_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_nuebar_pipm_weights;1'].axes[0].edges(), nu_e) - 1])
            elif int(pdg) == 14:
                hnom_pipm_weights.append(flux_g4numi['hnom_numu_pipm_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_numu_pipm_weights;1'].axes[0].edges(), nu_e) - 1])
            elif int(pdg) == -14:
                hnom_pipm_weights.append(flux_g4numi['hnom_numubar_pipm_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_numubar_pipm_weights;1'].axes[0].edges(), nu_e) - 1])
            else:
                hnom_pipm_weights.append(flux_g4numi['hnom_pipm_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_pipm_weights;1'].axes[0].edges(), nu_e) - 1])
            hnom_mu_weights.append(float("nan"))
            hnom_kpm_weights.append(float("nan"))
            hnom_k0l_weights.append(float("nan"))
        elif abs(int(parent_pdg)) == 13:  # mu
            if int(pdg) == 12:
                hnom_mu_weights.append(flux_g4numi['hnom_nue_mu_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_nue_mu_weights;1'].axes[0].edges(), nu_e) - 1])
            elif int(pdg) == -12:
                hnom_mu_weights.append(flux_g4numi['hnom_nuebar_mu_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_nuebar_mu_weights;1'].axes[0].edges(), nu_e) - 1])
            elif int(pdg) == 14:
                hnom_mu_weights.append(flux_g4numi['hnom_numu_mu_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_numu_mu_weights;1'].axes[0].edges(), nu_e) - 1])
            elif int(pdg) == -14:
                hnom_mu_weights.append(flux_g4numi['hnom_numubar_mu_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_numubar_mu_weights;1'].axes[0].edges(), nu_e) - 1])
            else:
                hnom_mu_weights.append(flux_g4numi['hnom_mu_weights;1'].values()[np.searchsorted(flux_g4numi['hnom_mu_weights;1'].axes[0].edges(), nu_e) - 1])
            hnom_pipm_weights.append(float("nan"))
            hnom_kpm_weights.append(float("nan"))
            hnom_k0l_weights.append(float("nan"))
        else:
            hnom_mu_weights.append(float("nan"))
            hnom_pipm_weights.append(float("nan"))
            hnom_kpm_weights.append(float("nan"))
            hnom_k0l_weights.append(float("nan"))

        pdg_to_nu_str = {12: 'nue', -12: 'nuebar', 14: 'numu', -14: 'numubar'}
        nu_str = pdg_to_nu_str.get(int(pdg))

        if nu_str is not None:
            def _bf(key):
                h = flux_beam_focus[f'{key}_{horn_current}_{nu_str};1']
                return h.values()[np.searchsorted(h.axes[0].edges(), nu_e) - 1]

            def _pca(idx):
                h = flux_pca[f'hpc_{idx}_{horn_current}_{nu_str};1']
                return h.values()[np.searchsorted(h.axes[0].edges(), nu_e) - 1]

            hysyst_beam_horn_2kA.append(list(map(lambda x: x + 1, abs_sigma[0:3] * _bf('hsyst_beam_Horn_m2kA'))))
            hysyst_beam_horn_2kA[-1].append(1.0)
            hysyst_beam_horn_2kA[-1].extend(list(map(lambda x: x + 1, abs_sigma[4:] * _bf('hsyst_beam_Horn_p2kA'))))
            hysyst_beam_horn_2kA_sigma.append(sigma)

            hysyst_beam_horn1_x_3mm.append(list(map(lambda x: x + 1, abs_sigma[0:3] * _bf('hsyst_beam_Horn1_x_m3mm'))))
            hysyst_beam_horn1_x_3mm[-1].append(1.0)
            hysyst_beam_horn1_x_3mm[-1].extend(list(map(lambda x: x + 1, abs_sigma[4:] * _bf('hsyst_beam_Horn1_x_p3mm'))))
            hysyst_beam_horn1_x_3mm_sigma.append(sigma)

            hysyst_beam_horn1_y_3mm.append(list(map(lambda x: x + 1, abs_sigma[0:3] * _bf('hsyst_beam_Horn1_y_m3mm'))))
            hysyst_beam_horn1_y_3mm[-1].append(1.0)
            hysyst_beam_horn1_y_3mm[-1].extend(list(map(lambda x: x + 1, abs_sigma[4:] * _bf('hsyst_beam_Horn1_y_p3mm'))))
            hysyst_beam_horn1_y_3mm_sigma.append(sigma)

            hysyst_beam_spot_1_3mm.append(list(map(lambda x: x + 1, sigma * _bf('hsyst_beam_Beam_spot_1_3mm'))))
            hysyst_beam_spot_1_3mm_sigma.append(sigma)

            hysyst_beam_spot_1_7mm.append(list(map(lambda x: x + 1, sigma * _bf('hsyst_beam_Beam_spot_1_7mm'))))
            hysyst_beam_spot_1_7mm_sigma.append(sigma)

            hysyst_beam_horn2_x_3mm.append(list(map(lambda x: x + 1, abs_sigma[0:3] * _bf('hsyst_beam_Horn2_x_m3mm'))))
            hysyst_beam_horn2_x_3mm[-1].append(1.0)
            hysyst_beam_horn2_x_3mm[-1].extend(list(map(lambda x: x + 1, abs_sigma[4:] * _bf('hsyst_beam_Horn2_x_p3mm'))))
            hysyst_beam_horn2_x_3mm_sigma.append(sigma)

            hysyst_beam_horn2_y_3mm.append(list(map(lambda x: x + 1, abs_sigma[0:3] * _bf('hsyst_beam_Horn2_y_m3mm'))))
            hysyst_beam_horn2_y_3mm[-1].append(1.0)
            hysyst_beam_horn2_y_3mm[-1].extend(list(map(lambda x: x + 1, abs_sigma[4:] * _bf('hsyst_beam_Horn2_y_p3mm'))))
            hysyst_beam_horn2_y_3mm_sigma.append(sigma)

            hysyst_beam_horns_0mm_water.append(list(map(lambda x: x + 1, sigma * _bf('hsyst_beam_Horns_0mm_water'))))
            hysyst_beam_horns_0mm_water_sigma.append(sigma)

            hysyst_beam_horns_2mm_water.append(list(map(lambda x: x + 1, sigma * _bf('hsyst_beam_Horns_2mm_water'))))
            hysyst_beam_horns_2mm_water_sigma.append(sigma)

            hysyst_beam_Beam_shift_x_1mm.append(list(map(lambda x: x + 1, abs_sigma[0:3] * _bf('hsyst_beam_Beam_shift_x_m1mm'))))
            hysyst_beam_Beam_shift_x_1mm[-1].append(1.0)
            hysyst_beam_Beam_shift_x_1mm[-1].extend(list(map(lambda x: x + 1, abs_sigma[4:] * _bf('hsyst_beam_Beam_shift_x_p1mm'))))
            hysyst_beam_Beam_shift_x_1mm_sigma.append(sigma)

            hysyst_beam_Beam_shift_y_1mm.append(list(map(lambda x: x + 1, abs_sigma[0:3] * _bf('hsyst_beam_Beam_shift_y_m1mm'))))
            hysyst_beam_Beam_shift_y_1mm[-1].append(1.0)
            hysyst_beam_Beam_shift_y_1mm[-1].extend(list(map(lambda x: x + 1, abs_sigma[4:] * _bf('hsyst_beam_Beam_shift_y_p1mm'))))
            hysyst_beam_Beam_shift_y_1mm_sigma.append(sigma)

            hysyst_beam_Target_z_7mm.append(list(map(lambda x: x + 1, abs_sigma[0:3] * _bf('hsyst_beam_Target_z_m7mm'))))
            hysyst_beam_Target_z_7mm[-1].append(1.0)
            hysyst_beam_Target_z_7mm[-1].extend(list(map(lambda x: x + 1, abs_sigma[4:] * _bf('hsyst_beam_Target_z_p7mm'))))
            hysyst_beam_Target_z_7mm_sigma.append(sigma)

            for i in range(N_PCA):
                hpc[i].append(list(map(lambda x: x + 1, sigma * _pca(i))))
                hpc_sigma[i].append(sigma)
        else:
            print(f"Warning: PDG {pdg} not recognized for PPFX correction.")

        if int(pdg) == 14:
            ppfx_cv_weight.append(ppfx_hweights_numu.values()[np.searchsorted(ppfx_hweights_numu.axes[0].edges(), nu_e) - 1])
        elif int(pdg) == -14:
            ppfx_cv_weight.append(ppfx_hweights_numubar.values()[np.searchsorted(ppfx_hweights_numubar.axes[0].edges(), nu_e) - 1])
        elif int(pdg) == 12:
            ppfx_cv_weight.append(ppfx_hweights_nue.values()[np.searchsorted(ppfx_hweights_nue.axes[0].edges(), nu_e) - 1])
        elif int(pdg) == -12:
            ppfx_cv_weight.append(ppfx_hweights_nuebar.values()[np.searchsorted(ppfx_hweights_nuebar.axes[0].edges(), nu_e) - 1])
        else:
            ppfx_cv_weight.append(float('nan'))

    return {
        "run": run, "events": events, "subrun": subrun,
        "hysyst_beam_horn_2kA": hysyst_beam_horn_2kA,
        "hysyst_beam_horn1_x_3mm": hysyst_beam_horn1_x_3mm,
        "hysyst_beam_horn1_y_3mm": hysyst_beam_horn1_y_3mm,
        "hysyst_beam_spot_1_3mm": hysyst_beam_spot_1_3mm,
        "hysyst_beam_spot_1_7mm": hysyst_beam_spot_1_7mm,
        "hysyst_beam_horn2_x_3mm": hysyst_beam_horn2_x_3mm,
        "hysyst_beam_horn2_y_3mm": hysyst_beam_horn2_y_3mm,
        "hysyst_beam_horns_0mm_water": hysyst_beam_horns_0mm_water,
        "hysyst_beam_horns_2mm_water": hysyst_beam_horns_2mm_water,
        "hysyst_beam_Beam_shift_x_1mm": hysyst_beam_Beam_shift_x_1mm,
        "hysyst_beam_Beam_shift_y_1mm": hysyst_beam_Beam_shift_y_1mm,
        "hysyst_beam_Target_z_7mm": hysyst_beam_Target_z_7mm,
        **{f"hpc_{i}": hpc[i] for i in range(N_PCA)},
        "hysyst_beam_horn_2kA_sigma": hysyst_beam_horn_2kA_sigma,
        "hysyst_beam_horn1_x_3mm_sigma": hysyst_beam_horn1_x_3mm_sigma,
        "hysyst_beam_horn1_y_3mm_sigma": hysyst_beam_horn1_y_3mm_sigma,
        "hysyst_beam_spot_1_3mm_sigma": hysyst_beam_spot_1_3mm_sigma,
        "hysyst_beam_spot_1_7mm_sigma": hysyst_beam_spot_1_7mm_sigma,
        "hysyst_beam_horn2_x_3mm_sigma": hysyst_beam_horn2_x_3mm_sigma,
        "hysyst_beam_horn2_y_3mm_sigma": hysyst_beam_horn2_y_3mm_sigma,
        "hysyst_beam_horns_0mm_water_sigma": hysyst_beam_horns_0mm_water_sigma,
        "hysyst_beam_horns_2mm_water_sigma": hysyst_beam_horns_2mm_water_sigma,
        "hysyst_beam_Beam_shift_x_1mm_sigma": hysyst_beam_Beam_shift_x_1mm_sigma,
        "hysyst_beam_Beam_shift_y_1mm_sigma": hysyst_beam_Beam_shift_y_1mm_sigma,
        "hysyst_beam_Target_z_7mm_sigma": hysyst_beam_Target_z_7mm_sigma,
        **{f"hpc_{i}_sigma": hpc_sigma[i] for i in range(N_PCA)},
        "ppfx_cv_weight": ppfx_cv_weight,
    }


def ensure_dir(rootdir, path):
    cur = rootdir
    for part in path.strip("/").split("/"):
        d = cur.GetDirectory(part) or cur.mkdir(part)
        cur = d
    return cur


CTOR_VEC = {
    "f": lambda: ROOT.std.vector('float')(),
    "d": lambda: ROOT.std.vector('double')(),
    "i": lambda: ROOT.std.vector('int')(),
    "l": lambda: ROOT.std.vector('long long')(),
}
CAST = {"f": float, "d": float, "i": int, "l": int}
LEAF = {"f": "F", "d": "D", "i": "I", "l": "L"}
INT_SENTINEL = -9999


def write_flux_tree(root_file, directory, treeName, weights):
    spec = {
        "Run":    ("i", weights["run"]),
        "Subrun": ("i", weights["subrun"]),
        "Evt":    ("i", weights["events"]),
        "hysyst_beam_horn_2kA_v":        ("f", weights["hysyst_beam_horn_2kA"]),
        "hysyst_beam_horn1_x_3mm_v":     ("f", weights["hysyst_beam_horn1_x_3mm"]),
        "hysyst_beam_horn1_y_3mm_v":     ("f", weights["hysyst_beam_horn1_y_3mm"]),
        "hysyst_beam_spot_1_3mm_v":      ("f", weights["hysyst_beam_spot_1_3mm"]),
        "hysyst_beam_spot_1_7mm_v":      ("f", weights["hysyst_beam_spot_1_7mm"]),
        "hysyst_beam_horn2_x_3mm_v":     ("f", weights["hysyst_beam_horn2_x_3mm"]),
        "hysyst_beam_horn2_y_3mm_v":     ("f", weights["hysyst_beam_horn2_y_3mm"]),
        "hysyst_beam_horns_0mm_water_v": ("f", weights["hysyst_beam_horns_0mm_water"]),
        "hysyst_beam_horns_2mm_water_v": ("f", weights["hysyst_beam_horns_2mm_water"]),
        "hysyst_beam_Beam_shift_x_1mm_v": ("f", weights["hysyst_beam_Beam_shift_x_1mm"]),
        "hysyst_beam_Beam_shift_y_1mm_v": ("f", weights["hysyst_beam_Beam_shift_y_1mm"]),
        "hysyst_beam_Target_z_7mm_v":    ("f", weights["hysyst_beam_Target_z_7mm"]),
        **{f"hysyst_hpc_{i}_v": ("f", weights[f"hpc_{i}"]) for i in range(N_PCA)},
        "hysyst_beam_horn_2kA_sigma":        ("f", weights["hysyst_beam_horn_2kA_sigma"]),
        "hysyst_beam_horn1_x_3mm_sigma":     ("f", weights["hysyst_beam_horn1_x_3mm_sigma"]),
        "hysyst_beam_horn1_y_3mm_sigma":     ("f", weights["hysyst_beam_horn1_y_3mm_sigma"]),
        "hysyst_beam_spot_1_3mm_sigma":      ("f", weights["hysyst_beam_spot_1_3mm_sigma"]),
        "hysyst_beam_spot_1_7mm_sigma":      ("f", weights["hysyst_beam_spot_1_7mm_sigma"]),
        "hysyst_beam_horn2_x_3mm_sigma":     ("f", weights["hysyst_beam_horn2_x_3mm_sigma"]),
        "hysyst_beam_horn2_y_3mm_sigma":     ("f", weights["hysyst_beam_horn2_y_3mm_sigma"]),
        "hysyst_beam_horns_0mm_water_sigma": ("f", weights["hysyst_beam_horns_0mm_water_sigma"]),
        "hysyst_beam_horns_2mm_water_sigma": ("f", weights["hysyst_beam_horns_2mm_water_sigma"]),
        "hysyst_beam_Beam_shift_x_1mm_sigma": ("f", weights["hysyst_beam_Beam_shift_x_1mm_sigma"]),
        "hysyst_beam_Beam_shift_y_1mm_sigma": ("f", weights["hysyst_beam_Beam_shift_y_1mm_sigma"]),
        "hysyst_beam_Target_z_7mm_sigma":    ("f", weights["hysyst_beam_Target_z_7mm_sigma"]),
        **{f"hysyst_hpc_{i}_sigma": ("f", weights[f"hpc_{i}_sigma"]) for i in range(N_PCA)},
    }

    def is_2d_of_7(seq):
        if not seq or not hasattr(seq, '__iter__'):
            return False
        first = seq[0]
        if not hasattr(first, '__iter__'):
            return False
        return all((row is None) or (hasattr(row, '__iter__') and len(row) == 7) for row in seq)

    keys_2d7 = [k for k, (_, v) in spec.items() if is_2d_of_7(v)]
    keys_1d  = [k for k, (_, v) in spec.items() if not is_2d_of_7(v)]

    if keys_2d7:
        N = len(spec[keys_2d7[0]][1])
    elif keys_1d:
        N = len(spec[keys_1d[0]][1])
    else:
        raise ValueError("spec is empty")

    for k in keys_2d7:
        if len(spec[k][1]) != N:
            raise ValueError(f"{k}: outer length {len(spec[k][1])} != {N}")
    for k in keys_1d:
        if len(spec[k][1]) != N:
            raise ValueError(f"{k}: length {len(spec[k][1])} != {N}")

    tdir = ensure_dir(root_file, directory)
    tdir.cd()

    t = ROOT.TTree(f"{treeName}_NuMIfluxsimTree", "per-entry vectors (len=7) plus scalars")
    branch_order = list(spec.keys())
    tdir.WriteObject(ROOT.TObjString(json.dumps(branch_order)), f"branch_labels_json_{treeName}")

    vec_br = {}
    for name in keys_2d7:
        code, rows = spec[name]
        vec = CTOR_VEC[code]()
        t.Branch(name, vec)
        vec_br[name] = (vec, code, rows)

    sca_br = {}
    for name in keys_1d:
        code, seq = spec[name]
        buf = array('f' if code in ('f', 'd') else 'i', [0.0] if code in ('f', 'd') else [0])
        t.Branch(name, buf, f"{name}/{LEAF[code]}")
        sca_br[name] = (buf, code, seq)

    for i in range(N):
        for name, (vec, code, rows) in vec_br.items():
            vec.clear()
            row = rows[i]
            if row is None:
                if code in ('f', 'd'):
                    for _ in range(7): vec.push_back(nan)
                else:
                    for _ in range(7): vec.push_back(INT_SENTINEL)
            else:
                for x in row:
                    if x is None:
                        vec.push_back(nan if code in ('f', 'd') else INT_SENTINEL)
                    else:
                        vec.push_back(CAST[code](x))
        for name, (buf, code, seq) in sca_br.items():
            x = seq[i]
            buf[0] = (nan if code in ('f', 'd') else INT_SENTINEL) if x is None else CAST[code](x)
        t.Fill()

    tdir.WriteTObject(t, t.GetName(), "Overwrite")

    main_tree = root_file.Get(f"{directory}/{treeName}")
    existing_branch = main_tree.GetBranch("ppfx_cv_weight")
    if existing_branch:
        main_tree.GetListOfBranches().Remove(existing_branch)
    ppfx_buf = array('d', [0.0])
    ppfx_branch = main_tree.Branch("ppfx_cv_weight", ppfx_buf, "ppfx_cv_weight/D")
    ppfx_cv_weight = weights["ppfx_cv_weight"]
    for i in range(main_tree.GetEntries()):
        ppfx_buf[0] = float(ppfx_cv_weight[i]) if not np.isnan(ppfx_cv_weight[i]) else 1.0
        ppfx_branch.Fill()
    root_file.cd(directory)
    main_tree.Write("", ROOT.TObject.kOverwrite)


f = ROOT.TFile(file_name, "UPDATE")

for treeName in tree_names:
    print(f"\n--- Processing tree: {treeName} ---")
    nu_df = file_nu[f'{numi_dir}/{treeName};1'].arrays(library='pd')
    weights = compute_weights(nu_df)
    write_flux_tree(f, numi_dir, treeName, weights)

f.Close()
