# imports
import uproot
from ROOT import TFile, TEfficiency, TH1D, TGraphAsymmErrors, RDataFrame, TCanvas
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('/exp/icarus/app/users/lkashur/medulla_dev/srcs/medulla/spineplot/styles/pi0ana.mplstyle')
import pandas as pd
import toml
import re
import argparse
from array import array

def main(args):
    
    ######################
    ### Load configuration
    ######################
    config = toml.load(args.config)

    ######################
    ### Load input
    ######################
    rf = uproot.open(args.input)
    all_samples = []
    for k,v in config['samples'].items():
        
        sample_str = v['key']
        sample_tree = v['trees'][0]
        
        sample_dict = {}
        sample_dict['sample'] = k
        sample_dict['tree'] = rf[f'events/{sample_str}/{sample_tree}']
        sample_dict['df'] = sample_dict['tree'].arrays(library='pd')
        sample_dict['pot'] = rf[f'events/{sample_str}/POT'].to_numpy()[0][0]
        sample_dict['livetime'] = rf[f'events/{sample_str}/Livetime'].to_numpy()[0][0]
        all_samples.append(sample_dict)
        
    # Add signal, which typically isn't included in TOML plotting config
    # first check for events/mc/signal, followed by events/cv/signal
    if any('events/cv/signal' in s for s in rf.keys()):
        sig_loc = 'events/cv'
    if any('events/mc/signal' in s for s in rf.keys()):
        sig_loc = 'events/mc'
    all_samples.append({'sample' : 'signal',
                        'tree' : rf[f'{sig_loc}/signal'],
                        'df' : rf[f'{sig_loc}/signal'].arrays(library='pd'),
                        'pot' : rf[f'{sig_loc}/POT'].to_numpy()[0][0],
                        'livetime' : rf[f'{sig_loc}/Livetime'].to_numpy()[0][0]})

    # Add purity info, which typicaly isn't includedin TOML plotting config
    # To-do

    ######################
    ### Tasks
    ######################
    
    # Efficiency calculation
    calculate_efficiency(all_samples)
    
    # Efficiency calculation (by cut)
    cuts = {'Flash Cut': ['reco_flash_icarus_satisfied == 1'],
            'Fiducial Cut' : ['reco_flash_icarus_satisfied == 1', 
                              'reco_fiducial_satisfied == 1'],
            'Base Topology Cut' : ['reco_flash_icarus_satisfied == 1', 
                                   'reco_fiducial_satisfied == 1', 
                                   'reco_base_topology_satisfied == 1'],
            'Leading Shower Energy Cut' : ['reco_flash_icarus_satisfied == 1', 
                                           'reco_fiducial_satisfied == 1', 
                                           'reco_base_topology_satisfied == 1', 
                                           'reco_leading_shower_energy_satisfied == 1'],
            'pi0 Mass Cut' : ['reco_flash_icarus_satisfied == 1', 
                              'reco_fiducial_satisfied == 1', 
                              'reco_base_topology_satisfied == 1', 
                              'reco_leading_shower_energy_satisfied == 1', 
                              'reco_valid_pi0_mass_satisfied == 1']}
    calculate_efficiency_by_cut(all_samples, cuts, true_category = config['analysis']['category_branch'])
    
    # Purity calculation
    calculate_purity(all_samples, true_category = config['analysis']['category_branch'])

    # Purity calculation (by cut, "purity" tree required)
    #calculate_purity_by_cut(all_samples, cuts, true_category = config['analysis']['category_branch'])

    # Efficiency claculation (as a funtion of true variable)
    plot_eff_by_var(all_samples, 'true_muon_momentum_mag', config, ana_eff_bins_equal=True, nbins_eff=40)


def add_plot_labels(ax, pot, adj_y=0.025, title=str()):
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.025, #0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = xrange[0]
    prelim_label = r'ICARUS Work-in-Progress'
    ax.text(x=usex, y=usey, s=prelim_label, fontsize=14, color='#d67a11')
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.02, 0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = (xrange[1] + xrange[0]) / 2
    ax.text(x=usex, y=usey, s=f'{title}', fontsize=14, fontweight='bold', color='black')
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.02, 0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = xrange[1]
    mag = int(np.floor(np.log10(pot)))
    usepot = pot/10**mag
    s = f'{usepot:.2f}'+f'$\\times 10^{{{mag}}}$ POT'
    ax.text(x=usex, y=usey, s=s, fontsize=13, color='black', horizontalalignment='right')


'''
Function to calculate selection effiency

'''
def calculate_efficiency(all_samples):
    
    # Retrieve the signal dataframe
    sig_df = [d for d in all_samples if d['sample'] == 'signal'][0]['df']
    
    # Numerator
    total_signal_events = len(sig_df)

    # Denominator
    selected_signal_events = len(sig_df[sig_df['reco_all_cuts_satisfied'] == 1])
    
    eff =  selected_signal_events / total_signal_events
    print(f'{eff:.4f}')
    return eff

'''
Function to calculate selection effiency for each cut, sequentially
'''
def calculate_efficiency_by_cut(all_samples, cuts, true_category = 'true_category'):

    # Retrieve signal dataframe
    sig_df = [d for d in all_samples if d['sample'] == 'signal'][0]['df']

    print(f'No Cut Efficiency: 1.000')
    for key, value in cuts.items():
        cond_string = str()
        for i,c in enumerate(value):
            cond_string += c
            if i != len(value) - 1: cond_string += ' & '

        sel_by_cut = sig_df.query(cond_string)
        eff = len(sel_by_cut) / len(sig_df)
        #eff = len(sel_by_cut[sel_by_cut[true_category] == 0]) / len(sig_df)
        print(f'{key} Efficiency: {eff:.4f}')

'''
Function to calculate selection purity
'''
def calculate_purity(all_samples, true_category='true_category'):
    
    # Retrieve selected dataframes
    sel_nu_df = [d for d in all_samples if d['sample'] == 'mc_nu'][0]['df']
    sel_cos_df = [d for d in all_samples if d['sample'] == 'mc_cos'][0]['df']
    sel_mc_df = pd.concat([sel_nu_df, sel_cos_df])
    sel_offbeam_df = [d for d in all_samples if d['sample'] == 'offbeam'][0]['df']

    # Retrieve scaling
    livetime_onbeam = [d for d in all_samples if d['sample'] == 'onbeam'][0]['livetime']
    livetime_offbeam = [d for d in all_samples if d['sample'] == 'offbeam'][0]['livetime']
    pot_onbeam = [d for d in all_samples if d['sample'] == 'onbeam'][0]['pot']
    pot_mc = [d for d in all_samples if d['sample'] == 'mc_nu'][0]['pot']

    # Numerator
    matched_selected_events = (pot_onbeam/pot_mc)*len(sel_mc_df[sel_mc_df[true_category] == 0])

    # Denominator
    total_selected_events = (pot_onbeam/pot_mc)*(len(sel_mc_df)) + (livetime_onbeam/livetime_offbeam)*len(sel_offbeam_df)

    pur = matched_selected_events / total_selected_events
    print(f'{pur:.4f}')
    return pur

'''
Function to calculate selection purity for each cut, sequentially
'''
def calculate_purity_by_cut(all_samples, cuts, true_category = 'true_category'):

    # Retrieve purity dataframes
    sel_nu_pur_df = [d for d in all_samples if d['sample'] == 'mc_nu_pur'][0]['df']
    sel_cos_pur_df = [d for d in all_samples if d['sample'] == 'mc_cos_pur'][0]['df']
    sel_mc_pur_df = pd.concat([sel_nu_pur_df, sel_cos_pur_df])
    sel_offbeam_pur_df = [d for d in all_samples if d['sample'] == 'offbeam_pur'][0]['df']

    # Retrieve scaling
    livetime_onbeam = [d for d in all_samples if d['sample'] == 'onbeam'][0]['livetime']
    livetime_offbeam = [d for d in all_samples if d['sample'] == 'offbeam'][0]['livetime']
    pot_onbeam = [d for d in all_samples if d['sample'] == 'onbeam'][0]['pot']
    pot_mc = [d for d in all_samples if d['sample'] == 'mc_nu'][0]['pot']

    # No cut purity
    no_cut_pur = len(sel_mc_pur_df[sel_mc_pur_df[true_category] == 0]) / (len(sel_mc_pur_df) + ((livetime_onbeam/livetime_offbeam)*len(sel_offbeam_pur_df)))
    print(f'No Cut Purity: {no_cut_pur:.4f}')

    # Loop over cuts
    for key, value in cuts.items():
        cond_string = str()
        for i,c in enumerate(value):
            cond_string += c
            if i != len(value) - 1: cond_string += ' & '

        sel_by_cut_mc = sel_mc_pur_df.query(cond_string)
        sel_by_cut_offbeam = sel_offbeam_pur_df.query(cond_string)
        sel_by_cut = pd.concat([sel_by_cut_mc, sel_by_cut_offbeam])
 
        cut_pur = len(sel_by_cut_mc[sel_by_cut_mc[true_category] == 0]) / (len(sel_by_cut_mc) + (livetime_onbeam/livetime_offbeam)*len(sel_by_cut_offbeam))
        print(f'{key} Purity: {cut_pur:.4f}')


'''
Function to plot efficiency (purity) as a function of a true (reconstucted) variable
'''
def plot_eff_by_var(all_samples, var, config, ana_eff_bins_equal=True, nbins_eff=25):
    bcs = []
    bxerr0s = []
    bxerr1s = []
    byerr0s = []
    byerr1s = []
    vals = []

    # Retrieve signal dataframe
    sig_df = [d for d in all_samples if d['sample'] == 'signal'][0]['df']
    
    # Retrieve POT
    pot = [d for d in all_samples if d['sample'] == 'signal'][0]['pot']

    # Filtering (keep in range of analyis bins)
    sig_df = sig_df[(sig_df[var] > config['variables'][var]['range'][0]) & (sig_df[var] < config['variables'][var]['range'][-1])]

    # Apply analyis binning
    ana_bins_df = pd.DataFrame(config['variables'][var]['range'], columns=['bin']) # prepare overflow
    sig_df['var_q_ana'] = pd.cut(sig_df[var], ana_bins_df['bin'])

    # Apply eff binning
    _nbins_eff = 0
    if ana_eff_bins_equal:
        _nbins_eff = config['variables'][var]['nbins']
        sig_df['var_q_eff'] = pd.cut(sig_df[var], np.linspace(config['variables'][var]['range'][0], config['variables'][var]['range'][-1], _nbins_eff+1))
    else:
        _nbins_eff = nbins_eff
        sig_df['var_q_eff'] = pd.cut(sig_df[var], np.linspace(config['variables'][var]['range'][0], config['variables'][var]['range'][-1], _nbins_eff+1))
    
    #######
    ### ana
    #######
    hedges_ana = sorted([i.left for i in sig_df.var_q_ana.unique().tolist()])
    lastedge_ana = max([i.right for i in sig_df.var_q_ana.unique().tolist()])
    hedges_ana.append(lastedge_ana)
    hedges_ana = np.array(hedges_ana)
    counts = []
    for i,(name,group) in enumerate(sig_df.groupby('var_q_ana')):
        counts.append(len(group))
    
    ##########
    ### eff
    ##########
    hedges_eff = sorted([i.left for i in sig_df.var_q_eff.unique().tolist()])
    lastedge_eff = max([i.right for i in sig_df.var_q_eff.unique().tolist()])
    hedges_eff.append(lastedge_eff)
    hedges_eff = np.array(hedges_eff)

    hpass = TH1D('hpass', '', _nbins_eff, hedges_eff)
    htotal = TH1D('htotal', '', _nbins_eff, hedges_eff)

    for i,(name,group) in enumerate(sig_df.groupby('var_q_eff')):

        bxerr0s.append(name.left)
        bcs.append(name.mid)
        bxerr1s.append(name.right)
        vals.append(len(group) / len(group[group['reco_all_cuts_satisfied'] == 1]))
        hpass.SetBinContent(i+1, len(group[group['reco_all_cuts_satisfied'] == 1]))
        htotal.SetBinContent(i+1, len(group))

    gr = TGraphAsymmErrors()
    gr.Divide(hpass, htotal, 'cl=0.683 b(1,1) mode')                                             
    for i in range(_nbins_eff):
        byerr0s.append(gr.GetErrorYlow(i))
        byerr1s.append(gr.GetErrorYhigh(i))

    fig, ax1 = plt.subplots(figsize=(8,8))
    ax1.bar(x=hedges_ana[:-1], height=counts, width=np.diff(hedges_ana), align='edge', fc='C1', alpha=0.5, ec='none')
    ax1.set_xlabel(config['variables'][var]['xlabel'])
    ax1.set_ylabel('Entries')
    ax2 = ax1.twinx()
    ax2.set_ylim([0,1])
    ax2.errorbar(bcs, vals, xerr=[np.array(bcs) - np.array(bxerr0s), np.array(bxerr1s) - np.array(bcs)], yerr=[np.array(byerr0s), np.array(byerr1s)], fmt='o', capsize=2)
    ylabel = 'Efficiency'
    ax2.set_ylabel(ylabel)
    ax2.spines['right'].set_color('C0')
    ax2.yaxis.label.set_color('C0')
    ax2.tick_params(axis='y', colors='C0')
    #plt.title(r'Signal $\nu_{\mu}$ CC $\pi^{0}$', fontsize=16)
    plottitle=r'Signal $\nu_{\mu}$ CC $\pi^{0}$'
    add_plot_labels(ax1,pot, adj_y=0.030, title=str())
    plt.savefig('eff_vs_' + var + '.png')
    #plt.savefig('eff_vs_' + var + '.pdf', format='pdf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    main(args)
