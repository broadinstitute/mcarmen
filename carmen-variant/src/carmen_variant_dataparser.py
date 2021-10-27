#!/usr/bin/env python
#
# Copied largely from the Jupyter notebook 'CFL048_Biomark_csv_v13ab.ipynb'
# Run like:
#  ```script_name.py --fcount ADAPTF001 --chip-dim 96
#   --in-csv data/20200728_1362564078_rawdata.csv --in-layout
#   data/ADAPTF001_96_assignment.xlsx --timepoints 25 --out-dir output/```

import pandas as pd
import numpy as np
from os import listdir,path
import warnings
import math
import csv
import argparse
import gzip
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
import csv
import os
import os
warnings.filterwarnings('ignore')

# Read args
parser = argparse.ArgumentParser()
parser.add_argument('--fcount',
                    required=True,
                    help="Experiment count (e.g., 'CFL091')")
parser.add_argument('--chip-dim',
                    required=True,
                    choices=['96', '192'],
                    help=("Fluidigm chip type (number of samples/targets)"))
parser.add_argument('--timepoints',
                    type=int,
                    default=37,
                    help="Number of timepoints")
parser.add_argument('--in-csv',
                    required=True,
                    help=("CSV file with data"))
parser.add_argument('--in-layout',
                    required=True,
                    help=("Excel sheet giving layout"))
parser.add_argument('--out-dir',
                    required=True,
                    help=("Output directory"))
args = parser.parse_args()


# Set up
ifc = args.chip_dim # 96 or 192
instrument_type = 'BM' # EP1 or BM (Biomark)
fcount = args.fcount
tgap = 1 # time gap between mixing of reagents (end of chip loading) and t0 image in minutes

if instrument_type == 'EP1':
    count_tp = 1 # End point run
else:
    count_tp = args.timepoints # number of timepoints, standard for 2h techdev is 25 atm

# Define variables based on inputs
exp_name = fcount + '_' + ifc
#csv_file = fcount + '_' + ifc + '_' + instrument_type +'_rawdata.csv' # csv file with raw data from Fluidigm RT PCR software
csv_file = args.in_csv

# Write a path to that excel sheet
layout_file = args.in_layout

out_folder = args.out_dir

if path.exists(layout_file) == False:
    raise Exception((layout_file + " doesn't exist"))

if path.exists(csv_file) == False:
    raise Exception((csv_file + " doesn't exist"))

# Definition of functions
# Create a dictonary for timepoints
time_assign = {}
for cycle in range(1,38):
    tpoint = "t" + str(cycle)
    time_assign[tpoint] = tgap + 3 + (cycle-1) * 5
# calculate the real timing of image
# used for image and axis labeling
def gettime(tname):
    realt = time_assign[tname]
    return (realt)

if instrument_type == 'BM':
    if ifc == '96':
        probe_df = pd.read_csv(csv_file,header = 18449, nrows = 9216) # FAM
        reference_df = pd.read_csv(csv_file, header = 9231, nrows = 9216) # ROX
        bkgd_ref_df = pd.read_csv(csv_file, header = 27667, nrows = 9216)
        bkgd_probe_df = pd.read_csv(csv_file,header = 36885, nrows = 9216)
    if ifc == '192':
        probe_df = pd.read_csv(csv_file,header = 9233, nrows = 4608)
        reference_df = pd.read_csv(csv_file, header = 4623, nrows = 4608)
        bkgd_ref_df = pd.read_csv(csv_file, header = 13843, nrows = 4608)
        bkgd_probe_df = pd.read_csv(csv_file,header = 18453, nrows = 4608)
elif instrument_type == 'EP1':
    if ifc == '192':
        probe_df = pd.read_csv(csv_file,header = 9238, nrows = 4608)
        reference_df = pd.read_csv(csv_file, header = 4628, nrows = 4608)
        bkgd_ref_df = pd.read_csv(csv_file, header = 18458, nrows = 4608)
        bkgd_probe_df = pd.read_csv(csv_file,header = 23068, nrows = 4608)

# Get rid of stuff
c_to_drop = 'Unnamed: ' + str(count_tp+1)
probe_df = probe_df.set_index("Chamber ID").drop(c_to_drop,1)
reference_df = reference_df.set_index("Chamber ID").drop(c_to_drop,1)
bkgd_ref_df = bkgd_ref_df.set_index("Chamber ID").drop(c_to_drop,1)
bkgd_probe_df = bkgd_probe_df.set_index("Chamber ID").drop(c_to_drop,1)

probe_df.columns = probe_df.columns.str.lstrip() # remove spaces from beginning of column names
reference_df.columns = reference_df.columns.str.lstrip()
bkgd_ref_df.columns = bkgd_ref_df.columns.str.lstrip()
bkgd_probe_df.columns = bkgd_probe_df.columns.str.lstrip()

# rename column names
probe_df.columns = ['t' + str(col) for col in probe_df.columns]
reference_df.columns = ['t' + str(col) for col in reference_df.columns]
bkgd_ref_df.columns = ['t' + str(col) for col in bkgd_ref_df.columns]
bkgd_probe_df.columns = ['t' + str(col) for col in bkgd_probe_df.columns]

# if an error like "Passed header=9233 but only 4618 lines in file" comes up, you probably exported the wrong csv file from the RT-PCR software (you need 'table results with raw data')
# if an error like ""['Unnamed: 26'] not found in axis" "comes up, look at the probe_df and look at the last column. Likely, the number of timepoints is wrong

# Substract the background from the probe and reference data
probe_bkgd_substracted = probe_df.subtract(bkgd_probe_df)
ref_bkgd_substracted = reference_df.subtract(bkgd_ref_df)

# Normalize the probe signal with the reference dye signal
signal_df = pd.DataFrame(probe_bkgd_substracted/ref_bkgd_substracted)

# reset index
signal_df = signal_df.reset_index()
# split Column ID into SampleID and AssayID
splitassignment = signal_df['Chamber ID'].str.split("-",n=1,expand=True)
signal_df["sampleID"] = splitassignment[0]
signal_df["assayID"] = splitassignment[1]

#set index again to Chamber ID
signal_df = signal_df.set_index('Chamber ID')

sampleID_list = signal_df.sampleID.unique()
assayID_list = signal_df.assayID.unique()

# Save csv
signal_df.to_csv(path.join(out_folder, exp_name+ '_' +instrument_type + '_1_signal_bkgdsubtracted_norm_' + str(count_tp) +'.csv'))

# Create two dictionaries that align the IFC wells to the sample and assay names
samples_layout = pd.read_excel(path.join('',layout_file),sheet_name='layout_samples').applymap(str)
assays_layout = pd.read_excel(path.join('',layout_file),sheet_name='layout_assays').applymap(str)

# Create a dictionary with assay/sample numbers and their actual crRNA / target name
assays = pd.read_excel(path.join('',layout_file),sheet_name='assays')
assay_dict = dict(zip(assays_layout.values.reshape(-1),assays.values.reshape(-1)))

samples = pd.read_excel(path.join('',layout_file),sheet_name='samples')
samples_dict = dict(zip(samples_layout.values.reshape(-1),samples.values.reshape(-1)))

# Map assay and sample names
signal_df['assay'] = signal_df['assayID'].map(assay_dict)
signal_df['sample'] = signal_df['sampleID'].map(samples_dict)

# Save csv
signal_df.to_csv(path.join(out_folder, exp_name+'_' +instrument_type +'_2_signal_bkgdsubtracted_norm_named_' + str(count_tp) +'.csv'))

# # Transform and summarize data for plotting

#create list with timepoints
# count_tp = 23 # in case you were wrong before
t_names = []
for i in range(1,count_tp+1):
    t_names.append(('t' + str(i)))

# Create a list of all assays and samples
# only indicate columns with unique assays. np.unique could be used on the list, but messes up our prefered the order
if ifc == '96':
    new_array = np.stack(assays[['C1','C2','C3','C4','C5']].values,axis=-1)
if ifc == '192':
    new_array = np.stack(assays[['C1','C2', 'C3']].values,axis=-1)

assay_list = np.concatenate(new_array).tolist()
print('identified crRNAs: ',len(assay_list))

# Do this for the samples
if ifc == '96':
    new_array = np.stack(samples[['C1','C2','C3','C4','C5','C6','C7', 'C8','C9','C10','C11','C12']].values,axis=-1)
if ifc == '192':
    new_array = np.stack(samples[['C1','C2','C3','C4','C5','C6', 'C7','C8','C9','C10','C11','C12', 'C13','C14','C15','C16','C17','C18', 'C19','C20','C21','C22','C23','C24']].values,axis=-1)
#     new_array = np.stack(samples[['C1','C2','C3','C4','C5','C6',\
#                                   'C7','C8','C9','C10','C11','C12']].values,axis=-1)

sample_list = np.concatenate(new_array).tolist()
print('identified samples: ',len(sample_list))

# Grouped medians
medians = signal_df.groupby(['assay','sample']).median()

# Creating dataframes in a loop:
# https://stackoverflow.com/questions/30635145/create-multiple-dataframes-in-loop
med_frames = {}
for name in t_names:
    time_med = signal_df.groupby(['assay','sample']).median()[name].unstack()
    time_med.index.names=['']
    time_med.columns.names=['']
    med_frames[name] = time_med


# Write results

# Write TSV, with one row per (timepoint, sample (target), assay (guide)) triplet
# The value written is: {median across replicates[(probe signal - probe background) / (reference signal - reference background)]}
# for different time points
with gzip.open(path.join(out_folder, exp_name+ '_'+ instrument_type + '_merged.tsv.gz'), 'wt') as fw:
    def write_row(row):
        fw.write('\t'.join(str(x) for x in row) + '\n')
    header = ['timepoint', 'minute', 'guide', 'target', 'value']
    write_row(header)

    for tp in med_frames.keys():
        rt = gettime(tp)
        for target in sample_list:
            for guide in assay_list:
                v = med_frames[tp][target][guide]
                row = [tp, rt, guide, target, v]
                write_row(row)


# Plot heatmap
def plt_heatmap(df_dict, samplelist, assaylist, tp):
    frame = df_dict[tp][samplelist].reindex(assaylist)
    fig, axes = plt.subplots(1,1,figsize=(len(frame.columns.values)*0.5,len(frame.index.values)*0.5))
    ax = sns.heatmap(frame,cmap='Reds',square = True,cbar_kws={'pad':0.002}, annot_kws={"size": 20})
    rt = gettime(tp)
    plt.title(exp_name+' '+str(rt)+'min - median values', size = 28)
    plt.xlabel('Samples', size = 14)
    plt.ylabel('Assays', size = 14)
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom + 0.5, top - 0.5)
    ax.tick_params(axis="y", labelsize=16)
    ax.tick_params(axis="x", labelsize=16)
    plt.yticks(rotation=0)

    tgt_num = len(sample_list)
    gd_num = len(assay_list)
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom + 0.5, top - 0.5)
    h_lines = np.arange(3,gd_num,3)
    v_lines = np.arange(3,tgt_num,3)
    axes.hlines(h_lines, colors = 'silver',alpha=0.9,linewidths = 0.35,*axes.get_xlim())
    axes.vlines(v_lines, colors = 'silver',alpha=0.9,linewidths = 0.35,*axes.get_ylim())

    plt.savefig(path.join(out_folder, exp_name + '_'+instrument_type +'2_heatmap_'+str(tp)+'.png'), format='png', bbox_inches='tight', dpi=400)
# plt_heatmap(med_frames, sample_list, assay_list, 't7')
# plt_heatmap(med_frames, sample_list, assay_list, 't12')
# plt_heatmap(med_frames, sample_list, assay_list, 't18')
# plt_heatmap(med_frames, sample_list, assay_list, 't24')
# plt_heatmap(med_frames, sample_list, assay_list, 't30')
# plt_heatmap(med_frames, sample_list, assay_list, 't37')

# Create a dataframe for the time series
seen = set()
sample_list_ = [x for x in sample_list if x not in seen and not seen.add(x)]
seen = set()
assay_list_ = [x for x in assay_list if x not in seen and not seen.add(x)]
med_series = pd.DataFrame(0,\
                          columns = sample_list_,\
                          index = assay_list_).astype('object')
# Pull out signal for each guide, sample, and timepoint
for guide in assay_list_:
    for target in sample_list_:
        med_list = []
        for tp in t_names:
            # pull out the right combination for each timepoint
            curr_med = med_frames[tp].loc[guide,target]
            med_list.append(curr_med)
        # add into df as an array
#         print(guide,target,med_series.loc[guide,target],'shape',len(med_list),np.asarray(med_list).shape)
        med_series.loc[guide,target] = np.asarray(med_list)

# Constructs the control and variant guide pairs from assignment sheet tab
# tab "guide_map" which should be constructed manually
variant_gmap_ = pd.read_excel(path.join('',layout_file),sheet_name='guide_map')
        ## IGNORE - NEED MORE STANDARDIZED NAMING TO FULLY AUTOMATE
        # variant_gmap_['guide'] = variant_gmap_['guide_name'].apply(lambda x: x.split('_vs')[0])
        # variant_gmap_['type_'] = 'var'
        # variant_gmap_.loc[variant_gmap_.guide.str.contains('Ref'),'type_'] = 'ctrl'
        # guide_ctrls = variant_gmap_[variant_gmap_.type_=='ctrl']['guide'].values
        # guide_vars = variant_gmap_[variant_gmap_.type_=='var']['guide'].values
        # variant_gmap_.loc[(variant_gmap_.)]
variant_gmap = pd.DataFrame(columns=['control_guide','variant_guide'])
k = 0
for i in variant_gmap_.group.unique():
    subgroup = variant_gmap_[variant_gmap_.group==i]
    guide_ctrls = subgroup[subgroup.type=='ctrl']['guide_name'].values
    guide_vars = subgroup[subgroup.type=='var']['guide_name'].values
    for ctrl in guide_ctrls:
        for var in guide_vars:
            variant_gmap.loc[k,'control_guide'] = ctrl
            variant_gmap.loc[k,'variant_guide'] = var
            k+=1
variant_gmap.index = variant_gmap['control_guide']+' : '+variant_gmap['variant_guide']

# Mapping SNP crRNA guide pairs to variants using variant_SNPmap tab
variant_snpmap = pd.read_excel(path.join('',layout_file),sheet_name='variant_SNPmap')
variant_gmap = pd.concat([variant_gmap,pd.DataFrame(columns=variant_snpmap.columns.values)])
variant_gmap['var_SNP'] = variant_gmap['variant_guide'].apply(lambda x: x.split('_vs_')[0])

for var in variant_snpmap:
    snps = variant_snpmap[var].dropna().values
    for snp in snps:
        variant_gmap.loc[variant_gmap.var_SNP.str.contains(snp),var] = 1
variant_gmap = variant_gmap.fillna(0)
snpmap_dict = variant_snpmap.to_dict('list')

# CUSTOM MODIFICATIONS
# CHECK THAT YOU ARE HAPPY WITH THIS MAP!!!
crRNA_remove = ['WGAN_SpikeRefE484K_0.1uM : wgan.E484Q_vs_Ref.1427','wgan.Ref_vs_P681R.2013 : GEN_P681H_vs_Spike_Ref-0.1uM']
variant_gmap = variant_gmap[~variant_gmap.index.isin(crRNA_remove)]

# dictionary for seedstock names
seedstock_dict = {
'B.1.1.7' : 'Alpha',
'B.1.351' : 'Beta',
'P.1' : 'Gamma',
'B.1.617.2' : 'Delta',
'B.1.429' : 'Epsilon',
'B.1.427' : 'Epsilon',
'B.1.427/429' : 'Epsilon'}

# plotting ground truth variant snp info
variant_gmap_ = variant_gmap[['Alpha', 'Beta', 'Gamma',
       'Delta', 'Epsilon', 'Kappa', 'Eta', 'Lamda', 'var_SNP']]
variant_gmap_.loc['sum',:] = variant_gmap_.sum()
variant_gmap_.iloc[-1,-1] = 'sum'
variant_gmap_.set_index('var_SNP',inplace=True)
variant_gmap_.sort_values(by='sum',axis=1,ascending=False,inplace=True)
variant_gmap_.drop(index='sum',inplace=True)
variant_gmap_.rename(columns={'B.1.427/429':'Epsilon'},inplace=True)
sns.heatmap(variant_gmap_.T,linewidths=0.05,linecolor='silver',cmap='Blues',cbar=False)
plt.title('Ground Truth SNPs for Variants')
plt.savefig('./groundtruth_var_SNPb.png',bbox_inches='tight',dpi=300)
plt.savefig(path.join(out_folder,exp_name+'_'+instrument_type+'_groundtruth_var_SNP.png'), format='png', bbox_inches='tight', dpi=400)
plt.savefig(path.join(out_folder,exp_name+'_'+instrument_type+'_groundtruth_var_SNP.pdf'), format='pdf', bbox_inches='tight')

# Takes the control:variant guide ratio for the med_series df for each sample
# Populates ratio time series into new df called ctrlvar_series_ratio
series_T = med_series.T
# series_T = med_series.applymap(lambda x: np.diff(x)).T

ctrlvar_series_ratio = pd.DataFrame(index=series_T.index,columns=variant_gmap.index)
for i in variant_gmap.index:
    ctrlg = variant_gmap.loc[i,'control_guide']
    varg = variant_gmap.loc[i,'variant_guide']
    ctrlvar_series_ratio[i] = series_T[ctrlg]/series_T[varg]
ctrlvar_series_ratio = ctrlvar_series_ratio.T
ctrlvar_series_ratio.to_pickle(path.join(out_folder, exp_name+ '_' +instrument_type + '_ctrlvar_ratio_tpseries.pkl'))

varctrl_series_ratio = pd.DataFrame(index=series_T.index,columns=variant_gmap.index)
for i in variant_gmap.index:
    ctrlg = variant_gmap.loc[i,'control_guide']
    varg = variant_gmap.loc[i,'variant_guide']
    varctrl_series_ratio[i] = series_T[varg]/series_T[ctrlg]
varctrl_series_ratio = varctrl_series_ratio.T
varctrl_series_ratio.to_pickle(path.join(out_folder, exp_name+ '_' +instrument_type + '_varctrl_ratio_tpseries.pkl'))

# get the max fold change for control:variant ratio and earliest timepoint to reach that max fold change
ctrlvar_maxFC = pd.DataFrame(index=ctrlvar_series_ratio.index,columns=ctrlvar_series_ratio.columns)
ctrlvar_maxFCtp = pd.DataFrame(index=ctrlvar_series_ratio.index,columns=ctrlvar_series_ratio.columns)
for i in ctrlvar_series_ratio.columns[:]:
    for j in ctrlvar_series_ratio.index:
        vals = ctrlvar_series_ratio.loc[j,i]
        ctrlvar_maxFC.loc[j,i] = max(vals)
        max_tp = np.argwhere(vals==max(vals))
        ctrlvar_maxFCtp.loc[j,i] = min(max_tp)[0]
ctrlvar_maxFC = ctrlvar_maxFC.astype(float)
ctrlvar_maxFCtp = ctrlvar_maxFCtp.astype(float)
ctrlvar_maxFC.to_csv(path.join(out_folder, exp_name+ '_' +instrument_type + '_ctrlvar_maxFC.csv'))
ctrlvar_maxFCtp.to_csv(path.join(out_folder, exp_name+ '_' +instrument_type + '_ctrlvar_maxFCtp.csv'))

# get the max fold change for variant:control ratio and earliest timepoint to reach that max fold change
varctrl_maxFC = pd.DataFrame(index=varctrl_series_ratio.index,columns=varctrl_series_ratio.columns)
varctrl_maxFCtp = pd.DataFrame(index=varctrl_series_ratio.index,columns=varctrl_series_ratio.columns)
for i in varctrl_series_ratio.columns[:]:
    for j in varctrl_series_ratio.index:
        vals = varctrl_series_ratio.loc[j,i]
        varctrl_maxFC.loc[j,i] = max(vals)
        max_tp = np.argwhere(vals==max(vals))
        varctrl_maxFCtp.loc[j,i] = min(max_tp)[0]
varctrl_maxFC = varctrl_maxFC.astype(float)
varctrl_maxFCtp = varctrl_maxFCtp.astype(float)
varctrl_maxFC.to_csv(path.join(out_folder, exp_name+ '_' +instrument_type + '_varctrl_maxFC.csv'))
varctrl_maxFCtp.to_csv(path.join(out_folder, exp_name+ '_' +instrument_type + '_varctrl_maxFCtp.csv'))

# collect all maxFC adn their tp's, decide if var:ctrl or ctrl:var FC is higher
# note which is higher under maxFC_crRNA
# new columns with maxFC, its invert, and its tp
varctrl_maxFC_melt = varctrl_maxFC.reset_index().melt(id_vars=['index']).rename(columns={'index':'crRNA_pair','variable':'sample_name','value':'maxFC_varctrl'})
ctrlvar_maxFC_melt = ctrlvar_maxFC.reset_index().melt(id_vars=['index']).rename(columns={'index':'crRNA_pair','variable':'sample_name','value':'maxFC_ctrlvar'})

varctrl_maxFCtp_melt = varctrl_maxFCtp.reset_index().melt(id_vars=['index']).rename(columns={'index':'crRNA_pair','variable':'sample_name','value':'maxFC_varctrl_tp'})
ctrlvar_maxFCtp_melt = ctrlvar_maxFCtp.reset_index().melt(id_vars=['index']).rename(columns={'index':'crRNA_pair','variable':'sample_name','value':'maxFC_ctrlvar_tp'})

merged_maxFC = varctrl_maxFC_melt.merge(ctrlvar_maxFC_melt)
merged_maxFC = merged_maxFC.merge(varctrl_maxFCtp_melt)
merged_maxFC = merged_maxFC.merge(ctrlvar_maxFCtp_melt)
merged_maxFC['maxFC_crRNA'] = merged_maxFC[['maxFC_varctrl','maxFC_ctrlvar']].idxmax(axis=1)
merged_maxFC['maxFC'] = merged_maxFC.apply(lambda x: x[x.maxFC_crRNA],axis=1)
merged_maxFC['maxFC_invert'] = 1/merged_maxFC['maxFC']
merged_maxFC['maxFC_tp'] = merged_maxFC.apply(lambda x: x[x.maxFC_crRNA+'_tp'],axis=1)
merged_maxFC['WT'] = np.nan
merged_maxFC['Mut'] = np.nan
merged_maxFC.loc[(merged_maxFC.maxFC_crRNA=='maxFC_ctrlvar'),'WT'] = merged_maxFC.loc[(merged_maxFC.maxFC_crRNA=='maxFC_ctrlvar'),'maxFC']
merged_maxFC.loc[(merged_maxFC.maxFC_crRNA=='maxFC_ctrlvar'),'Mut'] = merged_maxFC.loc[(merged_maxFC.maxFC_crRNA=='maxFC_ctrlvar'),'maxFC_invert']
merged_maxFC.loc[(merged_maxFC.maxFC_crRNA=='maxFC_varctrl'),'WT'] = merged_maxFC.loc[(merged_maxFC.maxFC_crRNA=='maxFC_varctrl'),'maxFC_invert']
merged_maxFC.loc[(merged_maxFC.maxFC_crRNA=='maxFC_varctrl'),'Mut'] = merged_maxFC.loc[(merged_maxFC.maxFC_crRNA=='maxFC_varctrl'),'maxFC']
merged_maxFC.to_csv(path.join(out_folder,exp_name+'_'+instrument_type+'_crRNApair_maxFC_melted_merged.csv'))
