#!/usr/bin/env python
#
# Copied largely from the Jupyter notebook 'CFL048_Biomark_csv_v13c.ipynb'
# Run like:
#  ```script_name.py --fcounta CFL181 --fcountb CFL182
#   --FCthresh-file "FCthresholds_v4.csv" --in-csva "../data/CFL181/1691455285.csv"
#   --in-csvb "../data/CFL182/1691278294.csv"
#   --in-layouta "../data/CFL181/192_assignment.xlsx"
#   --in-layoutb "../data/CFL182/192_assignment.xlsx"
#   --out-dir "./mCARMEN/"```

# .py --fcounta CFL181 --fcountb CFL182 --FCthresh-file "FCthresholds_v4.csv" --in-csva "../data/CFL181/1691455285.csv" --in-csvb "../data/CFL182/1691278294.csv" --in-layouta "../data/CFL181/192_assignment.xlsx" --in-layoutb "../data/CFL182/192_assignment.xlsx" --out-dir "./mCARMEN"
import pandas as pd
import numpy as np
from os import listdir,path
import warnings
import math
import csv
import argparse
import gzip
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from tqdm import tqdm
import csv
import os
import glob
warnings.filterwarnings('ignore')
plt.rcParams['pdf.fonttype'] = 42

# takes data from full variant panel spanning two 192.24 IFCs
# therefore, requiring two datasets to be read in for variant calling analysis

# Read args
parser = argparse.ArgumentParser()
parser.add_argument('--fcounta',
                    required=True,
                    help="Experiment count (e.g., 'CFL091')")
parser.add_argument('--fcountb',
                    required=True,
                    help="Experiment count (e.g., 'CFL092')")
# parser.add_argument('--chip-dim',
#                     required=True,
#                     choices=['96', '192'],
#                     help=("Fluidigm chip type (number of samples/targets)"))
parser.add_argument('--FCthresh-file',
                    required=True,
                    help=("csv with threshold values"))
parser.add_argument('--timepoints',
                    type=int,
                    default=37,
                    help="Number of timepoints")
parser.add_argument('--in-csva',
                    required=True,
                    help=("CSV file with data"))
parser.add_argument('--in-csvb',
                    required=True,
                    help=("CSV file with data"))
parser.add_argument('--in-layouta',
                    required=True,
                    help=("Excel sheet giving layout"))
parser.add_argument('--in-layoutb',
                    required=True,
                    help=("Excel sheet giving layout"))
parser.add_argument('--out-dir',
                    required=True,
                    help=("Output directory, \
                    must contain outputs from variant_parse_data.py"))
args = parser.parse_args()


# Set up
# args for notebook a
a_fcount = args.fcounta
a_chip_dim = "192"#args.chip_dim
a_timepoints = args.timepoints
a_in_csv = args.in_csva
a_in_layout = args.in_layouta
a_out_dir = args.out_dir
# os.makedirs(a_out_dir,exist_ok=True)

# args for notebook b
b_fcount = args.fcountb
b_chip_dim = "192"#args.chip_dim
b_timepoints = args.timepoints
b_in_csv = args.in_csvb
b_in_layout = args.in_layoutb
b_out_dir = args.out_dir
# os.makedirs(b_out_dir,exist_ok=True)


# define combined dataset name, path and make out folder
# ab_expname = a_fcount+'_'+b_fcount
# os.makedirs('../data/'+ab_expname,exist_ok=True)
out_path = args.out_dir
timepoints = max([a_timepoints,b_timepoints])

# Set up a
a_ifc = a_chip_dim # 96 or 192
a_instrument_type = 'BM' # EP1 or BM (Biomark)
a_fcount = a_fcount
a_tgap = 1 # time gap between mixing of reagents (end of chip loading) and t0 image in minutes

if a_instrument_type == 'EP1':
    a_count_tp = 1 # End point run
else:
    a_count_tp = a_timepoints # number of timepoints, standard for 2h techdev is 25 atm

# Define variables based on inputs
a_exp_name = a_fcount + '_' + a_ifc
# csv_file = fcount + '_' + ifc + '_' + instrument_type +'_rawdata.csv' # csv file with raw data from Fluidigm RT PCR software
a_csv_file = a_in_csv

# Write a path to that excel sheet
a_layout_file = a_in_layout

a_out_folder = a_out_dir

if path.exists(a_layout_file) == False:
    raise Exception((a_layout_file + " doesn't exist"))

if path.exists(a_csv_file) == False:
    raise Exception((a_csv_file + " doesn't exist"))

# Set up b
b_ifc = b_chip_dim # 96 or 192
b_instrument_type = 'BM' # EP1 or BM (Biomark)
b_fcount = b_fcount
b_tgap = 1 # time gap between mixing of reagents (end of chip loading) and t0 image in minutes

if b_instrument_type == 'EP1':
    b_count_tp = 1 # End point run
else:
    b_count_tp = b_timepoints # number of timepoints, standard for 2h techdev is 25 atm

# Define variables based on inputs
b_exp_name = b_fcount + '_' + b_ifc
# csv_file = fcount + '_' + ifc + '_' + instrument_type +'_rawdata.csv' # csv file with raw data from Fluidigm RT PCR software
b_csv_file = b_in_csv

# Write a path to that excel sheet
b_layout_file = b_in_layout

b_out_folder = b_out_dir

if path.exists(b_layout_file) == False:
    raise Exception((b_layout_file + " doesn't exist"))

if path.exists(b_csv_file) == False:
    raise Exception((b_csv_file + " doesn't exist"))

# read in varctrl and ctrlvar ratio series
a_varctrl_series_ratio = pd.read_pickle(path.join(a_out_folder, a_exp_name+ '_' +a_instrument_type + '_varctrl_ratio_tpseries.pkl'))
a_ctrlvar_series_ratio = pd.read_pickle(path.join(a_out_folder, a_exp_name+ '_' +a_instrument_type + '_ctrlvar_ratio_tpseries.pkl'))

b_varctrl_series_ratio = pd.read_pickle(path.join(b_out_folder, b_exp_name+ '_' +b_instrument_type + '_varctrl_ratio_tpseries.pkl'))
b_ctrlvar_series_ratio = pd.read_pickle(path.join(b_out_folder, b_exp_name+ '_' +b_instrument_type + '_ctrlvar_ratio_tpseries.pkl'))

# concatenate dfs
varctrl_series_ratio = pd.concat([a_varctrl_series_ratio,b_varctrl_series_ratio])
ctrlvar_series_ratio = pd.concat([a_ctrlvar_series_ratio,b_ctrlvar_series_ratio])

varctrl_series_ratio.to_pickle(path.join(out_path,'varctrl_ratio_tpseries.pkl'))
ctrlvar_series_ratio.to_pickle(path.join(out_path,'ctrlvar_ratio_tpseries.pkl'))

# Constructs the control and variant guide pairs from assignment sheet tab
# tab "guide_map" which should be constructed manually
variant_gmap_a = pd.read_excel(path.join('',a_layout_file),sheet_name='guide_map')
variant_gmap_b = pd.read_excel(path.join('',b_layout_file),sheet_name='guide_map')
variant_gmap_a['group'] = '1'+variant_gmap_a['group']
variant_gmap_b['group'] = '2'+variant_gmap_b['group']

variant_gmap_ = pd.concat([variant_gmap_a,variant_gmap_b])
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
variant_snpmap = pd.read_excel(path.join('',a_layout_file),sheet_name='variant_SNPmap')
variant_gmap = pd.concat([variant_gmap,pd.DataFrame(columns=variant_snpmap.columns.values)])
variant_gmap['var_SNP'] = variant_gmap['variant_guide'].apply(lambda x: x.split('_vs_')[0])

for var in variant_snpmap:
    snps = variant_snpmap[var].dropna().values
    for snp in snps:
        variant_gmap.loc[variant_gmap.var_SNP.str.contains(snp),var] = 1
variant_gmap = variant_gmap.fillna(0)
snpmap_dict = variant_snpmap.to_dict('list')

# remove crRNAs
crRNA_remove = ['WGAN_SpikeRefE484K_0.1uM : wgan.E484Q_vs_Ref.1427','wgan.Ref_vs_P681R.2013 : GEN_P681H_vs_Spike_Ref-0.1uM']
variant_gmap = variant_gmap[~variant_gmap.index.isin(crRNA_remove)]
# change crRNA to work for both SNPs at same position
variant_gmap.loc['wgan.Ref_vs_P681R.2013 : wgan.P681R_vs_Ref.2013','var_SNP'] = 'wgan.P681R/H'

# add unique snp and variants columns
variants = ['Alpha','Beta','Gamma','Delta','Epsilon','Kappa','Eta','Lamda']
# variants = ['Alpha','Beta','Gamma','Delta','Epsilon']

variant_gmap['unique'] = variant_gmap[variants].sum(axis=1)
variant_gmap['variants'] = variant_gmap[variants].multiply(variants).values.tolist()
variant_gmap['variants'] = variant_gmap['variants'].apply(lambda y: [x for x in y if x])

# dictionary for seedstock names
seedstock_dict = {
'B.1.1.7' : 'Alpha',
'B.1.351' : 'Beta',
'P.1' : 'Gamma',
'B.1.617.2' : 'Delta',
'B.1.429' : 'Epsilon',
'B.1.427' : 'Epsilon',
'B.1.427/429' : 'Epsilon',
'Delta':'Delta'}

# plotting ground truth variant snp info
variant_gmap_ = variant_gmap[['Alpha', 'Beta', 'Gamma',
       'Delta', 'Epsilon', 'Kappa', 'Eta', 'Lamda', 'var_SNP']]
# variant_gmap_ = variant_gmap[['Alpha', 'Beta', 'Gamma',
#        'Delta', 'Epsilon', 'var_SNP']]
variant_gmap_.loc['sum',:] = variant_gmap_.sum()
variant_gmap_.iloc[-1,-1] = 'sum'
variant_gmap_.set_index('var_SNP',inplace=True)
variant_gmap_.sort_values(by='sum',axis=1,ascending=False,inplace=True)
variant_gmap_.drop(index='sum',inplace=True)
variant_gmap_.rename(columns={'B.1.427/429':'Epsilon'},inplace=True)
variant_gmap_ = variant_gmap_[~variant_gmap_.index.duplicated(keep='first')]

plt.figure(figsize=(10,3))
sns.heatmap(variant_gmap_.T,linewidths=0.05,linecolor='silver',cmap='Blues',cbar=False,square=True)
plt.title('Ground Truth SNPs for Variants');plt.xlabel('')
plt.savefig(path.join(out_path,'groundtruth_var_SNPab.png'),bbox_inches='tight',dpi=300)

#### VARIANT CALLING FUNCTIONS #####
def thresh_maxFC(seriesdf,thresh_dict,crsnp_dict):
    '''
    :seriesdf: FC time series df
    :thresh_dict: threshold dictionary with crRNA_pair as keys
    :crsnp_dict: SNP dictionary with crRNA_pair as keys
    returns :thresh_maxFC_df: with maxFC, var_SNP, threshold, threshFC_tp (min tp FC>=thresh)
    '''
    thresh_maxFC_df = seriesdf.reset_index().melt(id_vars='index',value_vars=seriesdf.columns.values)
    thresh_maxFC_df.rename(columns={'index':'crRNA_pair','variable':'sample_name','value':'FC_tpseries'},inplace=True)
    thresh_maxFC_df['maxFC'] = thresh_maxFC_df['FC_tpseries'].apply(lambda x: np.max(x))
    thresh_maxFC_df['var_SNP'] = thresh_maxFC_df['crRNA_pair'].map(crsnp_dict)
    thresh_maxFC_df['threshold'] = thresh_maxFC_df['crRNA_pair'].map(thresh_dict)
    thresh_maxFC_df['threshFC_idx'] = thresh_maxFC_df.apply(lambda x: np.append(np.argwhere(x.FC_tpseries>=x.threshold),timepoints+1),axis=1)
    thresh_maxFC_df['threshFC_mintp'] = thresh_maxFC_df['threshFC_idx'].apply(lambda x: np.min(x))
    return thresh_maxFC_df

def cxcomp(threshFC_df_,cxcomp_1a_,cxcomp_1b_):
    '''
    cross compare crRNA at the same position call the SNP that
    is has the highest FC value takes in :threshFC_df_: with
    crRNA SNP guide pairs :cxcomp_1a_: and :cxcomp_1b_: and
    returns :threshFC_df_: with modified SNP calls based on
    the cx comparison
    '''
    threshFC_df_cc = threshFC_df_.copy()
    cxcomp_1a_ct = threshFC_df_cc[(threshFC_df_cc.crRNA_pair==cxcomp_1a_)&(threshFC_df_cc.mintp_crRNA==\
                                              'var')].groupby('sample_name').count().reset_index()
    cxcomp_1b_ct = threshFC_df_cc[(threshFC_df_cc.crRNA_pair==cxcomp_1b_)&(threshFC_df_cc.mintp_crRNA==\
                                              'var')].groupby('sample_name').count().reset_index()
    cxcomp_1_ct = pd.concat([cxcomp_1a_ct,cxcomp_1b_ct]).groupby('sample_name').count()
    cxcomp_samples = cxcomp_1_ct[cxcomp_1_ct.crRNA_pair==2].index.values
    for i in cxcomp_samples:
        a = threshFC_df_cc.loc[(threshFC_df_cc.sample_name==i)&(threshFC_df_cc.crRNA_pair==cxcomp_1a_),'maxFC_varctrl'].values
        b = threshFC_df_cc.loc[(threshFC_df_cc.sample_name==i)&(threshFC_df_cc.crRNA_pair==cxcomp_1b_),'maxFC_varctrl'].values
        if a>b:
            threshFC_df_cc.loc[(threshFC_df_cc.sample_name==i)&(threshFC_df_cc.crRNA_pair==cxcomp_1b_),'mintp_crRNA'] = 'ctrl'
#             print(i,a,'>',b)
        elif a<b:
            threshFC_df_cc.loc[(threshFC_df_cc.sample_name==i)&(threshFC_df_cc.crRNA_pair==cxcomp_1a_),'mintp_crRNA'] = 'ctrl'
#             print(i,a,'<',b)
    return threshFC_df_cc

def seedstock_check(threshFC_sum):
    # Pull out seedstock samples
    threshFC_sum = threshFC_sum.reset_index()
    nsingle = threshFC_sum[(threshFC_sum.sample_name.str.contains('NSingle3'))]
    nsingle['variant'] = nsingle['sample_name'].apply(lambda x: x.split('_NSingle3_')[-1].split('_')[0])
    nsingle['variant'] = nsingle['variant'].map(seedstock_dict)
    nsingle['sample_dilution'] = nsingle['sample_name'].apply(lambda x: x.split('_')[-1])
    ref = threshFC_sum[(threshFC_sum.sample_name.str.contains('Nicole_reference'))]
    ref['variant'] = 'Ancestral'
    ref['sample_dilution'] = ref['sample_name'].apply(lambda x: x.split('_')[-1])
    delta = threshFC_sum[threshFC_sum.sample_name.str.contains('Delta')]
    delta['variant'] = 'Delta'
    delta['sample_dilution'] = delta['sample_name'].apply(lambda x: x.split('_')[-1].split('to')[0]+':'+x.split('_')[-1].split('to')[-1])
    delta.loc[delta.sample_name=='Delta','sample_dilution'] = '1'
    return pd.concat([ref,nsingle,delta])

def sen_spe(seedstock_ref_sum,variant_gmap_):
    # create varsnp_ref for sen and spe calc
    varsnp_ref = variant_gmap_.T
    # calc sen and spe
    senspe_df = pd.DataFrame(columns=['threshold','sen','spe'],index=varsnp_ref.columns)
    for i in varsnp_ref.columns:
        vardf = pd.DataFrame(varsnp_ref[i])
        vars_ = vardf[vardf[i]==1].index.values
        var_num = seedstock_ref_sum[(seedstock_ref_sum.variant.isin(vars_))].shape[0]
        senspe_df.loc[i,'threshold'] = variant_gmap.loc[variant_gmap.var_SNP==i,'FCthresh'].values[0]
        if (len(vars_)>0) & (var_num>0):
            tp = seedstock_ref_sum[(seedstock_ref_sum.variant.isin(vars_))&(seedstock_ref_sum[i]==1)].shape[0]
            fp = seedstock_ref_sum[~(seedstock_ref_sum.variant.isin(vars_))&(seedstock_ref_sum[i]==1)].shape[0]
            fn = seedstock_ref_sum[(seedstock_ref_sum.variant.isin(vars_))&(seedstock_ref_sum[i]!=1)].shape[0]
            tn = seedstock_ref_sum[~(seedstock_ref_sum.variant.isin(vars_))&(seedstock_ref_sum[i]!=1)].shape[0]
            sen = tp/(tp+fn)
            spe = tn/(tn+fp)
            senspe_df.loc[i,'sen'] = sen
            senspe_df.loc[i,'spe'] = spe
        else:
            senspe_df.loc[i,'sen'] = 'x'
            senspe_df.loc[i,'spe'] = 'x'
    return senspe_df

variants_ = ['Gamma', 'Beta', 'Alpha', 'Delta', \
             'Epsilon', 'Kappa', 'Eta','Lamda']
snps_ = variant_gmap_.index.values
unique_snp_map_ = variant_gmap[variant_gmap.unique==1].set_index('var_SNP')[['Gamma', 'Beta', 'Alpha', 'Delta', 'Epsilon', 'Kappa', 'Eta',
       'Lamda']]

def check_var(snp_pos_,var):
    # corr snp for var
    corr_vars_all = variant_gmap_.index[variant_gmap_[var].to_numpy().nonzero()[0]].tolist()
    corr_vars_uni = unique_snp_map_.index[unique_snp_map_[var].to_numpy().nonzero()[0]].tolist()
    corr_vars_shared = list(set(corr_vars_all).difference(corr_vars_uni))
    if var == 'Epsilon':
        corr_vars_all.remove('genetic.P26S')
        corr_vars_shared.remove('genetic.P26S')
    # sample calls
    uni_var_pos = list(set(corr_vars_uni) & set(snp_pos_))
    shared_var_pos = list(set(corr_vars_shared) & set(snp_pos_))
    off_var_pos = list(set(snp_pos_).difference(corr_vars_all))
    if len(uni_var_pos)==0:
        var_call_ = 0
    else:
        if len(off_var_pos)>2:
            var_call_ = 'Offtarget_'+str(len(off_var_pos))
        else:
            if (len(shared_var_pos)+len(uni_var_pos))>=(len(corr_vars_all)-2):
                var_call_ = len(uni_var_pos)
            else:
                var_call_ = 'Ontarget_'+str(len(uni_var_pos))+'_'+str(len(shared_var_pos))+'_'+str(len(corr_vars_all))
    return var_call_

# make dicts
cr_snp = dict(zip(variant_gmap.index,variant_gmap.var_SNP))
var_snpmap_total = variant_gmap_.sum().to_dict()

# import FC threshold dictionary
FCthresh = pd.read_csv(args.FCthresh_file)
thresh_dict = dict(zip(FCthresh.snp,FCthresh.thresh))

variant_gmap['FCthresh'] = variant_gmap['var_SNP'].map(thresh_dict)
thresholds = dict(zip(variant_gmap.index,variant_gmap.FCthresh))

# get varctrl and ctrlvar dfs for thresh_maxFC
varctrl_df = thresh_maxFC(varctrl_series_ratio,thresholds,cr_snp)
ctrlvar_df = thresh_maxFC(ctrlvar_series_ratio,thresholds,cr_snp)
# merge varctrl_df and ctrlvar_df
threshFC = varctrl_df.merge(ctrlvar_df,on=['crRNA_pair','sample_name'],suffixes=['_varctrl','_ctrlvar'])
# determine if varctrl or ctrlvar reached threshold first
threshFC['mintp_crRNA'] = 0
threshFC.loc[threshFC['threshFC_mintp_varctrl']<threshFC['threshFC_mintp_ctrlvar'],'mintp_crRNA'] = 'var'
threshFC.loc[threshFC['threshFC_mintp_varctrl']>threshFC['threshFC_mintp_ctrlvar'],'mintp_crRNA'] = 'ctrl'

# cross comparisons for same position SNPs
cxcomp_1a =  'GEN-Spike_Ref_vs_T478K-1 : GEN-T478K_vs_Spike_Ref-1'
cxcomp_1b =  'WGAN_SpikeRefE484K_0.1uM : WGAN_E484K_1uM'
cxcomp_2a =  'GEN_Spike_Ref_vs_K417N-K417T-5 : GEN_K417N_vs_K417T-Spike_Ref-1'
cxcomp_2b =  'GEN_Spike_Ref_vs_K417N-K417T-5 : GEN_K417T_vs_K417N-Spike_Ref-1'
cxcomp_3a =  'GEN_Spike_Ref_vs_N501T-N501Y-2.5 : GEN_N501T_vs_N501Y-Spike_Ref-2.5'
cxcomp_3b =  'GEN_Spike_Ref_vs_N501T-N501Y-2.5 : GEN_N501Y_vs_N501T-Spike_Ref-1'
cxcomp_4a =  'WGAN-Spike_Ref_vs_Q677H-1 : WGAN-Q677P_vs_Spike_Ref-1'
cxcomp_4b =  'WGAN-Spike_Ref_vs_Q677H-1 : GEN-Q677H_vs_Spike_Ref-1'

threshFC = cxcomp(threshFC,cxcomp_1a,cxcomp_1b)
threshFC = cxcomp(threshFC,cxcomp_2a,cxcomp_2b)
threshFC = cxcomp(threshFC,cxcomp_3a,cxcomp_3b)
threshFC = cxcomp(threshFC,cxcomp_4a,cxcomp_4b)

# adding snp columns and marking 0/1
var_snps = variant_gmap.var_SNP.unique()
threshFC[var_snps] = 0
# change value = 1 in variant column if var crRNA SNP max FC > ctrl crRNA SNP max FC
for snp in var_snps:
    threshFC.loc[(threshFC.var_SNP_varctrl==snp)&(threshFC.mintp_crRNA=='var'),snp] = 1

# remove sample_names that didn't have sample in them
others_to_remove=['Delta_1to1000','Delta_1to10000','Delta_1to100000']
threshFC = threshFC[threshFC.sample_name!='USA-RI-CDCBI-RIDOH_10437-2021']
threshFC = threshFC[~threshFC.sample_name.isin(others_to_remove)]
# summing with groupby sample_name
threshFC_sum = threshFC.groupby(['sample_name']).sum()[var_snps]

# Variant potential calling
variant_calling = pd.DataFrame(index=threshFC_sum.index.values,columns=variants_)
for sample in threshFC_sum.index:
    for var in variants_:
        sample_df = threshFC_sum[threshFC_sum.index==sample][snps_]
        sample_dfT = sample_df.T
        snp_pos_ = sample_dfT.index[sample_dfT.to_numpy().nonzero()[0]].tolist()
        var_call_ = check_var(snp_pos_,var)
        variant_calling.loc[sample,var] = var_call_

# Variant final calling
variant_calling['variant_call'] = 'VNI'
var_call_melt = variant_calling.reset_index().melt(id_vars='index').rename(columns={'index':'sample_name'})
for sample in variant_calling.index:
    sampledf = var_call_melt[var_call_melt.sample_name==sample]
    if sampledf[(sampledf.value==2)].shape[0]==1:
        variant_calling.loc[sample,'variant_call'] = sampledf[sampledf.value==2].variable.values[0]
    elif sampledf[(sampledf.value==1)].shape[0]==1:
        variant_calling.loc[sample,'variant_call'] = sampledf[sampledf.value==1].variable.values[0]
    elif sampledf[sampledf.value==0].shape[0]==len(variants_):
        variant_calling.loc[sample,'variant_call'] = 'VNI'
    elif sampledf[(sampledf.value==1)].shape[0]>1:
        vars_ = sampledf[(sampledf.value==1)].variable.values.tolist()
        vars_df = pd.DataFrame(index=vars_,columns=['shared','unique','offtarget'])
        sample_dfT = threshFC_sum[threshFC_sum.index==sample].T
        snp_pos_ = sample_dfT.index[sample_dfT.to_numpy().nonzero()[0]].tolist()
        for var_ in vars_df.index:
            unique_snp_map_sub = variant_gmap[variant_gmap.unique==1].set_index('var_SNP')[vars_]
            corr_vars_all = variant_gmap_.index[variant_gmap_[var_].to_numpy().nonzero()[0]].tolist()
            corr_vars_uni = unique_snp_map_sub.index[unique_snp_map_sub[var_].to_numpy().nonzero()[0]].tolist()
            corr_vars_shared = list(set(corr_vars_all).difference(corr_vars_uni))
            if var == 'Epsilon':
                corr_vars_all.remove('genetic.P26S')
                corr_vars_shared.remove('genetic.P26S')
            # sample calls
            uni_var_pos = list(set(corr_vars_uni) & set(snp_pos_))
            shared_var_pos = list(set(corr_vars_shared) & set(snp_pos_))
            off_var_pos = list(set(snp_pos_).difference(corr_vars_all))
            vars_df.loc[var_,'shared'] = len(uni_var_pos)
            vars_df.loc[var_,'unique'] = len(shared_var_pos)
            vars_df.loc[var_,'offtarget'] = len(off_var_pos)
        max_df = vars_df[vars_df['unique']==vars_df['unique'].max()]
        if max_df.shape[0]==1:
            variant_calling.loc[sample,'variant_call'] = max_df.index.values[0]
        else:
            vars_ = sampledf[(sampledf.value==1)].variable.values
            vars_str = ''
            for var__ in vars_:
                vars_str+=(var__+'_')
            variant_calling.loc[sample,'variant_call'] = vars_str
    else:
        variant_calling.loc[sample,'variant_call'] = 'VNI'

var_snp_call = variant_calling.merge(threshFC_sum[snps_],left_index=True,right_index=True)
seedstock_ref_sum = seedstock_check(threshFC_sum)
senspe_df = sen_spe(seedstock_ref_sum,variant_gmap_)

threshFC.to_csv(path.join(out_path,'threshFC.csv'))
threshFC.to_pickle(path.join(out_path,'threshFC.pkl'))
threshFC_sum.to_csv(path.join(out_path,'threshFC_sum.csv'))
threshFC_sum.to_pickle(path.join(out_path,'threshFC_sum.pkl'))
seedstock_ref_sum.to_csv(path.join(out_path,'seedstock_ref_sum.csv'))
seedstock_ref_sum.to_pickle(path.join(out_path,'seedstock_ref_sum.pkl'))
var_snp_call.to_csv(path.join(out_path,'variant_calls.csv'))
senspe_df.to_csv(path.join(out_path,'senspe_df.csv'))

seed_varcall = variant_calling[(variant_calling.index.str.contains('Delta'))\
                |(variant_calling.index.str.contains('Nicole'))|(variant_calling.index.str.contains('NSingle'))]
seed_df = seedstock_ref_sum.set_index('sample_name').merge(seed_varcall,left_index=True,right_index=True)
seed_df['sample_varcall'] = list(zip(seed_df.index,seed_df['variant_call']))
hm = seed_df.set_index('sample_varcall')[var_snps]
plt.figure(figsize=(0.5*hm.shape[1],0.5*hm.shape[0]))
sns.heatmap(hm,vmax=1,vmin=.99,cmap='Blues',square=True,xticklabels=True,yticklabels=True,\
                cbar=False,annot=True,linewidths=0.1,linecolor='silver')
plt.ylabel('')
plt.savefig(path.join(out_path,'seedstock_SNP_calls.pdf',bbox_inches='tight',dpi=300))
plt.savefig(path.join(out_path,'seedstock_SNP_calls.png',bbox_inches='tight',dpi=300))

# plot seedstock sen and spe
senspe_df_melt = senspe_df.reset_index().melt(id_vars='var_SNP')
senspe_df_melt = senspe_df_melt[(senspe_df_melt.variable!='threshold')&(senspe_df_melt.value!='x')]
plt.figure()
sns.barplot(x="var_SNP", y="value", hue="variable", data=senspe_df_melt)
plt.xticks(rotation=90)
plt.ylabel('Sensitivity / Specificity')
plt.legend(bbox_to_anchor=(1,1));plt.xlabel('')
plt.savefig(path.join(out_path,'seedstock_senspe.pdf',bbox_inches='tight',dpi=300))
plt.savefig(path.join(out_path,'seedstock_senspe.png',bbox_inches='tight',dpi=300))
