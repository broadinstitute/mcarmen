import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

def PlotHeatmap(median_frame, hit_df, tp, assaylist, samplelist, prefix):    
    frame = median_frame[samplelist].reindex(assaylist)
    fig, axes = plt.subplots(1, 1, figsize=(len(frame.columns.values)*0.5, len(frame.index.values)*0.5))
        
    if hit_df is not None:
        # Change positive, negative and invalid in hit dataframe to "+", "" and "!"
        hit_df = hit_df.replace({'positive': "+", 'negative':'', 'invalid':'!'})
        ax = sns.heatmap(frame, annot = hit_df, cmap='Reds', fmt = "s", square = True, annot_kws={"size": 12, "color": 'black'}) # cbar_kws={'pad':0.002} 
    else:
        ax = sns.heatmap(frame, cmap='Reds', fmt = "s", square = True)
    
    plt.title('Median values and hits for {}'.format(tp), size = 14)
    
    plt.xlabel('Samples', size = 12)
    plt.ylabel('Assays', size = 12)
    
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom + 0.5, top - 0.5)
    ax.tick_params(axis="y", labelsize=12)
    ax.tick_params(axis="x", labelsize=12)
    plt.yticks(rotation=0) 
    
    h_lines = np.arange(3, len(assaylist), 3)
    v_lines = np.arange(3, len(samplelist), 3)
    axes.hlines(h_lines, colors = 'silver', alpha=0.9, linewidths = 0.35, *axes.get_xlim())
    axes.vlines(v_lines, colors = 'silver', alpha=0.9, linewidths = 0.35, *axes.get_ylim())

    plt.savefig("{}_heatmap_{}.png".format(prefix,tp), format = 'png', bbox_inches='tight', dpi=300)
    plt.close()         
