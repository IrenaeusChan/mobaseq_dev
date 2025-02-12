import traceback, os, sys
import mobaseq.utils.logger as log
import mobaseq.process.tools as tools
import mobaseq.qc.summarize as summarize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing as mp
import regex as re
from adjustText import adjust_text
from math import sqrt
from itertools import cycle

# Function to extract tissue and number
def sort_key(key):
    # Define tissue sorting order
    tissue_order = {"Liver": 1, "Liv": 1, "Lung": 2, "Lug" : 2, "Brain": 3, "Brn" : 3, "Blood": 4, "WB": 5, "BM": 6, "preTran": 7}
    tissue = re.match(r"(\w+)", key)(1)  # Extract tissue type
    number = int(re.match(r"(\d+)", key)(1))  # Extract number
    if tissue != None or number != None:
        tissue_rank = tissue_order.get(tissue, 99)  # Assign tissue rank, default to 99 if not found
        return (tissue_rank, number)  # Sort by tissue, then by number
    return (99, 99999)  # Default sorting for unexpected values

def mapped_reads(mapped_percentages_csv, out_dir, debug):
    log.logit(f"Plotting % mapped reads", color = "green")
    df = pd.read_csv(mapped_percentages_csv)
    # Plotting
    plt.figure(figsize=(12, 10))
    bar_width = 0.35
    mapped_bar = plt.bar(df["Sample"], df["Mapped %"], label="Mapped %", color="blue", alpha=0.7, width=bar_width)
    # Use the bottom parameter to stack the bars, in this case, the Unmapped % bars will be stacked on top of the Mapped % bars
    unmapped_bar = plt.bar(df["Sample"], df["Unmapped %"], label="Unmapped %", color="red", alpha=0.7, bottom=df["Mapped %"], width=bar_width)
    # Outline bars with mapped percentage under 30% in red
    for bar in mapped_bar:
        if bar.get_height() < 30:
            bar.set_edgecolor('red')
            bar.set_linewidth(2)
    # Add percentage numbers on top of the bars
    for bar in mapped_bar:
        height = bar.get_height()
        # The text() function takes the x and y coordinates of the text, and the text itself
        plt.text(bar.get_x() + bar.get_width() / 2.0, height, f'{height:.1f}%', ha='center', va='bottom')
    for bar1, bar2 in zip(mapped_bar, unmapped_bar):
        height = bar2.get_height()
        plt.text(bar2.get_x() + bar2.get_width() / 2.0, (height+bar1.get_height()+ 0.1), f'{height:.1f}%', ha='center', va='bottom')
    plt.xlabel("Sample")
    plt.ylabel("Percentage")
    plt.title("Mapped % and Unmapped % per Sample")
    plt.xticks(rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{out_dir}/mapped_percentages.png")
    plt.show()

def spike_ins(spike_count, slope, rsq, plot_color, sample_name, out_dir, redraw = False):
    log.logit(f"Plotting spike-in counts vs. expected cell numbers", color = "green")
    fig, ax = plt.subplots(figsize=(12, 10))
    # Scatter plot
    ax.scatter(spike_count['count'], spike_count['expected_cellnum'], color=plot_color)
    # Regression line
    x = spike_count['count']
    y_pred = slope * x
    ax.plot(x, y_pred, color=plot_color)
    
    # Add labels with z-scores
    texts = []
    for _, row in spike_count.iterrows():
        texts.append(ax.text(row['count'], row['expected_cellnum'], 
                           f"{row['name']} (z={row['z_scores']:.2f})"))
    
    # Adjust text positions to avoid overlap
    adjust_text(texts)
    
    # Add regression equation and R² text
    eq_text = f"y = {slope:.5f}x  Rsq = {rsq:.3f}"
    ax.text(min(x), max(spike_count['expected_cellnum']), eq_text,
            horizontalalignment='left', verticalalignment='top')
    
    # Set title and labels
    if redraw:
        ax.set_title(f'Spike-in Counts vs. Expected Cell Numbers - {sample_name} - After Dropping Highest Z-Score')
    else:
        ax.set_title(f'Spike-in Counts vs. Expected Cell Numbers - {sample_name}')
    ax.set_xlabel('Spike-in Counts')
    ax.set_ylabel('Expected Cell Numbers')
        
    # Show plot
    plt.grid(True)
    plt.tight_layout()
    if redraw:
        plt.savefig(f"{out_dir}/{sample_name}_R-SquaredPlot_Redraw.png")
    else:
        plt.savefig(f"{out_dir}/{sample_name}_R-SquaredPlot.png")
    plt.show()
    return fig

def rsq_per_sample(Rsq, out_dir):
    log.logit(f"Plotting R² values per sample", color = "green")
    # Plotting
    plt.figure(figsize=(12, 10))
    plt.barh(Rsq["Sample_ID"], Rsq["Spike_Rsq"], color="blue", alpha=0.7, height=0.9)
    plt.axvline(x=0.9, color='red', linestyle='--', linewidth=1)
    
    plt.xlim(0.6, 1.0)
    plt.ylabel("Sample ID")
    plt.xlabel("R²")
    plt.title("Spike-In R² Values per Sample")
    #plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/R-squared_per_sample.png")
    plt.show()

def reading_depth_per_sample(df, out_dir):
    log.logit(f"Plotting reading depth per sample", color = "green")
    # Plotting
    plt.figure(figsize=(12, 10))
    plt.barh(df["Sample_ID"], df["reading_depth"], color="blue", alpha=0.7, height=0.9)
    plt.axvline(x=0.003, color='red', linestyle='--', linewidth=1)

    plt.xlim(0, 0.1)
    plt.ylabel("Sample ID")
    plt.xlabel("Reading Depths")
    plt.title("Reading Depth per Sample")
    #plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/reading_depth_per_sample.png")
    plt.show()

def mean_reading_depth_distribution(aggr_df, out_dir):
    log.logit(f"Plotting mean reading depth per sample", color = "green")
    sns.set_style("whitegrid")
    # Plotting
    plt.figure(figsize=(12, 10))
    max_reading_depth = round(aggr_df['mean_ReadingDepth'].max())
    sns.kdeplot(data=aggr_df, 
                x="mean_ReadingDepth",
                fill=True, 
                color="blue",
                alpha=0.5,
                clip = (0, max_reading_depth),
                bw_adjust=0.3)
    plt.axvline(x=2, color='red', linestyle='--', linewidth=1)
    plt.text(3, plt.ylim()[1]*0.95, 'x=2', color='red')
    
    plt.xlabel('Reading Depth')
    plt.ylabel('Density')
    plt.title('Density Plot of Reading Depth per Sample')
    plt.tight_layout()
    plt.savefig(f"{out_dir}/mean_reading_depth_density.png")
    plt.show()

#ggplot(agg_data, aes(x = mean_ReadingDepth)) +
#  geom_density(fill = "blue", alpha = 0.5) +
#  geom_vline(xintercept = 2, linetype = "dashed", color = "red", size = 1) +
#  annotate("text", x = 3, y = 0.087, label = "x=2", color = "red") +
#  labs(title = "Density Plot of Reading Depth per Sample", x = "Reading Depth", y = "Density")

def mean_rsq_distribution(aggr_df, out_dir):
    log.logit(f"Plotting R² of cell count distribution", color = "green")
    # Plotting
    plt.figure(figsize=(12, 10))
    plt.hist(aggr_df['mean_Rsq'], 
        bins=np.arange(0.1, 1.1, 0.01),  # binwidth = 0.01
        color='red',
        edgecolor='black',
        alpha=0.5)
    lowest_Rsq = round(aggr_df[aggr_df['mean_Rsq'] > 0]['mean_Rsq'].min(), 2)
    plt.xlim(lowest_Rsq, 1)
    
    plt.xlabel("R²")
    plt.ylabel("Cell Count")
    plt.title("Histogram of R² of Cell Count Distribution Calculated from Spike In Regression")
    #plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/rsq_of_cell_count_distribution.png")
    plt.show()

def check_all_sgIDs(df1, df2, df3, df4):
    # Check if all sgIDs are the same
    if not df1['sgID'].equals(df2['sgID']) or not df2['sgID'].equals(df3['sgID']) or not df3['sgID'].equals(df4['sgID']):
        log.logit("The sgIDs in the Dataframes are not equal")
        raise ValueError("Dataframes are not equal")
    return True

def per_sgid(total_reads_df, unique_barcodes_df, relative_reads_df, relative_barcodes_df, out_dir, no_dummy, debug):
    log.logit(f"Plotting the relative distribution of reads and barcodes per sgID by sample.")
    check_all_sgIDs(total_reads_df, unique_barcodes_df, relative_reads_df, relative_barcodes_df)
    dynamic_threshold = 1/(2*sqrt(len(total_reads_df)))
    all_sgIDs = total_reads_df['sgID'].unique()
    all_samples = total_reads_df.columns[1:].values
    if debug: 
        log.logit(f"Dynamic threshold: {dynamic_threshold}")
        log.logit(f"all_sgIDs: {all_sgIDs}")
        log.logit(f"all_samples: {all_samples}")
    
    # Make dataframes long format and preserve column names
    dfs = {
        'total_reads': total_reads_df,
        'barcodes': unique_barcodes_df,
        'rel_reads': relative_reads_df,
        'rel_barcodes': relative_barcodes_df
    }
    
    for key in dfs:
        # Set sgID as index before transpose
        dfs[key] = dfs[key].set_index('sgID').T

    # Plotting all figures side by side
    fig, axs = plt.subplots(2, 2, figsize=(24, 16))  # 2 rows, 2 columns

    # Initialize bottom position for stacking
    bottoms = [None] * 4

     # Plot stacked bars
    for idx, sgID in enumerate(all_sgIDs):
        if bottoms[0] is None:
            # First bars
            for i, (key, ax) in enumerate(zip(dfs.keys(), axs.flat)):
                bottoms[i] = dfs[key][sgID]
                ax.bar(all_samples, dfs[key][sgID], label=sgID)
        else:
            # Stack subsequent bars
            for i, (key, ax) in enumerate(zip(dfs.keys(), axs.flat)):
                ax.bar(all_samples, dfs[key][sgID], bottom=bottoms[i], label=sgID)
                bottoms[i] += dfs[key][sgID]
        # Add labels for significant values
        for i, (key, ax) in enumerate(zip(dfs.keys(), axs.flat)):
            if ax == axs[1, 0] or ax == axs[1, 1]:
                for j, value in enumerate(dfs[key][sgID]):
                    if value > dynamic_threshold:
                        ax.text(j, bottoms[i][j] - value / 2, sgID, ha='center', va='center', fontsize=8, color='black')

    # Customize plot for total reads per sgID
    axs[0, 0].set_xlabel('Samples')
    axs[0, 0].set_ylabel('Total Reads')
    axs[0, 0].set_title('Total Reads per sgID by Sample')
    axs[0, 0].tick_params(axis='x', rotation=75)
    # Customize plot for unique barcodes per sgID
    axs[0, 1].set_xlabel('Samples')
    axs[0, 1].set_ylabel('Counts')
    axs[0, 1].set_title('Unique Barcodes per sgID by Sample')
    axs[0, 1].tick_params(axis='x', rotation=75)
    # Customize plot for relative reads per sgID
    axs[1, 0].set_xlabel('Samples')
    axs[1, 0].set_ylabel('Relative Counts')
    axs[1, 0].set_title('Relative Reads per sgID by Sample')
    axs[1, 0].tick_params(axis='x', rotation=75)
    # Customize plot for relative barcodes per sgID
    axs[1, 1].set_xlabel('Samples')
    axs[1, 1].set_ylabel('Relative Counts')
    axs[1, 1].set_title('Relative Barcodes per sgID by Sample')
    axs[1, 1].tick_params(axis='x', rotation=75)
    
    # handles, labels = axs[0, 0].get_legend_handles_labels()
    # num_cols = min(18, len(labels))
    # fig.legend(handles, labels, loc='upper center', 
    #           bbox_to_anchor=(0.5, 1.), ncol=num_cols,
    #           title='sgIDs', fontsize=10, title_fontsize=13)

    # Adjust layout and save figure
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to make space for the legend
    suffix = '_noDummy' if no_dummy else ''
    save_name = os.path.join(out_dir, f'sgID_distribution_with_relative_reads_and_barcodes_filtered{suffix}.pdf')
    plt.savefig(save_name, format='pdf', bbox_inches='tight')
    plt.show()