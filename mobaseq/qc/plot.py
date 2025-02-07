import traceback
import mobaseq.utils.logger as log
import mobaseq.process.tools as tools
import mobaseq.qc.summarize as summarize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
import regex as re

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
