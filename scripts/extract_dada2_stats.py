import os
import glob
import pandas as pd
import argparse
from tabulate import tabulate

# ANSI color codes
RED = "\033[91m"
RESET = "\033[0m"

# Argument parser
parser = argparse.ArgumentParser(description="Parse dada2_stats.tsv files and compute statistics.")
parser.add_argument("--working-directory", "-w", required=True, help="Base directory to search for dada2_stats.tsv files.")
args = parser.parse_args()

# Construct file pattern
base_dir = os.path.expanduser(args.working_directory)
pattern = os.path.join(base_dir, "**", "dada2_stats.tsv")
files = glob.glob(pattern, recursive=True)

# Columns to analyze
metrics = [
    'percentage of input passed filter',
    'percentage of input merged',
    'percentage of input non-chimeric'
]

# Prepare data for tabulate
headers = ["Sample_Path"] + [f"{m} (Min%)\t{m} (Max%)\t{m} (Mean%)" for m in metrics]
table_data = []

# Process each file
for file in files:
    try:
        df = pd.read_csv(file, sep="\t", comment="#")
        result_line = [os.path.dirname(file)]

        for col in metrics:
            if col in df.columns:
                min_val = df[col].min()
                max_val = df[col].max()
                mean_val = df[col].mean()

                # Highlight mean in red if < 70
                mean_str = f"{mean_val:.2f}"
                if mean_val < 70:
                    mean_str = f"{RED}{mean_str}{RESET}"

                result_line += [f"{min_val:.2f}", f"{max_val:.2f}", mean_str]
            else:
                result_line += ["MISSING", "MISSING", "MISSING"]
        table_data.append(result_line)
    except Exception as e:
        table_data.append([os.path.dirname(file), f"ERROR: {e}"])

# Print the table using tabulate
print(tabulate(table_data, headers=headers, tablefmt="grid"))


