# Usage: python plot_controlled_meta_all_cells.py <amaranth_meta_dir> <aletsch_dir> <amaranth_meta_beaver_dir> <aletsch_beaver_dir>
import matplotlib.pyplot as plt
import re
import os
import glob
from sys import argv
import numpy as np

def natural_sort_key(s):
    """Convert a string into a list of mixed strings and integers for natural sorting."""
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

def read_roc_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    correct_values = []
    precision_values = []
    for line in lines:
        if line.startswith("#"):
            continue
        correct = int(line.split("correct = ")[1].split()[0])
        precision = float(line.split("precision = ")[1].split()[0])
        correct_values.append(correct)
        precision_values.append(precision)
    return correct_values, precision_values

# axis controlled: 0 -> recall, 1 -> precision
def controlled_1_cell(axis, methods, simple=False):
    # get lowest value
    lowest_precision = min(max(x[axis]) for x in methods)
    print('lowest_precision', lowest_precision)
    print('methods',    methods)
    # get controlled dependent value
    controlled_dependent_values = []
    for method in methods:
        k,v = -1,-1
        for row in zip(*method):
            if row[axis] >= lowest_precision:
                v = row[1-axis]
                k = row[axis]
                print('k', k, 'v', v)
            elif row[axis] == lowest_precision:
                controlled_dependent_values.append(row[1-axis])
                break
            else:
                assert k != -1 and v != -1
                assert k >= lowest_precision
                if not simple:
                    vv = ((lowest_precision - row[axis]) * row[1-axis] + (k - lowest_precision) * v) / (k - row[axis])
                    controlled_dependent_values.append(vv)
                else:
                    controlled_dependent_values.append((v + row[1-axis]) / 2)
                break
        else:
            print(method, f"no hitting of lowest {lowest_precision}")
            # If we went through all rows without breaking, use the last values
            # This happens when all values are >= lowest_precision
            assert k != -1 and v != -1
            controlled_dependent_values.append(v)
    print('controlled_dependent_values', controlled_dependent_values)
    assert len(controlled_dependent_values) == len(methods)
    return controlled_dependent_values

def get_method_data(method_dir, pattern):
    """Get data for a specific method by parsing .roc files."""
    matching_files = glob.glob(os.path.join(method_dir, pattern))

    # Filter out files containing '.nr.' in the basename
    filtered_files = []
    for f in matching_files:
        basename = os.path.basename(f)
        if '.nr.' not in basename and 'meta' not in basename:
            filtered_files.append(f)

    matching_files = sorted(filtered_files, key=natural_sort_key)

    if not matching_files:
        print(f"Warning: No files found in directory {method_dir} with pattern {pattern}")
        return []

    print(f"\nDirectory: {method_dir}")
    print(f"Found {len(matching_files)} files")
    print("First 10 files:")
    for i, f in enumerate(matching_files):
        print(f"  {i+1}. {os.path.basename(f)}")

    data = []
    for filepath in matching_files:
        correct_values, precision_values = read_roc_file(filepath)
        # Return as (correct_values, precision_values) tuple for each file
        data.append((correct_values, precision_values))

    return data

def plot_scatter(axis, method_configs, method_display_names, colors, simple=False):
    plt.figure(figsize=(5, 5))

    methods = list(method_configs.keys())

    # Collect all data per method
    all_method_data = {}
    for method, config in method_configs.items():
        data = get_method_data(config['dir'], config['pattern'])
        all_method_data[method] = data

    # Verify all methods have same number of files
    file_counts = [len(data) for data in all_method_data.values()]
    if len(set(file_counts)) > 1:
        print(f"Error: Methods have different number of files: {dict(zip(methods, file_counts))}")
        for method in methods:
            print(f"{method}: {file_counts[methods.index(method)]} files")
        raise AssertionError("All methods must have the same number of input files")

    num_cells = min(file_counts) if file_counts else 0

    all_controlled_dependent_values = []
    for i in range(num_cells):
        # For each cell, collect (correct_values, precision_values) data from all methods
        cell_data = []
        for method in methods:
            correct_values, precision_values = all_method_data[method][i]
            # Store as (correct_values, precision_values) to match controlled_1_cell expected format
            cell_data.append((correct_values, precision_values))

        controlled_dependent_values = controlled_1_cell(axis, cell_data, simple)
        all_controlled_dependent_values.append(controlled_dependent_values)

    # Sort by first method (amaranth_meta)
    all_controlled_dependent_values.sort(key=lambda x: x[0])

    # Plot each method
    for idx, method in enumerate(methods):
        y = [x[idx] for x in all_controlled_dependent_values]
        x_vals = list(range(len(y)))
        plt.scatter(x_vals, y, s=10, zorder=2, c=[colors[method]], label=method_display_names[method], alpha=0.6)
        # Draw horizontal line for the average, formatted to 2 decimals
        avg = np.mean(y)
        avg_fmt = f"{avg:.2f}"
        print("avg_fmt", avg_fmt, "method", method)

    plt.xlabel(f'cells (n={num_cells})', fontsize=14)
    plt.ylabel('adj precision (%)' if axis == 0 else 'adj # matching transcripts', fontsize=14)
    plt.xticks(fontsize=12, rotation=45, ha='right')
    plt.yticks(fontsize=12, rotation=90, va='center')
    # plt.grid(True, alpha=0.3)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    # Save main figure without legend
    plt.savefig(f'controlled_all_cells_{axis}_1.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'controlled_all_cells_{axis}_1.pdf', dpi=300, bbox_inches='tight')

    # Save legend separately
    handles, labels = plt.gca().get_legend_handles_labels()
    fig_legend = plt.figure(figsize=(8, 0.5))
    legend = fig_legend.legend(handles, labels, loc='center', ncol=4, fontsize=12, frameon=True, fancybox=True, markerscale=2)
    fig_legend.savefig(f'controlled_all_cells_{axis}_1_legend.png', dpi=300, bbox_inches='tight')
    fig_legend.savefig(f'controlled_all_cells_{axis}_1_legend.pdf', dpi=300, bbox_inches='tight')
    plt.close(fig_legend)
    plt.clf()

if __name__ == "__main__":
    # Get directories from command line
    if len(argv) != 4:
        print("Usage: python plot_controlled_meta_all_cells.py <amaranth_meta_dir> <aletsch_dir> <psiclass_dir>")
        exit(1)

    amaranth_meta_dir = argv[1]
    aletsch_dir = argv[2]
    psiclass_dir = argv[3]

    # Define file patterns and directories for the methods
    method_configs = {
        'amaranth_meta': {
            'dir': amaranth_meta_dir,
            'pattern': '*.roc'
        },
        'aletsch': {
            'dir': aletsch_dir,
            'pattern': '*.roc'
        },
        'psiclass': {
            'dir': psiclass_dir,
            'pattern': '*.roc'
        }
    }

    # Method display names
    method_display_names = {
        'amaranth_meta': 'Amaranth-meta',
        'aletsch': 'Aletsch',
        'psiclass': 'PsiClass'
    }

    # Colors for each method
    colors = {
        'amaranth_meta': '#FF9500',
        'aletsch': '#52b8d6',
        'psiclass': "#b788e4"
    }

    plot_scatter(0, method_configs, method_display_names, colors)
