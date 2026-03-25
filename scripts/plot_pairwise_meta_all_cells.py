"""
Pairwise controlled precision comparison for meta-assembly methods.

Accepts pairs of <dir> <method> arguments. For each pair of methods, generates a
scatter plot of adjusted precision values (controlled at the lower sensitivity)
across all cells. Also generates an unsorted version with sensitivity as x-axis.

Usage:
  python plot_pairwise_meta_all_cells.py <dir1> <method1> <dir2> <method2> ... \\
      [--output-prefix pairwise]

Known method names (for automatic file pattern matching):
  amaranth_meta, aletsch, psiclass
"""
import argparse
import matplotlib.pyplot as plt
import re
import os
import glob
import sys
import numpy as np


def natural_sort_key(s):
    """Convert a string into a list of mixed strings and integers for natural sorting."""
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]


def read_roc_file(filename):
    """Read ROC-style file with correct and precision values."""
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


def controlled_pairwise_precision(method1_data, method2_data, simple=False):
    """
    Given two methods' (correct, precision) data, control for sensitivity (correct values).
    Returns adjusted precision values for both methods at the lower sensitivity.
    """
    correct1, precision1 = method1_data
    correct2, precision2 = method2_data

    max_correct1 = max(correct1)
    max_correct2 = max(correct2)
    target_correct = min(max_correct1, max_correct2)

    adjusted_precision1 = None
    for i in range(len(correct1)):
        if correct1[i] >= target_correct:
            adjusted_precision1 = precision1[i]
        elif correct1[i] < target_correct:
            if i == 0:
                adjusted_precision1 = precision1[i]
            else:
                prev_correct = correct1[i-1]
                prev_precision = precision1[i-1]
                curr_correct = correct1[i]
                curr_precision = precision1[i]
                if not simple:
                    adjusted_precision1 = prev_precision + (target_correct - prev_correct) * \
                                         (curr_precision - prev_precision) / (curr_correct - prev_correct)
                else:
                    adjusted_precision1 = (prev_precision + curr_precision) / 2
            break

    if adjusted_precision1 is None:
        print("no hitting prec1")
        adjusted_precision1 = precision1[-1]

    adjusted_precision2 = None
    for i in range(len(correct2)):
        if correct2[i] >= target_correct:
            adjusted_precision2 = precision2[i]
        elif correct2[i] < target_correct:
            if i == 0:
                adjusted_precision2 = precision2[i]
            else:
                prev_correct = correct2[i-1]
                prev_precision = precision2[i-1]
                curr_correct = correct2[i]
                curr_precision = precision2[i]
                if not simple:
                    adjusted_precision2 = prev_precision + (target_correct - prev_correct) * \
                                         (curr_precision - prev_precision) / (curr_correct - prev_correct)
                else:
                    adjusted_precision2 = (prev_precision + curr_precision) / 2
            break

    if adjusted_precision2 is None:
        print("no hitting prec2")
        adjusted_precision2 = precision2[-1]

    return adjusted_precision1, adjusted_precision2, target_correct


def get_method_data(method_dir, pattern):
    """Get data for a specific method by parsing ROC files."""
    matching_files = glob.glob(os.path.join(method_dir, pattern))

    filtered_files = [f for f in matching_files
                      if '.nr.' not in os.path.basename(f) and 'meta' not in os.path.basename(f)]

    matching_files = sorted(filtered_files, key=natural_sort_key)

    if not matching_files:
        print(f"Warning: No files found in directory {method_dir} with pattern {pattern}")
        return []

    print(f"\nDirectory: {method_dir}")
    print(f"Found {len(matching_files)} files")

    data = []
    for filepath in matching_files:
        correct_values, precision_values = read_roc_file(filepath)
        data.append((correct_values, precision_values))

    return data


def plot_pairwise_comparison(method1_name, method2_name, method1_data, method2_data,
                             method_display_names, colors, output_prefix):
    """Create pairwise scatter plot comparing two methods."""
    fig, ax = plt.subplots(figsize=(5, 5))

    if len(method1_data) != len(method2_data):
        print(f"Error: {method1_name} has {len(method1_data)} files, {method2_name} has {len(method2_data)} files")
        return

    num_cells = len(method1_data)

    adjusted_precision1_list = []
    adjusted_precision2_list = []
    controlled_sensitivity_list = []

    for i in range(num_cells):
        correct1, precision1 = method1_data[i]
        correct2, precision2 = method2_data[i]
        adj_precision1, adj_precision2, target_sensitivity = controlled_pairwise_precision(
            (correct1, precision1),
            (correct2, precision2)
        )
        adjusted_precision1_list.append(adj_precision1)
        adjusted_precision2_list.append(adj_precision2)
        controlled_sensitivity_list.append(target_sensitivity)

    sorted_indices = np.argsort(adjusted_precision1_list)
    adjusted_precision1_sorted = [adjusted_precision1_list[i] for i in sorted_indices]
    adjusted_precision2_sorted = [adjusted_precision2_list[i] for i in sorted_indices]
    controlled_sensitivity_sorted = [controlled_sensitivity_list[i] for i in sorted_indices]

    x_vals = list(range(num_cells))

    ax.scatter(x_vals, adjusted_precision1_sorted,
               color=colors[method1_name], s=10, zorder=2,
               label=method_display_names[method1_name], alpha=0.6)
    ax.scatter(x_vals, adjusted_precision2_sorted,
               color=colors[method2_name], s=10, zorder=2,
               label=method_display_names[method2_name], alpha=0.6)

    avg1 = np.mean(adjusted_precision1_sorted)
    avg2 = np.mean(adjusted_precision2_sorted)
    avg_sensitivity = np.mean(controlled_sensitivity_sorted)

    method1_higher = sum(1 for i in range(num_cells) if adjusted_precision1_list[i] > adjusted_precision2_list[i])
    method2_higher = sum(1 for i in range(num_cells) if adjusted_precision2_list[i] > adjusted_precision1_list[i])
    equal = sum(1 for i in range(num_cells) if adjusted_precision1_list[i] == adjusted_precision2_list[i])

    print(f"\n{method1_name} avg adjusted precision: {avg1:.2f}%")
    print(f"{method2_name} avg adjusted precision: {avg2:.2f}%")
    print(f"Controlled sensitivity: {avg_sensitivity:.2f}")
    print(f"\nSample counts (n={num_cells}):")
    print(f"  {method1_name} > {method2_name}: {method1_higher} samples ({method1_higher/num_cells*100:.1f}%)")
    print(f"  {method2_name} > {method1_name}: {method2_higher} samples ({method2_higher/num_cells*100:.1f}%)")
    if equal > 0:
        print(f"  Equal: {equal} samples ({equal/num_cells*100:.1f}%)")

    ax.set_xlabel(f'cells (n={num_cells})', fontsize=14)
    ax.set_ylabel('adj precision (%)', fontsize=14)
    ax.set_ylim(0, 100)
    plt.xticks(fontsize=12, rotation=45, ha='right')
    plt.yticks(fontsize=12, rotation=90, va='center')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    plt.savefig(f"{output_prefix}.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches='tight')
    plt.close(fig)

    # Unsorted version with sensitivity as x-axis
    fig2, ax2 = plt.subplots(figsize=(5, 5))

    ax2.scatter(controlled_sensitivity_list, adjusted_precision1_list,
                color=colors[method1_name], s=10, zorder=2,
                label=method_display_names[method1_name], alpha=0.6)
    ax2.scatter(controlled_sensitivity_list, adjusted_precision2_list,
                color=colors[method2_name], s=10, zorder=2,
                label=method_display_names[method2_name], alpha=0.6)

    ax2.set_xlabel('adj # matching transcripts', fontsize=14)
    ax2.set_ylabel('adj precision (%)', fontsize=14)
    ax2.set_ylim(0, 100)
    plt.xticks(fontsize=12, rotation=45, ha='right')
    plt.yticks(fontsize=12, rotation=90, va='center')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    plt.tight_layout()

    plt.savefig(f"{output_prefix}_unsorted.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}_unsorted.png", dpi=300, bbox_inches='tight')
    plt.close(fig2)


# Predefined file patterns for known methods
PATTERNS = {
    'amaranth_meta': '*.roc',
    'aletsch': '*.roc',
    'psiclass': '*.roc',
}

METHOD_DISPLAY_NAMES = {
    'amaranth_meta': 'Amaranth-meta',
    'aletsch': 'Aletsch',
    'psiclass': 'PsiClass',
}

COLORS = {
    'amaranth_meta': "#FF9500",
    'aletsch': '#52b8d6',
    'psiclass': "#b788e4",
}


if __name__ == "__main__":
    args_list = sys.argv[1:]
    output_prefix = 'pairwise'
    if '--output-prefix' in args_list:
        idx = args_list.index('--output-prefix')
        output_prefix = args_list[idx + 1]
        args_list = args_list[:idx] + args_list[idx+2:]

    if len(args_list) < 4 or len(args_list) % 2 != 0:
        print("Usage: python plot_pairwise_meta_all_cells.py <dir1> <method1> <dir2> <method2> ... [--output-prefix <prefix>]")
        sys.exit(1)

    num_methods = len(args_list) // 2
    methods = []
    directories = []
    for i in range(num_methods):
        directories.append(args_list[2*i])
        methods.append(args_list[2*i + 1])

    print("Loading data for all methods...")
    all_method_data = {}
    for i, method in enumerate(methods):
        pattern = PATTERNS.get(method, '*.roc')
        data = get_method_data(directories[i], pattern)
        all_method_data[method] = data

    # Combined legend for all methods
    fig_combined_legend = plt.figure(figsize=(12, 0.5))
    ax_legend = fig_combined_legend.add_subplot(111)
    for spine in ax_legend.spines.values():
        spine.set_visible(False)
    ax_legend.set_xticks([])
    ax_legend.set_yticks([])

    seen = {}
    dummy_handles = []
    dummy_labels = []
    for method in methods:
        display_name = METHOD_DISPLAY_NAMES.get(method, method)
        if display_name not in seen:
            handle = plt.scatter([], [], s=10, color=COLORS.get(method, '#000000'), alpha=0.6)
            dummy_handles.append(handle)
            dummy_labels.append(display_name)
            seen[display_name] = method

    ax_legend.legend(dummy_handles, dummy_labels, loc='center',
                     fontsize=12, ncol=len(dummy_labels), frameon=True, markerscale=2,
                     edgecolor='black', fancybox=True)
    fig_combined_legend.tight_layout()
    fig_combined_legend.savefig(f"{output_prefix}_combined_legend.pdf", dpi=300, bbox_inches='tight')
    fig_combined_legend.savefig(f"{output_prefix}_combined_legend.png", dpi=300, bbox_inches='tight')
    plt.close(fig_combined_legend)

    # Generate all pairwise comparisons
    print(f"\nGenerating {len(methods) * (len(methods) - 1) // 2} pairwise comparisons...")

    for i in range(len(methods)):
        for j in range(i + 1, len(methods)):
            method1 = methods[i]
            method2 = methods[j]
            pair_prefix = f"{output_prefix}_{method1}_vs_{method2}"

            print(f"\n{'='*60}")
            print(f"Comparing {METHOD_DISPLAY_NAMES.get(method1, method1)} vs {METHOD_DISPLAY_NAMES.get(method2, method2)}")
            print(f"{'='*60}")

            plot_pairwise_comparison(
                method1, method2,
                all_method_data[method1], all_method_data[method2],
                METHOD_DISPLAY_NAMES, COLORS, pair_prefix
            )

    print("\n" + "="*60)
    print("All pairwise comparisons completed!")
    print("="*60)
