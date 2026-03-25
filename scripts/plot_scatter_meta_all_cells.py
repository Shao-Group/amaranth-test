"""
Scatter plot of per-cell assembly performance for meta-assembly comparison.

Accepts pairs of <dir> <method> arguments specifying directories of GFFCompare
compareIC.stats files and their method labels. Plots precision vs matching
transcripts for each method.

Usage:
  python plot_scatter_meta_all_cells.py <dir1> <method1> [<dir2> <method2> ...] [--output-prefix all.scatter_meta]

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


def parse_compareIC_stats(filepath):
    """Parse compareIC.stats file and extract matching transcripts and precision."""
    matching_transcripts = 0
    precision = 0.0

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Matching intron chains:'):
                matching_transcripts = int(line.split(':')[1].strip())
            elif line.startswith('Intron chain level:'):
                parts = line.split('|')
                if len(parts) >= 2:
                    precision = float(parts[1].strip())

    return matching_transcripts, precision


# Predefined file patterns for known methods
PATTERNS = {
    'amaranth_meta': '*.stats',
    'aletsch': '*.stats',
    'psiclass': '*.stats',
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
    # Parse positional pairs + optional --output-prefix
    args_list = sys.argv[1:]
    output_prefix = 'all.scatter_meta'
    if '--output-prefix' in args_list:
        idx = args_list.index('--output-prefix')
        output_prefix = args_list[idx + 1]
        args_list = args_list[:idx] + args_list[idx+2:]

    if len(args_list) < 2 or len(args_list) % 2 != 0:
        print("Usage: python plot_scatter_meta_all_cells.py <dir1> <method1> [<dir2> <method2> ...] [--output-prefix <prefix>]")
        sys.exit(1)

    num_methods = len(args_list) // 2
    methods = []
    directories = []
    for i in range(num_methods):
        directories.append(args_list[2*i])
        methods.append(args_list[2*i + 1])

    method_configs = {}
    for i, method in enumerate(methods):
        pattern = PATTERNS.get(method, '*.compareIC.stats')
        method_configs[method] = {'dir': directories[i], 'pattern': pattern}

    fig, ax = plt.subplots(figsize=(5, 5))

    for method in list(methods)[::-1]:
        config = method_configs[method]
        matching_files = glob.glob(os.path.join(config['dir'], config['pattern']))

        filtered_files = [f for f in matching_files
                          if 'nr' not in os.path.basename(f) and 'meta' not in os.path.basename(f)]

        if not filtered_files:
            print(f"Warning: No files found for method {method} in {config['dir']}")
            continue

        correct_values = []
        precision_values = []
        for filepath in filtered_files:
            mt, prec = parse_compareIC_stats(filepath)
            correct_values.append(mt)
            precision_values.append(prec)

        print(f"\n{'='*60}\nMethod: {method}\n{'='*60}")
        print(f"Found {len(filtered_files)} files")
        print(f"Recall: {np.mean(correct_values):.2f} ± {np.std(correct_values):.2f}")
        print(f"Precision: {np.mean(precision_values):.2f} ± {np.std(precision_values):.2f}")

        display_name = METHOD_DISPLAY_NAMES.get(method, method)
        ax.scatter(correct_values, precision_values,
                   color=COLORS.get(method, '#000000'), s=10, zorder=2,
                   label=display_name, alpha=0.6)

    ax.set_xlabel('# matching transcripts', fontsize=14)
    ax.set_ylabel('precision (%)', fontsize=14)
    ax.tick_params(axis='x', labelsize=12, rotation=45)
    ax.tick_params(axis='y', labelsize=12, rotation=90)
    ax.set_ylim(0, 100)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    plt.savefig(f"{output_prefix}.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches='tight')
    plt.close(fig)

    # Combined legend
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
    fig_combined_legend.savefig(f"{output_prefix}_legend.pdf", dpi=300, bbox_inches='tight')
    fig_combined_legend.savefig(f"{output_prefix}_legend.png", dpi=300, bbox_inches='tight')
    plt.close(fig_combined_legend)
