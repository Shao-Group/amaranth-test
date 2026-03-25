"""
Controlled comparison of assembler precision across cells.

For each cell, reads ROC files (.roc) and computes precision controlled at the
lowest sensitivity level (correct transcript count) across all methods, enabling
fair comparison at a common operating point.

Usage:
  python plot_controlled_all_cells.py --roc-dir <dir_containing_roc_files> [--output-prefix controlled_all_cells_0]
"""
import argparse
import matplotlib.pyplot as plt
import numpy as np
import sys
import os


def get_methods_files(methods_list, dir):
    methods_files = {}
    for method in methods_list:
        method_files = []
        for root, _, files in os.walk(dir):
            for file in files:
                if file.endswith(method + ".roc"):
                    method_files.append(os.path.join(root, file))
        method_files.sort()
        methods_files[method] = method_files

    file_number = -1
    for k, v in methods_files.items():
        if file_number == -1:
            file_number = len(v)
        else:
            assert file_number == len(v), f"Method {k} has {len(v)} files, expected {file_number}"
    return methods_files

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

def controlled_1_cell(axis, methods, simple=False):
    lowest_precision = min(max(x[axis]) for x in methods)
    controlled_dependent_values = []
    for method in methods:
        k, v = -1, -1
        for row in zip(*method):
            if row[axis] >= lowest_precision:
                v = row[1-axis]
                k = row[axis]
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
    assert len(controlled_dependent_values) == len(methods)
    return controlled_dependent_values


def plot_scatter(axis, method_files, amr_key, output_prefix, simple=False):
    plt.figure(figsize=(5, 5))

    methods = []
    files = []
    for method, f in method_files.items():
        methods.append(method)
        files.append(f)

    colors = ['#a89c87', '#468a71', '#FF9500']

    all_controlled_dependent_values = []
    for i in range(len(files[0])):
        controlled_dependent_values = controlled_1_cell(axis, [read_roc_file(x[i]) for x in files], simple)
        all_controlled_dependent_values.append(controlled_dependent_values)
    all_controlled_dependent_values.sort(key=lambda x: x[methods.index(amr_key + ".covTX")])

    display_names = ['Amaranth', 'Scallop2', 'StringTie2']
    desired_order = [amr_key + '.covTX', 'sc2.covTX', 'stg.covTX']

    method_to_display = dict(zip(desired_order, display_names))
    method_to_color = {method: colors[methods.index(method)] for method in methods}

    for idx, method in enumerate(desired_order):
        j = methods.index(method)
        y = [x[j] for x in all_controlled_dependent_values]
        x_vals = list(range(len(y)))
        plt.scatter(x_vals, y, s=10, zorder=2, c=[method_to_color[method]], label=method_to_display[method], alpha=0.6)
        avg = np.mean(y)
        print(f"avg {method}: {avg:.2f}")

    plt.xlabel(f'cells (n={len(files[0])})', fontsize=14)
    plt.ylabel('adj precision (%)' if axis == 0 else 'adj # matching transcripts', fontsize=14)
    if axis == 0:
        plt.ylim(0, 100)
    plt.xticks(fontsize=12, rotation=45, ha='right')
    plt.yticks(fontsize=12, rotation=90, va='center')
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    plt.savefig(f'{output_prefix}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}.pdf', dpi=300, bbox_inches='tight')

    handles, labels = plt.gca().get_legend_handles_labels()
    fig_legend = plt.figure(figsize=(3, 0.5))
    fig_legend.legend(handles, labels, loc='center', ncol=3, fontsize=12, frameon=True, fancybox=True, markerscale=2)
    fig_legend.savefig(f'{output_prefix}_legend.png', dpi=300, bbox_inches='tight')
    fig_legend.savefig(f'{output_prefix}_legend.pdf', dpi=300, bbox_inches='tight')
    plt.close(fig_legend)
    plt.clf()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--roc-dir', required=True,
                        help='Directory containing .roc files from GFFCompare evaluation')
    parser.add_argument('--output-prefix', default='controlled_all_cells_0',
                        help='Output file prefix (default: controlled_all_cells_0)')
    args = parser.parse_args()

    amr = "amr"
    sc2 = "sc2"
    stg = "stg"
    methods_list = [x + ".covTX" for x in [stg, sc2, amr]]
    methods_files = get_methods_files(methods_list, args.roc_dir)

    plot_scatter(0, methods_files, amr, args.output_prefix)
