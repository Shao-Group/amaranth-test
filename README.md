# Amaranth-test

This repository contains the code and instructions for reproducing experimental results for [Amaranth](https://github.com/Shao-Group/amaranth), a single-cell RNA transcript assembler. It benchmarks Amaranth against StringTie2 and Scallop2 using Nextflow pipelines, evaluating assembly quality with GFFCompare.

## Installation

Most dependencies can be installed using [Pixi](https://pixi.sh) and then activate the environment:

```bash
pixi install
pixi shell
```
Additionally, [rnaseqtools](https://github.com/Shao-Group/rnaseqtools) need to be installed according to the instruction in the repo.

## Quick Start

1. Copy the template config and edit the paths:

```bash
cp params.template.config params.config
# Edit params.config: set bam_files and reference_gtf
```

2. Run a workflow:

```bash
nextflow run main.nf -c params.config -resume
```

## Workflows

| Workflow | Purpose |
|----------|---------|
| `main.nf` | Single-cell benchmarking: runs Amaranth, StringTie2, and Scallop2 on per-cell BAMs, evaluates with GFFCompare|
| `main-meta.nf` | Meta-assembly benchmarking: runs Amaranth in `--meta` mode |
| `main-ablation.nf` | Ablation study: systematically varies one Amaranth parameter at a time (permissive-baseline isolation) |
| `explore.nf` | Data exploration: library type detection, read distribution, BAM splitting by UMI tag |

## Configuration

All workflows read parameters from a Nextflow config file. See `params.template.config` for the full list. Required parameters:

- `bam_files` — glob pattern for aligned BAM files (one per cell)
- `reference_gtf` — reference gene annotation in GTF format

Optional:

- `output_dir` — results directory (default: `results`)
- `num_ref_transcripts` — for ROC curves; auto-detected from GTF if omitted

## Datasets

All data preprocessing follows the instructions provided in the [scallop2-test](https://github.com/Shao-Group/scallop2-test) repository.

**HEK293T**: 192 human cells from the Smart-seq3 project (Hagemann-Jensen et al., 2020). Sequenced with strand-specific, paired-end protocol using barcoding technology. Raw data are demultiplexed and preprocessed using zUMIs (with STAR for alignment). Data is available at ArrayExpress [E-MTAB-8735](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8735/).

**Mouse-Fibroblast**: 369 mouse fibroblast cells from the Smart-seq3 project (Hagemann-Jensen et al., 2020). Same sequencing and preprocessing protocol. Data is available at ArrayExpress [E-MTAB-8735](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8735/).
