#!/usr/bin/env nextflow

// global variable
def current_time = new Date().format('yyyy-MM-dd_HH-mm-ss')


process GETBEDFROMGTF {
    input:
    path gtf

    output:
    path "${gtf.baseName}.bed", emit: out

    publishDir "${params.output_dir}/xref_files", mode: 'copy', overwrite: true

    script:
    """
    gxf2bed -i ${gtf} -o "${gtf.baseName}.bed"
    """
}

// split a bam to single-cell bam and internal bam w.r.t UB tag
process SPLIT_BAM_SC_AND_INTERNAL {
    input:
    tuple val(id), path(bam)

    output:
    path "*.bam", emit: out

    // publishDir "${params.output_dir}/xsplit-${current_time}", mode: 'copy', overwrite: true

    script:
    """
    #!/usr/bin/env bash

    samtools view -@8 -h ${bam} > ${id}.full.sam

    awk -v fileA="${id}.internal.sam" -v fileB="${id}.umied.sam" '
    {
        if (/^@/) { print \$0 > fileA; print \$0 > fileB; next }
        if (/UB:Z:[^ \t]/) { print \$0 > fileB; next }
        { print \$0 > fileA }
    }
    ' ${id}.full.sam

    samtools view -@8 -bS ${id}.internal.sam > ${id}.internal.bam
    samtools view -@8 -bS ${id}.umied.sam > ${id}.umied.bam
    rm ${id}.full.sam
    rm ${id}.internal.sam
    rm ${id}.umied.sam
    """
}

process READ_DIST{
    input:
    tuple path(bam), path(bed)

    output:
    path "${bam.baseName}.dist", emit: out

    publishDir "${params.output_dir}/xread_dist-${current_time}", mode: 'copy', overwrite: true

    script:
    """
    #!/usr/bin/env bash
    read_distribution.py -i ${bam} -r ${bed} > ${bam.baseName}.dist
    """
}

process INFEREXPERIMENT {
    input:
    tuple path(bam), path(bed)

    output:
    tuple path("${bam.baseName}.libtype.txt"), emit: out
    
    publishDir "${params.output_dir}/xlibtype-${current_time}", mode: 'copy', overwrite: true

    script:
    """
    #!/bin/bash
    infer_experiment.py -r ${bed} -i ${bam} > ${bam.baseName}.libtype.txt
    """
    // Example:
    // This is PairEnd Data
    // Fraction of reads failed to determine: 0.0172
    // Fraction of reads explained by "1++,1--,2+-,2-+": 0.4903
    // Fraction of reads explained by "1+-,1-+,2++,2--": 0.4925
}

workflow {
    println "Starting exploration workflow at ${current_time}"
    println "launched workflow at ${launchDir}"

    bams_ch = Channel.fromPath(params.bam_files)
        .map { bam -> 
            def sample_id = bam.baseName
            tuple(sample_id, bam)
        }

    bed = GETBEDFROMGTF(params.reference_gtf)

    bams_ch_split = SPLIT_BAM_SC_AND_INTERNAL(bams_ch)

    x = bams_ch_split
        .flatten()
        .combine(bed.out)
    
    x1 = x
    x2 = x
    
    READ_DIST(x1)
    INFEREXPERIMENT(x2)
}