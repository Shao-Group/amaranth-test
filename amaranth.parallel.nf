#!/usr/bin/env nextflow

// global variable
def current_time = new Date().format('yyyy-MM-dd_HH-mm-ss')

process RUNAMARANTH {
    input:
    tuple val(id), path(bam)

    output:
    path "${id}.amaranth.gtf", emit: gtf

    publishDir "${params.output_dir}/amaranth-${current_time}", mode: 'copy', overwrite: true

    script:

    """
    amaranth -i ${bam} \
        --gene_name_prefix cell_${id}. \
        -o ${id}.amaranth
    """
}



workflow {
    println "Starting workflow at ${current_time}"

    bams_ch = Channel.fromPath(params.bam_files)
        .map { bam -> 
            def sample_id = bam.baseName
            tuple(sample_id, bam)
        }

    amr = RUNAMARANTH(bams_ch)
}