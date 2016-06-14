#!/usr/bin/env nextflow

/*
WGS Structural Variations Pipeline
*/


// 0. Pre-flight checks

if (params.help) {
    usage_message()
    exit 0
}

if (!params.bam) {
    exit 1, 'You need to specify a bam file, see --help for more information'
}

bamfile = file(params.bam)

if (! bamfile.exists()) {
    exit 1, "The bamfile, '$params.bam', does not exist"
}

if (!params.fastq) {
    params.fastq = infer_fastq_from_bam()
}

if (!params.fastq) {
    process create_fastq {
        input:
            file 'bamfile' from bamfile

        output:
            file 'fastq.gz' into fastqs

        publishDir 'results'

        module 'bioinfo-tools'
        module "$params.modules.samtools"

        """
        samtools bam2fq bamfile | gzip - > fastq.gz
        """
    }
}
else {
    fastqs = Channel.fromPath(params.fastq)
}


def usage_message() {
    log.info ''
    log.info 'Usage:'
    log.info '    nextflow main.nf --bam <bamfile> [more options]'
    log.info ''
    log.info 'Options:'
    log.info '    --help     Show this message and exit'
    log.info '    --bam      Input bamfile'
    log.info '    --fastq    Input fastqfile (default is bam but with fq as fileending)'
    log.info ''
}

def infer_fastq_from_bam() {
    path = params.bam.replaceAll(/.bam$/, '.fq')
    if (file(path).exists()) {
        return path
    }
    return false
}
