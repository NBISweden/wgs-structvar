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
            file 'fastq.fq.gz' into fastqs

        publishDir 'results'

        module 'bioinfo-tools'
        module "$params.modules.samtools"

        """
        samtools bam2fq bamfile | gzip - > fastq.fq.gz
        """
    }
}
else {
    fastqs = Channel.fromPath(params.fastq)
}


if ( params.run_manta ) {
    process index_bamfile {
        input:
            file 'bamfile' from bamfile
        output:
            file 'bamfile.bai' into bamfile_index

        module 'bioinfo-tools'
        module "$params.modules.samtools"

        """
        samtools index bamfile
        """
    }

    process run_manta {
        input:
            file 'bamfile_tmp' from bamfile
            file 'bamfile.bai' from bamfile_index
        output:
            file 'sample.manta.sv.vcf.gz' into manta_vcf_gz
            file 'sample.manta.sv.vcf' into manta_vcf
            file 'sample.manta.sv.bed' into manta_bed

        publishDir 'results'

        module 'bioinfo-tools'
        module "$params.modules.manta"

        """
        # Manta follows symlinks and expects the index to be with the original
        # file, so we copy it, and then clean up at script EXIT with a trap.
        # TODO, this is fixed in manta v0.29.5, https://github.com/Illumina/manta/issues/32
        DIR=`pwd`
        function cleanup() {
            cd \$DIR
            if [ -f bamfile ]; then
                rm bamfile
            fi
        }
        trap cleanup EXIT
        cp bamfile_tmp bamfile

        configManta.py --normalBam bamfile --referenceFasta $params.ref_fasta --runDir testRun
        cd testRun
        ./runWorkflow.py -m local -j $params.threads
        mv results/variants/diploidSV.vcf.gz ../sample.manta.sv.vcf.gz
        cd ..
        gunzip -c sample.manta.sv.vcf.gz > sample.manta.sv.vcf
        $params.programs.svvcf2bed sample.manta.sv.vcf > sample.manta.sv.bed
        """
    }
}


def usage_message() {
    log.info ''
    log.info 'Usage:'
    log.info '    nextflow main.nf --bam <bamfile> [more options]'
    log.info ''
    log.info 'Options:'
    log.info '    --help          Show this message and exit'
    log.info '    --bam           Input bamfile'
    log.info '    --fastq         Input fastqfile (default is bam but with fq as fileending)'
    log.info '    --run_manta     Run manta'
    log.info ''
}

def infer_fastq_from_bam() {
    path = params.bam.replaceAll(/.bam$/, '.fq.gz')
    if (file(path).exists()) {
        return path
    }
    return false
}
