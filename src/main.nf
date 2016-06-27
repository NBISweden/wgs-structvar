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

if (params.run_all) {
    params.run_fermikit = true
    params.run_manta = true
    params.run_intersections = true
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


if (params.run_fermikit) {
    process fermikit_calling {
        input:
            file 'sample.fq.gz' from fastqs
        output:
            file 'sample.fermikit.sv.vcf.gz' into fermi_vcf_gz
            file 'sample.fermikit.sv.vcf' into fermi_vcf
            file 'sample.fermikit.sv.bed' into fermi_bed

        publishDir 'results'

        module 'bioinfo-tools'
        module "$params.modules.fermikit"
        module "$params.modules.samtools"
        module "$params.modules.vcftools"
        module "$params.modules.tabix"

        """
        fermi2.pl unitig -s$params.genome_size -t$params.threads -l$params.readlen -p sample sample.fq.gz > sample.mak
        make -f sample.mak
        run-calling -t$params.threads $params.ref_fasta sample.mag.gz > calling.sh
        bash calling.sh
        vcf-sort -c sample.sv.vcf.gz > sample.fermikit.sv.vcf
        bgzip -c sample.fermikit.sv.vcf > sample.fermikit.sv.vcf.gz
        $params.programs.svvcf2bed sample.fermikit.sv.vcf > sample.fermikit.sv.bed
        """
    }
}

if (params.run_intersections) {
    process download_masks {
        output:
            file 'ceph18.b37.lumpy.exclude.2014-01-15.bed' into mask_ceph
            file 'LCR-hs37d5.bed.gz' into mask_lcr

        // We only need one core for this part
        clusterOptions = {
            "-A $params.project -p core"
        }

        """
        wget https://github.com/cc2qe/speedseq/raw/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed
        wget https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs37d5.bed.gz
        """
    }

    process run_intersections {
        input:
            file 'manta.bed' from manta_bed
            file 'fermi.bed' from fermi_bed
            file 'mask_ceph.bed' from mask_ceph
            file 'mask_lcr.bed' from mask_lcr
        output:
            file 'manta_masked*.bed'
            file 'fermi_masked*.bed'
            file 'combined*.bed'

        publishDir 'results'

        // We only need one core for this part
        clusterOptions = {
            '-A $params.project -p core'
        }

        module 'bioinfo-tools'
        module "$params.modules.bedtools"

        """
        cat manta.bed \
            | bedtools intersect -v -a stdin -b mask_lcr.bed -f 0.25 \
            | bedtools intersect -v -a stdin -b mask_ceph.bed -f 0.25 > manta_masked.bed
        cat fermi.bed \
            | bedtools intersect -v -a stdin -b mask_lcr.bed -f 0.25 \
            | bedtools intersect -v -a stdin -b mask_ceph.bed -f 0.25 > fermi_masked.bed

        ## In case grep doesn't find anything we want to continue
        set +e

        ## Create filtered bed files
        for FILE in manta_masked fermi_masked; do
            for WORD in DEL INS DUP; do
                grep -w \$WORD \$FILE.bed > \$FILE.\$WORD.bed
            done
        done

        ## But now we want to abort in case of errors
        set -e

        ## Create intersected bed file
        for WORD in DEL INS DUP; do
            intersectBed -a fermi_masked.\$WORD.bed -b manta_masked.DEL.bed -f 0.5 -r | sort -k1,1V -k2,2n > combined_masked.SV.\$WORD.bed
        done
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
    log.info '    --run_fermikit  Run fermikit'
    log.info '    --run_intersections  Run intersections'
    log.info '    --run_all       Run all'
    log.info ''
}

def infer_fastq_from_bam() {
    path = params.bam.replaceAll(/.bam$/, '.fq.gz')
    if (file(path).exists()) {
        return path
    }
    return false
}
