#!/usr/bin/env nextflow

/*
WGS Structural Variation Pipeline
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


if (!params.project) {
    exit 1, 'You need to specify what project to run under, see --help for more information'
}


workflowSteps = processWorkflowSteps(params.steps)


startup_message()

// 1. Run manta

// Try to guess location of bamindex file. If we can't find it create it
// else put that in the bamfile_index channel.

bamindex = infer_bam_index_from_bam()
if (!bamindex) {
    process index_bamfile {
        input:
            file 'bamfile' from bamfile
        output:
            file 'bamfile.bai' into bamfile_index

        executor choose_executor()
        queue 'core'
        time params.runtime.simple

        module 'bioinfo-tools'
        module "$params.modules.samtools"

        when: 'indexbam' in workflowSteps

        script:
        """
        samtools index bamfile
        """
    }
}
else {
    // The bamfile file already exists, put it in the channel.
    Channel.fromPath( bamindex ).set { bamfile_index }
}

process manta {
    input:
        file 'bamfile' from bamfile
        file 'bamfile.bai' from bamfile_index
    output:
        file 'manta.vcf' into manta_vcf

    publishDir params.outdir, mode: 'copy'

    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    time { params.runtime.caller * 2**(task.attempt-1) }
    maxRetries 3
    queue 'core'
    cpus 4

    module 'bioinfo-tools'
    module "$params.modules.manta"

    when: 'manta' in workflowSteps

    script:
    """
    configManta.py --normalBam bamfile --referenceFasta $params.ref_fasta --runDir testRun
    cd testRun
    ./runWorkflow.py -m local -j $params.threads
    mv results/variants/diploidSV.vcf.gz ../manta.vcf.gz
    cd ..
    gunzip -c manta.vcf.gz > manta.vcf
    """
}


// 2. Run fermikit

// Try to guess location of fastq file. If we can't find it create it
// else put that in the fastq channel.
if (!params.fastq) {
    params.fastq = infer_fastq_from_bam()
}

if (!params.fastq) {
    process create_fastq {
        input:
            file 'bamfile' from bamfile
        output:
            file 'fastq.fq.gz' into fastq

        executor choose_executor()
        queue 'core'
        time params.runtime.simple

        module 'bioinfo-tools'
        module "$params.modules.samtools"

        when: 'fastq' in workflowSteps

        script:
        """
        samtools bam2fq bamfile | gzip - > fastq.fq.gz
        """
    }
}
else {
    // The fastq file already exists, put it in the channel.
    Channel.fromPath( params.fastq ).set { fastq }
}

process fermikit {
    input:
        file 'sample.fq.gz' from fastq
    output:
        file 'fermikit.vcf' into fermi_vcf

    publishDir params.outdir, mode: 'copy'

    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    time { params.runtime.fermikit * 2**( task.attempt - 1 ) }
    maxRetries 3
    queue 'node'

    module 'bioinfo-tools'
    module "$params.modules.fermikit"
    module "$params.modules.samtools"
    module "$params.modules.vcftools"
    module "$params.modules.tabix"

    when: 'fermikit' in workflowSteps

    script:
    """
    fermi2.pl unitig -s$params.genome_size -t$params.threads -l$params.readlen -p sample sample.fq.gz > sample.mak
    make -f sample.mak
    run-calling -t$params.threads $params.ref_fasta sample.mag.gz > calling.sh
    bash calling.sh
    vcf-sort -c sample.sv.vcf.gz > fermikit.vcf
    bgzip -c fermikit.vcf > fermikit.vcf.gz
    """
}


// 3. Create summary files

// Collect vcfs and beds into one channel
vcfs = manta_vcf.mix( fermi_vcf )

mask_files = [
    "$baseDir/data/ceph18.b37.lumpy.exclude.2014-01-15.bed",
    "$baseDir/data/LCR-hs37d5.bed.gz"
]

masks = mask_files.collect { file(it) }.channel()
// Collect both bed files and combine them with the mask files
vcfs.tap { vcfs }.spread( masks.buffer(size: 2) ).set { mask_input }

process mask_beds {
    input:
        set file(svfile), file(mask1), file(mask2) from mask_input
    output:
        file '*_masked.vcf' into masked_vcfs

    executor choose_executor()
    queue 'core'
    time params.runtime.simple

    module 'bioinfo-tools'
    module "$params.modules.bedtools"

    """
    BNAME=\$( echo $svfile | cut -d. -f1 )
    MASK_FILE=\${BNAME}_masked.vcf
    cat $svfile \
        | bedtools intersect -header -v -a stdin -b $mask1 -f 0.25 \
        | bedtools intersect -header -v -a stdin -b $mask2 -f 0.25 > \$MASK_FILE
    """
}


// To make intersect files we need to combine them into one channel with
// toSortedList() (fermi is before manta in alphabet). And also figure out if we
// have one or two files, therefore the tap and count_vcfs.
masked_vcfs.tap { count_vcfs_tmp }
           .tap { masked_vcfs }
           .toSortedList().set { intersect_input }
count_vcfs_tmp.count().set { count_vcfs }

process intersect_files {
    input:
        set file(fermi_vcf), file(manta_vcf) from intersect_input
        val nvcfs from count_vcfs
    output:
        file "combined_masked.vcf" into intersections

    executor choose_executor()
    queue 'core'
    time params.runtime.simple

    module 'bioinfo-tools'
    module "$params.modules.bedtools"

    when: nvcfs == 2

    script:
    """
    ## Create intersected vcf files
    for WORD in DEL INS DUP; do
        intersectBed -a <( grep -w "^#.\\+\\|\$WORD" $fermi_vcf) \
                     -b <( grep -w "^#.\\+\\|\$WORD" $manta_vcf) \
            -f 0.5 -r \
            | sort -k1,1V -k2,2n > combined_masked_\${WORD,,}.vcf
    done

    cat <( grep -v -w '^#.\\+\\|DEL\\|INS\\|DUP' $fermi_vcf ) \
        <( grep -v -w '^#.\\+\\|DEL\\|INS\\|DUP' $manta_vcf ) \
        | cut -f 1-8 \
        | sort -k1,1V -k2,2n > combined_masked_OTHER.vcf

    sort -k1,1V -k2,2n combined_masked_*.vcf >> combined_masked.vcf
    """
}

annotate_files = intersections.flatten().mix( masked_vcfs.tap { masked_vcfs } )

process variant_effect_predictor {
    input:
        file infile from annotate_files.tap { annotate_files }
    output:
        file '*.vep.vcf'

    publishDir params.outdir, mode: 'copy'

    executor choose_executor()
    queue 'core'
    time params.runtime.simple

    module 'bioinfo-tools'
    module "$params.modules.vep"

    when: 'vep' in workflowSteps

    script:
    """
    INFILE="$infile"
    OUTFILE="\${INFILE%.vcf}.vep.vcf"
    VEP_CACHE="/sw/data/uppnex/vep/84"
    ASSEMBLY="GRCh37"

    case "\$INFILE" in
        *vcf) FORMAT="vcf" ;;
        *bed) FORMAT="ensembl" ;;
        *)    printf "Unrecognized format for '%s'" "\$INFILE" >&2
              exit 1;;
    esac

    ## If the input file is empty, just copy it
    if [[ -f "\$INFILE" && -s "\$INFILE" ]]; then
        cp "\$INFILE" "\$OUTFILE"
        exit
    fi

    variant_effect_predictor.pl \
        -i "\$infile"              \
        --format "\$FORMAT"        \
        -cache --dir "\$VEP_CACHE" \
        -o "\$outfile"             \
        --vcf                      \
        --merged                   \
        --regulatory               \
        --force_overwrite          \
        --sift b                   \
        --polyphen b               \
        --symbol                   \
        --numbers                  \
        --biotype                  \
        --total_length             \
        --canonical                \
        --ccds                     \
        --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \
        --assembly "\$ASSEMBLY" \
        --offline
    """
}

process snpEff {
    input:
        file infile from annotate_files.tap { annotate_files }
    output:
        file '*.snpeff.vcf'

    publishDir params.outdir, mode: 'copy'

    executor choose_executor()
    queue 'core'
    time params.runtime.simple

    module 'bioinfo-tools'
    module "$params.modules.snpeff"

    when: 'snpeff' in workflowSteps

    script:
    """
    INFILE="$infile" ## Use bash-semantics for variables
    OUTFILE="\${INFILE%.vcf}.snpeff.vcf"
    SNPEFFJAR=''

    for P in \$( tr ':' ' ' <<<"\$CLASSPATH" ); do
        if [ -f "\$P/snpEff.jar" ]; then
            SNPEFFJAR="\$P/snpEff.jar"
            break
        fi
    done
    if [ -z "\$SNPEFFJAR" ]; then
        printf "Can't find snpEff.jar in '%s'" "\$CLASSPATH" >&2
        exit 1
    fi

    sed 's/ID=AD,Number=./ID=AD,Number=R/' "\$INFILE" \
        | vt decompose -s - \
        | vt normalize -r $params.ref_fasta - \
        | java -Xmx7G -jar "\$SNPEFFJAR" -formatEff -classic GRCh37.75 \
        > "\$OUTFILE"
    """
}


// Utility functions

def usage_message() {
    log.info ''
    log.info 'Usage:'
    log.info '    nextflow main.nf --bam <bamfile> [more options]'
    log.info ''
    log.info 'Options:'
    log.info '  Required'
    log.info '    --bam           Input bamfile'
    log.info '    --project       Uppmax project to log cluster time to'
    log.info '  Optional'
    log.info '    --help          Show this message and exit'
    log.info '    --fastq         Input fastqfile (default is bam but with fq as fileending)'
    log.info '    --steps         Specify what steps to run, comma separated:'
    log.info '                Callers: manta, fermikit, cnvnator (choose one or many)'
    log.info '                Annotation: vep OR snpeff'
    log.info '    --outdir        Directory where resultfiles are stored'
    log.info ''
}

def startup_message() {
    revision = grab_git_revision()

    log.info "======================"
    log.info "WGS-structvar pipeline"
    log.info "======================"
    log.info "Bamfile    : $params.bam"
    log.info "Scriptdir  : $baseDir"
    log.info "Revision   : $revision"
    log.info "Work dir   : $workDir"
    log.info "Output dir : $params.outdir"
    log.info "Project    : $params.project"
    log.info "Will run   : " + workflowSteps.join(", ")
    log.info ""
}

def grab_git_revision() {
    if ( workflow.commitId ) { // it's run directly from github
        return workflow.commitId
    }

    // Try to find the revision directly from git
    head_pointer_file = file("${baseDir}/.git/HEAD")
    if ( ! head_pointer_file.exists() ) {
        return ''
    }
    ref = head_pointer_file.newReader().readLine().tokenize()[1]

    ref_file = file("${baseDir}/.git/$ref")
    if ( ! ref_file.exists() ) {
        return ''
    }
    revision = ref_file.newReader().readLine()

    return revision
}

def infer_bam_index_from_bam() {
    // If the ".bam.bai" file does not exist, try ".bai" without ".bam"
    return infer_filepath(params.bam, /$/, '.bai')
        ?: infer_filepath(params.bam, /.bam$/, '.bai')
}

def infer_fastq_from_bam() {
    return infer_filepath(params.bam, /.bam$/, '.fq.gz')
}

def infer_filepath(from, match, replace) {
    path = file( from.replaceAll(match, replace) )
    if (path.exists()) {
        return path
    }
    return false
}

def nextflow_running_as_slurmjob() {
    if ( System.getenv()["SLURM_JOB_ID"] ) {
        return true
    }
    return false
}

/* If the nextflow deamon is running as a slurm job, we can use the local CPU
 * for a lot of our work */
def choose_executor() {
    return nextflow_running_as_slurmjob() ? 'local' : 'slurm'
}

def processWorkflowSteps(steps) {
    if ( ! steps ) {
        return []
    }

    workflowSteps = steps.split(',').collect { it.trim().toLowerCase() }

    if ('vep' in workflowSteps && 'snpeff' in workflowSteps) {
        exit 1, 'You can only run one annotator, either "vep" or "snpeff"'
    }

    if ('manta' in workflowSteps) {
        workflowSteps.push( 'indexbam' )
    }

    if ('fermikit' in workflowSteps) {
        workflowSteps.push( 'fastq' )
    }

    return workflowSteps
}
