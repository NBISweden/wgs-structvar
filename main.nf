#!/usr/bin/env nextflow
import static java.util.UUID.randomUUID

/*
WGS Structural Variation Pipeline
*/

// 0. Pre-flight checks

if (params.help) {
    usage_message()
    exit 0
}

check_input_params()


/* Figure out what steps to run */
workflowSteps = processWorkflowSteps(params.steps)


/*** Setup input channels
 *  We create an array that will contain the file to work on, a unique
 * identifier and the output directory for this file.
 *  Last two elements of the array is ALWAYS `uuid` and `outdir` in that
 * order.
 */

if (params.bam) {
    ch_in = setup_input_channel_from_bam(params.bam)
} else if (params.runfile) {
    ch_in = setup_input_channel_from_runfile(params.runfile)
}

// One input for fermi and one for manta
ch_in_manta = ch_in.tap { ch_in_fermi }

/* Display startup message and go, go, go! */
startup_message()

// 1. Run manta

// Add bamindex path to the channel and generate indexes if they're missing
ch_bamfile = Channel.create()
ch_index_bam = Channel.create()
ch_already_indexed = Channel.create()

ch_in_manta.map {
        index = infer_bam_index_from_bam(it[0])
        [ it[0], file(index), it[1], it[2] ]
    }.choice(ch_bamfile, ch_already_indexed) { it[1].exists() ? 1 : 0 }

// Remove bam index filename from channel map so it doesnt get produced in the working dir
// by the process index_bamfile
ch_bamfile.map {
        [ it[0], it[2], it[3] ]
    }.into(ch_index_bam)

process index_bamfile {
    input:
        set file(bamfile), val(uuid), val(dir) from ch_index_bam
    output:
        set file(bamfile), file('*.bam.bai'), val(uuid), val(dir) into ch_indexed_bam

    tag "$uuid"

    executor choose_executor()

    when: 'manta' in workflowSteps

    script:
    """
    samtools index "$bamfile"
    """
}

ch_already_indexed.mix( ch_indexed_bam ).set { ch_manta }

process manta {
    input:
        set file(bamfile), file(bamindex), val(uuid), val(dir) from ch_manta
    output:
        set file("manta.vcf"), val(uuid), val(dir) into ch_manta_vcf

    tag "$uuid"
    publishDir "$dir", mode: 'copy', saveAs: { "$params.prefix$it" }

    when: 'manta' in workflowSteps

    script:
    """
    if [ -z "\$SLURM_CPUS_ON_NODE" ]; then
        CPUS=1
    else
        CPUS=\$SLURM_CPUS_ON_NODE
    fi

    echo "CPUS: \$CPUS"
    configManta.py --normalBam $bamfile --referenceFasta $params.ref_fasta --runDir testRun
    cd testRun
    ./runWorkflow.py -m local -j \$CPUS
    gunzip -c results/variants/diploidSV.vcf.gz > ../manta.vcf
    """
}

// 2. Run fermikit
ch_bamfile_2fastq = Channel.create()
ch_create_fastq = Channel.create()
ch_has_fastq = Channel.create()

ch_in_fermi.map {
        fastqfile = infer_fastq_from_bam(it[0])
        [ it[0], file(fastqfile), it[1], it[2] ]
    }.choice( ch_bamfile_2fastq, ch_has_fastq ) { it[1].exists() ? 1 : 0 }

// Remove fastq file from input so it is not produced in current (launch) dir
// by create_fastq, but instead it will be in a workdir
ch_bamfile_2fastq.map {
        [ it[0], it[2], it[3] ]
    }.into(ch_create_fastq)

process create_fastq {
    input:
        set file(bamfile), val(uuid), val(dir) from ch_create_fastq
    output:
        set file('fastq.fq.gz'), val(uuid), val(dir) into ch_created_fastq

    tag "$uuid"

    executor choose_executor()

    when: 'fermikit' in workflowSteps

    script:
    """
    samtools bam2fq "$bamfile" | gzip - > fastq.fq.gz
    """
}

ch_has_fastq.map { [it[1], it[2], it[3]] }.mix( ch_created_fastq ).set { ch_fermikit }

process fermikit {
    input:
        set file('sample.fq.gz'), val(uuid), val(dir) from ch_fermikit
    output:
        set file('*.vcf'), val(uuid), val(dir) into ch_fermi_vcf

    tag "$uuid"
    publishDir "$dir", mode: 'copy', saveAs: { "$params.prefix$it" }

    when: 'fermikit' in workflowSteps

    script:
    """
    if [ -z "\$SLURM_CPUS_ON_NODE" ]; then
        CPUS=1
    else
        CPUS=\$SLURM_CPUS_ON_NODE
    fi

    fermi2.pl unitig -s3g -t\$CPUS -l150 -p sample sample.fq.gz > sample.mak
    make -f sample.mak
    run-calling -t\$CPUS $params.ref_fasta sample.mag.gz > calling.sh
    bash calling.sh
    vcf-sort -c sample.sv.vcf.gz > fermikit.vcf
    bgzip -c fermikit.vcf > fermikit.vcf.gz
    """
}

// 3. Create summary files

ch_vcfs = ch_manta_vcf.mix( ch_fermi_vcf )

process artifact_mask_vcfs {
    input:
        set file(svfile), val(uuid), val(dir) from ch_vcfs
    output:
        set file(svfile), val(uuid), val(dir) into ch_artifact_masked_vcfs

    tag "$uuid $svfile"

    executor choose_executor()

    """
    BNAME=\$( echo $svfile | cut -d. -f1 )
    MASK_DIR=$params.mask_artifacts_dir

    # We don't want to change the filename in this process so we copy the
    # infile and remove the symbolic link. And then recreate the file at the
    # end.
    cp "$svfile" workfile
    rm "$svfile" # It's a link, should be ok :)
    for mask in \$MASK_DIR/*; do
        if [[ ! -f "\$mask" ]]; then
            continue
        fi
        cat workfile \
            | bedtools intersect -header -v -a stdin -b "\$mask" -f 0.25 \
            > tempfile
        mv tempfile workfile
    done
    mv workfile "$svfile"
    """
}


ch_nocohort_mask = ch_artifact_masked_vcfs.tap { ch_cohort_mask_in }
reciprocal = params.no_sg_reciprocal ? '': '-r'

process cohort_mask_vcfs {
    input:
        set file(svfile), val(uuid), val(dir) from ch_cohort_mask_in
    output:
        set file('*_cohort_masked.vcf'), val(uuid), val(dir) into ch_cohort_masked_vcfs

    tag "$uuid $svfile"

    executor choose_executor()
    when 'mask_cohort' in workflowSteps

    """
    BNAME=\$( echo $svfile | cut -d. -f1 )
    MASK_FILE=\${BNAME}_cohort_masked.vcf
    MASK_DIR=$params.mask_cohort_dir

    cp $svfile workfile
    for mask in \$MASK_DIR/*; do
        if [[ ! -f "\$mask" ]]; then
            continue
        fi
        cat workfile \
            | bedtools intersect -header -v -a stdin -b "\$mask" -f "$params.sg_mask_ovlp" "$reciprocal" \
            > tempfile
        mv tempfile workfile
    done
    mv workfile "\$MASK_FILE"
    """
}

ch_masked_vcfs = ch_nocohort_mask.mix(ch_cohort_masked_vcfs)

// To make intersect files we need to combine them into one channel with
// toList() and then sort in the map so that fermi is before manta in the
// channel. We can't use toSortedList here since the full pathname is used
// which includes the work directory (`3a/7c63f4...`).
ch_masked_vcfs.tap { ch_masked_vcfs_vep }
              .groupTuple(by: 1, size: 2)
              .set { ch_intersect_input }

process intersect_files {
    input:
        set file(vcfs), val(uuid), val(dir) from ch_intersect_input
    output:
        set file('combined_*.vcf'), val(uuid), val("${dir[0]}") into ch_intersections

    tag "$uuid"

    executor choose_executor()

    when: 'make_intersect' in workflowSteps

    script:
    """
    if head -n 10 ${vcfs[0]} | grep -q 'source=htsbox'; then
        fermi_vcf=${vcfs[0]}
        manta_vcf=${vcfs[1]}
    else
        fermi_vcf=${vcfs[1]}
        manta_vcf=${vcfs[0]}
    fi

    OUTNAME=`basename \$fermi_vcf|sed 's/fermikit_//;s/.vcf//'`
    ## Create intersected vcf files
    for WORD in DEL INS DUP; do
        intersectBed -a <( grep -w "^#.\\+\\|\$WORD" \$fermi_vcf) \
                     -b <( grep -w "^#.\\+\\|\$WORD" \$manta_vcf) \
            -f 0.5 -r \
            | sort -k1,1V -k2,2n > cmb_\${OUTNAME}_\${WORD,,}.vcf
    done

    cat <( grep -v -w '^#.\\+\\|DEL\\|INS\\|DUP' \$fermi_vcf ) \
        <( grep -v -w '^#.\\+\\|DEL\\|INS\\|DUP' \$manta_vcf ) \
        | cut -f 1-8 \
        | sort -k1,1V -k2,2n > cmb_\${OUTNAME}_OTHER.vcf

    ( grep '^#' \$fermi_vcf; \
        sort -k1,1V -k2,2n cmb_\${OUTNAME}_*.vcf ) >> combined_\${OUTNAME}.vcf
    """
}

if ( 'normalize' in workflowSteps ) {
    ch_masked_vcfs_vep.mix( ch_intersections ).set { ch_normalize_vcf }
}
else {
    ch_masked_vcfs_vep.mix( ch_intersections ).set { ch_annotate }

    // So we don't get stuck in an infinite loop
    Channel.empty().set { ch_normalize_vcf }
}

process normalize_vcf {
    input:
        set file(infile), val(uuid), val(dir) from ch_normalize_vcf
    output:
        set file("*.vt.vcf"), val(uuid), val(dir) into ch_normalized_vcf

    tag "$uuid - $infile"

    executor choose_executor()

    """
    INFILE="$infile"
    OUTFILE="\${INFILE%.vcf}.vt.vcf"

    ## If the input file is empty, just copy it
    if [[ -f "\$INFILE" && -s "\$INFILE" ]]; then
        cp "\$INFILE" "\$OUTFILE"
        exit
    fi


    ## Normalization
    sed 's/ID=AD,Number=./ID=AD,Number=R/' "\$INFILE" \
        | vt decompose -s - > "\$INFILE.vt_temp"

    if ! vt normalize -r "$params.ref_fasta" "\$INFILE.vt_temp" > "\$OUTFILE"
    then
        printf "VT normalisation Failed\n" >&2
        cp "\$INFILE.vt_temp" "\$OUTFILE"
    fi
    """
}

if ( 'normalize' in workflowSteps ) {
    ch_annotate_snpeff = ch_normalized_vcf.tap { ch_annotate_vep }
}
else {
    ch_annotate_snpeff = ch_annotate.tap { ch_annotate_vep }
}

process variant_effect_predictor {
    input:
        set file(infile), val(uuid), val(dir) from ch_annotate_vep
    output:
        file '*.vep.vcf'

    tag "$uuid - $infile"
    publishDir "$dir", mode: 'copy', saveAs: { "$params.prefix$it" }

    when: 'vep' in workflowSteps

    script:
    """
    INFILE="$infile"
    OUTFILE="\${INFILE%.vcf}.vep.vcf"
    VEP_CACHE="/sw/data/uppnex/vep/$params.vep_cache_version"
    ASSEMBLY="$params.assembly"

    case "\$INFILE" in
        *vcf) FORMAT="vcf" ;;
        *bed) FORMAT="ensembl" ;;
        *)    printf "Unrecognized format for '%s'\n" "\$INFILE" >&2
              exit 1;;
    esac

    # VEP failes files without variants, but the pipeline should still run
    if ! grep -qv '^#' "\$INFILE"; then
        cp "\$INFILE" "\$OUTFILE"
        exit 0
    fi

    if [ -z "\$SLURM_CPUS_ON_NODE" ]; then
        CPUS=1
    else
        CPUS=\$SLURM_CPUS_ON_NODE
    fi

    variant_effect_predictor.pl \
        -i "\$INFILE"              \
        --format "\$FORMAT"        \
        -cache --dir "\$VEP_CACHE" \
        -o "\$OUTFILE"             \
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
        \$( test "\$CPUS" -gt 1 && echo "--fork \$CPUS" ) \
        --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \
        --assembly "\$ASSEMBLY" \
        --offline
    """
}

process snpEff {
    input:
        set file(infile), val(uuid), val(dir) from ch_annotate_snpeff
    output:
        file '*.snpeff.vcf'

    tag "$uuid - $infile"
    publishDir "$dir", mode: 'copy', saveAs: { "$params.prefix$it" }

    executor choose_executor()

    when: 'snpeff' in workflowSteps

    script:
    """
    INFILE="$infile" ## Use bash-semantics for variables
    OUTFILE="\${INFILE%.vcf}.snpeff.vcf"

    snpEff -formatEff -classic ${params.assembly}.75 < "\$INFILE" > "\$OUTFILE"
    """
}


// Utility functions

def usage_message() {
    log.info 'WGS-structvar pipeline'
    log.info 'USAGE:'
    log.info 'Run a local copy of the wgs-structvar WF:'
    log.info '    nextflow main.nf --bam <bamfile> [more options]'
    log.info 'OR run from github:'
    log.info '    nextflow nbisweden/wgs-structvar --bam <bamfile> [more options]'
    log.info ''
    log.info 'The log file .nextflow.log will be produced when running and can be monitored'
    log.info 'by e.g. tail -f .nextflow.log'
    log.info 'More information about this pipeline can be found on:'
    log.info 'https://github.com/NBISweden/wgs-structvar'
    log.info ''
    log.info 'Options:'
    log.info '  Required'
    log.info '    --bam           Input bamfile'
    log.info '       OR'
    log.info '    --runfile       Input runfile for multiple bamfiles in the same run.'
    log.info '                    Whitespace separated, first column is bam file,'
    log.info '                    second column is output directory and an optional third column'
    log.info '                    with a run id to more easily keep track of the run (otherwise'
    log.info '                    it\'s autogenerated).'
    log.info '    --project       Uppmax project to log cluster time to'
    log.info '       OR'
    log.info '    -profile local  Run locally, only for very small datasets'
    log.info '  Optional'
    log.info '    --help          Show this message and exit'
    log.info '    --fastq         Input fastqfile (default is bam but with fq as fileending)'
    log.info '                    Used by fermikit, will be created from the bam file if'
    log.info '                    missing.'
    log.info '    --steps         Specify what steps to run, comma separated: (default: manta, vep)'
    log.info '                Callers: manta, fermikit'
    log.info '                Annotation: vep, snpeff'
    log.info '                Extra: normalize (with vt),'
    log.info '                       mask_cohort (with bed files in mask_cohort/)'
    log.info '    --sg_mask_ovlp  Fractional overlap for use with the filter option'
    log.info '    --no_sg_reciprocal  Don\'t use a reciprocal overlap for the filter option'
    log.info '    --outdir        Directory where resultfiles are stored (default: results)'
    log.info '    --prefix        Prefix for result filenames (default: no prefix)'
    log.info '    --mask_artifacts_dir'
    log.info '                    Directory with bed files for artifact filtering (default: mask_artifacts/)'
    log.info '    --mask_cohort_dir'
    log.info '                    Directory with bed files for cohort filtering (default: mask_cohort/)'
    log.info ''
}

def check_input_params() {
    error = false
    if (!params.project && workflow.profile != 'local') {
        log.info('You need to specify what project to run under')
        error = true
    }
    if (!params.bam && !params.runfile) {
        log.info('You need to specify a bam or runfile file')
        error = true
    } else if (params.bam && params.runfile) {
        log.info('You can only specify one of bam and runfile')
        error = true
    }
    if (error) {
        log.info('See --help for more information')
        exit 1
    }
}

def startup_message() {
    revision = grab_git_revision() ?: 'v0.2'

    log.info "======================"
    log.info "WGS-structvar pipeline"
    log.info "======================"
    log.info "Command line : $workflow.commandLine"
    log.info "Bamfile      : $params.bam"
    log.info "Scriptdir    : $baseDir"
    log.info "Revision     : $revision"
    log.info "Work dir     : $workDir"
    log.info "Output dir   : $params.outdir"
    log.info "Project      : $params.project"
    log.info "Will run     : " + workflowSteps.join(", ")
    log.info ""
}

def setup_input_channel_from_runfile(rf) {
    ch = Channel.from(
        file(rf).readLines().collect {
            res = it.split()
            if ( res.size() == 3 ) {
                (f, outdir, uuid) = res
            } else {
                (f, outdir) = res
                uuid = randomUUID() as String
            }
            [file(f), uuid, outdir]
        }
    )
    return ch
}

def setup_input_channel_from_bam(bf) {
    bamfile = file(bf)
    if (! bamfile.exists()) {
        exit 1, "The bamfile, '$bf', does not exist"
    }
    ch = Channel.from(bamfile).map {
        uuid = randomUUID() as String
        [file(it), uuid, params.outdir]
    }
    return ch
}

def grab_git_revision() {
	return workflow.revision ?: workflow.scriptId.substring(0,10)
}

def infer_bam_index_from_bam(f) {
    // If the ".bam.bai" file does not exist, try ".bai" without ".bam"
    return infer_filepath(f, /$/, '.bai')
        ?: infer_filepath(f, /.bam$/, '.bai')
        ?: filepath_from(f, /$/, '.bai') // Default filename if none exist
}

def infer_fastq_from_bam(f) {
    return infer_filepath(f, /.bam$/, '.fq.gz')
        ?: filepath_from(f, /.bam$/, '.fq.gz') // Default filename if none exist
}

def filepath_from(from, match, replace) {
    path = file( from.toString().replaceAll(match, replace) )
    return path
}

def infer_filepath(from, match, replace) {
    path = file( from.toString().replaceAll(match, replace) )
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
 * for a lot of our work, this overrides the slurm executor specified in
 * the -profile command line option */
def choose_executor() {
    if (workflow.profile == 'local') {
        return 'local'
    }
    return nextflow_running_as_slurmjob() ? 'local' : 'slurm'
}


def encode_directory(dirname) {
    parts = dirname.split("/")
    return parts.join("....")
}

def generate_filename(orig_string, prefix='') {
    parts = orig_string.split('....')
    if ( parts > 1 ) {
        return parts.join("/")
    }
    m = orig_string =~ /(.+)\.\.\.\.(.+)/
    if (m) {
        (full, f, d) = m[0]
        return "$d/$prefix$f"
    }
    return orig_string
}

def processWorkflowSteps(steps) {
    if ( ! steps ) {
        return []
    }

    workflowSteps = steps.split(',').collect { it.trim().toLowerCase() }

    if ('manta' in workflowSteps && 'fermikit' in workflowSteps) {
        workflowSteps.push( 'make_intersect' )
    }

    return workflowSteps
}

/* vim: set filetype=groovy */
