# Whole Genome Sequencing Structural Variation Pipeline

## Quick start

### Install nextflow

```bash
curl -fsSL get.nextflow.io | bash
mv ./nextflow ~/bin
```

### Run the pipeline

```bash
nextflow run NBISweden/wgs-structvar --project <uppmax_project_id> --bam <bamfile.bam> --steps manta,fermikit,vep
```

This will run both manta and fermikit, annotate the results with variant
effect predictor and create summary files for everything in the `results`
subdirectory.

To run everything locally (only for testing on very small datasets), pass `-profile local` and omit the
`--project` option.

It is recommended that you set the environment variable `NXF_WORK` to
something like

```bash
export NXF_WORK=$SNIC_NOBACKUP/work
```

Preferably in your `.bashrc`.

## General information

This is a pipeline for running the two structural variation callers fermikit
and manta on UPPMAX.

You can choose to run either of the two structural variation callers or both
(and generate summary files).

### Masking

The pipeline will use the following mask files to remove known artifacts:

* From [cc2qe/speedseq](https://github.com/cc2qe/speedseq): https://github.com/cc2qe/speedseq/raw/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed
* From [lh3/varcmp](https://github.com/lh3/varcmp): https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs37d5.bed.gz


## Detailed usage

### Command line options

```
Usage:
    nextflow main.nf --bam <bamfile> [more options]

Options:
  Required
    --bam           Input bamfile
    --project       Uppmax project to log cluster time to
  Optional:
    (default values in parenthesis where applicable)
    --help          Show this message and exit
    --fastq         Input fastqfile (default is bam but with fq as fileending)
                    Used by fermikit, will be created from the bam file if
                    missing.
    --steps         Specify what steps to run, comma separated (manta,vep):
                Callers: manta, fermikit
                Annotation: vep, snpeff
    --outdir        Directory where resultfiles are stored (results)
    --prefix        Prefix for result filenames ()
```


### Customization

The file `nextflow.config` can be used to make some further customizations to
the workflow.

It's probably only the `params` scope of the config file that is of interest
to customize.

The first part has the default values for the command line parameters, see the
usage message for information on them.

The next section has the reference assembly to use, both as fasta and assembly
name.

The `modules` section contains all modules used by the workflow and their
versions, change modules here not in the `main.nf` file.

Finally the `runtime` section has the different runtimes for the different
parts of the workflow. `fermikit` has it's own timespec since that is a very
long running program, otherwise the workflow differentiates between `callers`
and other supporting `simple` single-core jobs.


## External links

* [NextFlow website](http://www.nextflow.io)
* [NextFlow gitter chat](https://gitter.im/nextflow-io/nextflow)

[![Stories in Ready](https://badge.waffle.io/NBISweden/wgs-structvar.png?label=ready&title=Ready)](https://waffle.io/NBISweden/wgs-structvar)
