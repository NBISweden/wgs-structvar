# Whole Genome Sequenceing Structural Variation Pipelines

Pipelines for WGS Structural variation analysis.

Status: Work in progress.

## How to run reference piplines in Make

This command will run the pipelines implemented by @pallolason in Make (and
which are supposed to be re-implemented in Nextflow as part of this project),
on a small example dataset:

```bash
make -nf src/SV.makefile sample.SV.DEL.filt.bed VPATH=data
```

## External links

* [NextFlow website](http://www.nextflow.io)
* [NextFlow gitter chat](https://gitter.im/nextflow-io/nextflow)
