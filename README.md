# Whole Genome Sequencing Structural Variation Pipelines

Pipelines for WGS Structural variation analysis.

Status: Work in progress.

## How to run prototype pipeline in Make

This command will dry-run the prototype pipeline implemented by @pallolason in Make (and
which are supposed to be re-implemented in Nextflow as part of this project),
on a small example dataset:

```bash
make -nf src/SV.makefile sample.SV.DEL.filt.bed VPATH=data
```

## External links

* [NextFlow website](http://www.nextflow.io)
* [NextFlow gitter chat](https://gitter.im/nextflow-io/nextflow)

[![Stories in Ready](https://badge.waffle.io/NBISweden/wgs-structvar.png?label=ready&title=Ready)](https://waffle.io/NBISweden/wgs-structvar)
