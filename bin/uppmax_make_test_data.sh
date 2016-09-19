#!/bin/bash

module load bioinfo-tools
module load samtools

set -x

CHR="1 2"
FN="1-2"

BASEDIR=/proj/a2009002/webexport/opendata/HiSeqX_CEPH/CEP-1-7/03-BAM
BASEFILE=CEP-1-7.clean.dedup.recal

OUT="$BASEFILE.chr$FN"

if [ ! -f $OUT.bam ]; then
    samtools view -bh "$BASEDIR/$BASEFILE.bam" $CHR > $OUT.bam
fi

if [ ! -f $OUT.bam.bai ]; then
    samtools index $OUT.bam
fi

if [ ! -f $OUT.fq.gz ]; then
    samtools bam2fq $OUT.bam | gzip - > $OUT.fq.gz
fi
