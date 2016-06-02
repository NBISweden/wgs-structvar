###
# Top level makefile to consolidate "best practise" calls for SVs
#
# © pall.olason@scilifelab.se
#

###
# General stuff defs:
#

# TODO - fix install perl script somewhere central or get from repo
REF_FASTA:=/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37.fasta
SVVCF2BED:=/proj/b2014152/nobackup/private/pallol/src/SVvcf2bed.pl
INTERSECTBED:=/sw/apps/bioinfo/BEDTools/2.23.0/milou/bin/intersectBed

###
# General rules
#

%.SV.DEL.bed : %.fermikit.DEL.bed %.manta.DEL.bed
	$(INTERSECTBED) -a $*.fermikit.DEL.bed -b $*.manta.DEL.bed -f 0.5 -r | sort -k1,1V -k2,2n > $@.tmp && mv $@.tmp $@	

%.SV.DUP.bed : %.fermikit.DUP.bed %.manta.DUP.bed
	$(INTERSECTBED) -a $*.fermikit.DUP.bed -b $*.manta.DUP.bed -f 0.5 -r | sort -k1,1V -k2,2n > $@.tmp && mv $@.tmp $@	

%.SV.INS.bed : %.fermikit.INS.bed %.manta.INS.bed
	$(INTERSECTBED) -a $*.fermikit.INS.bed -b $*.manta.INS.bed -f 0.5 -r | sort -k1,1V -k2,2n > $@.tmp && mv $@.tmp $@	





###
# formatting targets
#

%.sort.vcf : %.vcf.gz
	vcf-sort -c $< > $@.tmp && mv $@.tmp $@

%.sort.vcf : %.vcf
	vcf-sort -c $< > $@.tmp && mv $@.tmp $@
##
# comfort targets ... very unspecific....
#

%.DEL.bed : %.sv.bed
	cat $< | grep -w DEL > $@.tmp && mv $@.tmp $@

%.INS.bed : %.sv.bed
	cat $< | grep -w INS > $@.tmp && mv $@.tmp $@

%.DUP.bed : %.sv.bed
	cat $< | grep -w DUP > $@.tmp && mv $@.tmp $@





%.bed : %.vcf
	$(SVVCF2BED) $<  > $@.tmp && mv $@.tmp $@

%.bed : %.vcf.gz
	zcat $< | $(SVVCF2BED) -  > $@.tmp && mv $@.tmp $@

%.vcf.gz : %.vcf
	bgzip $<

%.vcf.gz.tbi : %.vcf.gz
	tabix -p vcf $<


###
# Manta specific stuff:
#



NUMTHREADS:=16
MANTA_HOME:=/sw/apps/bioinfo/manta/0.27.1/milou/bin
CONFIG_MANTA:=$(MANTA_HOME)/configManta.py


.PRECIOUS: %.mantaDir

%.mantaDir : %.bam
	$(CONFIG_MANTA) --normalBam $< --referenceFasta $(REF_FASTA) --runDir $@

%.mantaRun %.mantaDir/results/variants/diploidSV.vcf.gz : %.mantaDir
	cd $< && ./runWorkflow.py -m local -j $(NUMTHREADS)

%.manta.sv.vcf.gz : %.mantaDir/results/variants/diploidSV.vcf.gz
	mv $< $@.tmp && mv $@.tmp $@


###
# Fermikit specific stuff:
#


GENOME_SIZE=3g
GENOME_REF_FILE=/sw/data/uppnex/reference/Homo_sapiens/g1k_v37/concat/human_g1k_v37.fasta
THREADS=16
READLEN=150
VCFSORT:=/sw/apps/bioinfo/vcftools/0.1.14/milou/bin/vcf-sort
FERMIKIT_HOME:=/sw/apps/bioinfo/fermikit/r178/milou/fermi.kit
RUN_CALLING:=$(FERMIKIT_HOME)/run-calling
FERMI2.PL:=$(FERMIKIT_HOME)/fermi2.pl
SAMTOOLS:=/sw/apps/bioinfo/samtools/1.3/milou/bin/samtools
BGZIP:=/sw/apps/bioinfo/tabix/0.2.6/milou/bin/bgzip


###
# in case we only have bam and no fq:
#
%.fq.gz : %.clean.dedup.recal.bam
	$(SAMTOOLS) bam2fq $< | gzip  - > $@.tmp && mv $@.tmp $@
%.fq.gz : %.bam
	$(SAMTOOLS) bam2fq $< | gzip  - > $@.tmp && mv $@.tmp $@


%.mak : %.fq.gz
	echo "set read length with READLEN=X - using default read length $(READLEN)"
	$(FERMI2.PL) unitig -s$(GENOME_SIZE) -t$(THREADS) -l$(READLEN) -p $* $< > $@.tmp && mv $@.tmp $@ 



%.run %.sv.vcf.gz : %.mak
	make -f $<
	$(RUN_CALLING) -t$(THREADS) $(GENOME_REF_FILE) $*.mag.gz | sh


%.fermikit.sv.vcf.gz : %.sv.vcf.gz
	$(VCFSORT) -c $< > $*.fermikit.sv.vcf
	$(BGZIP) -f $*.fermikit.sv.vcf





#.DUMMY : MODULES ECHOMOD
MODULES : 
	module  use /proj/b2013006/sw/modules
	module load fermikit

ECHOMOD :
	@echo -e "you need to load module fermi before running:\nmodule use /proj/b2013006/sw/modules\nmodule load fermikit"


###
# makefile to filter bed/vcf files against dbs of known/promiscuous regions
# and/or make track of dubious regions within the dataset under inspection
#




BEDTOOLS:=/glob/pallol/apps/bin/bedtools
LCRREGIONS:=mask_regions/LCR-hs37d5.bed.gz
ABNCOVREGIONS:=mask_regions/ceph18.b37.lumpy.exclude.2014-01-15.bed
1KGSVS:=mask_regions/ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.esv.vcf


$(ABNCOVREGIONS) :
	mkdir -p mask_regions
	wget -P mask_regions https://github.com/cc2qe/speedseq/raw/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed
$(LCRREGIONS) : 
	mkdir -p mask_regions
	wget -P mask_regions https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs37d5.bed.gz

%.vcf : %.vcf.gz
	gunzip -c $< > $@

%.sort.vcf : %.vcf
	vcf-sort -c $< > $@.tmp && mv $@.tmp $@

%.flcr.vcf : %.vcf $(LCRREGIONS)
	cat $< | $(BEDTOOLS) intersect -header -v -a stdin -b $(LCRREGIONS) > $@.tmp && mv $@.tmp $@
%.fcov.vcf : %.vcf $(ABNCOVREGIONS)
	cat $< | $(BEDTOOLS) intersect -header -v -a stdin -b $(ABNCOVREGIONS) > $@.tmp && mv $@.tmp $@
%.f1kg.vcf : %.vcf
	cat $< | $(BEDTOOLS) intersect -header -v -a stdin -b $(1KGSVS) > $@.tmp && mv $@.tmp $@


###
# BED support
#
# filter out those where 25% of "a" overlap LCR regions
#

# remove events with 25% overlap with LCR regions:
%.lcrfilt.bed : %.bed $(LCRREGIONS)
	$(BEDTOOLS) intersect -v -a $< -b $(LCRREGIONS) -f 0.25 > $@.tmp && mv $@.tmp $@

# remove events with 25% overlap with LUMPY abnormal coverage regions:
%.abncovfilt.bed : %.bed $(ABNCOVREGIONS)
	$(BEDTOOLS) intersect -v -a $< -b $(ABNCOVREGIONS) -f 0.25 > $@.tmp && mv $@.tmp $@

# and combined
%.filt.bed : %.lcrfilt.abncovfilt.bed
	mv $< $@

# remove "known" events from 1000g??? alternatively common events from a reference population
# eg swedish references or Umea300

%.popfilt.bed : %.bed $(POPBED)
	$(BEDTOOLS) intersect -v a $< -b $(POPBED) -f 0.25 > $@.tmp && mv $@.tmp $@

###
# DEPTH filter, this is experimental stuff and probably obsolete
#

INPUTLIST:=
DEPTHTH:=5

FERMI3LIST:=fermi.min3.promisc.bed # need to restructure this stuff
FERMI7LIST:=fermi.min7.promisc.bed # need to restructure this stuff


%.fermi3.vcf : %.vcf $(FERMI3LIST)
	$(BEDTOOLS) intersect -header -v -a $< -b $(FERMI3LIST) > $@.tmp && mv $@.tmp $@

%.fermi7.vcf : %.vcf $(FERMI7LIST)
	$(BEDTOOLS) intersect -header -v -a $< -b $(FERMI7LIST) > $@.tmp && mv $@.tmp $@

%.$(DEPTHTH).vcf: %.vcf $(INPUTLIST)
	$(BEDTOOLS) multiinter -i $(INPUTLIST) | perl -ane 'BEGIN{$$th = shift} print if ($$F[3] > $$th) ' | $(BEDTOOLS) merge > $@.tmp && mv $@.tmp $@

%.filt.$(DEPTHTH).vcf : %.vcf %.$(DEPTHTH).vcf
	$(BEDTOOLS) intersect -header -v -a $< -b %.$(DEPTHTH).vcf > $@.tmp && mv $@.tmp $@

# bedtools multiinter -i *.bed | perl -ane 'print if ($F[3] > 2) ' | bedtools merge | bedtools intersect -v -a stdin -b ../mask_regions/ceph18.b37.lumpy.exclude.2014-01-15.bed  | bedtools intersect -v -a stdin -b ../mask_regions/LCR-hs37d5.bed.gz > bd.min3.promisc.bed
# bedtools multiinter -i *.sv.sort.vcf.gz| perl -ane 'print if ($F[3] > 2) ' | bedtools merge |  bedtools intersect -v -a stdin -b ../mask_regions/ceph18.b37.lumpy.exclude.2014-01-15.bed  | bedtools intersect -v -a stdin -b ../mask_regions/LCR-hs37d5.bed.gz > fermi.min3.promisc.bed
