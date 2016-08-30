# setup

REF:=/glob/pallol/apps/bioinfo/variant_effect_predictor/.vep/homo_sapiens/79_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
SNPEFFJAR:=/sw/apps/bioinfo/snpEff/4.1/milou/snpEff.jar

# decompose, normalize and annotate VCF with snpEff.
# NOTE: can also swap snpEff with VEP
#NOTE: -classic and -formatEff flags needed with snpEff >= v4.1
%.sort.vcf.gz : %.vcf.gz
	vcf-sort -c $< | bgzip > $@.tmp && mv $@.tmp $@ && tabix $@

#snpEff does not like INS where stop=start-1
%.sort.fix.vcf.gz : %.sort.vcf.gz
	 zcat $< |  perl -ane 'if (/^\#/){print; next} if (/END\=(\d+)\;/){if ($$1 < $$F[1]){$$F[7] =~ s/$$1/$$F[1]/; print join("\t", @F), "\n"} else {print ;}} else {print;}' | bgzip -c > $@.tmp && mv $@.tmp $@

# in case there are sorted/unsorted, go for the sorted one first
%.snpEff.vcf.bgz %.snpEff.vcf.bgz.tbi : %.sort.vcf.gz
	zcat $< | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r $(REF) - | java -Xmx7G -jar $(SNPEFFJAR) -formatEff -classic GRCh37.75 | bgzip -c > $@.tmp && mv $@.tmp $@
	tabix $@

%.snpEff.vcf.bgz : %.vcf.gz
	zcat $< | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r $(REF) - | java -Xmx7G -jar $(SNPEFFJAR) -formatEff -classic GRCh37.75 | bgzip -c > $@.tmp && mv $@.tmp $@
	tabix $@
# in case there are sorted/unsorted, go for the sorted one first
%.snpEff.vcf.bgz %.snpEff.vcf.bgz.tbi : %.sort.vcf.bgz
	zcat $< | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r $(REF) - | java -Xmx7G -jar $(SNPEFFJAR) -formatEff -classic GRCh37.75 | bgzip -c > $@.tmp && mv $@.tmp $@
	tabix $@

%.snpEff.vcf.bgz : %.vcf.bgz
	zcat $< | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r $(REF) - | java -Xmx7G -jar $(SNPEFFJAR) -formatEff -classic GRCh37.75 | bgzip -c > $@.tmp && mv $@.tmp $@
	tabix $@

#w reg
%.snpEff.reg.vcf.bgz %.snpEff.vcf.bgz.tbi : %.sort.vcf.gz
	zcat $< | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r $(REF) - | java -Xmx7G -jar $(SNPEFFJAR) -formatEff -classic GRCh37.75 -reg HeLa-S3 -reg NHEK | bgzip -c > $@.tmp && mv $@.tmp $@
	tabix $@

%.snpEff.reg.vcf.bgz : %.vcf.gz
	zcat $< | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r $(REF) - | java -Xmx7G -jar $(SNPEFFJAR) -formatEff -classic GRCh37.75 -reg HeLa-S3 -reg NHEK | bgzip -c > $@.tmp && mv $@.tmp $@
	tabix $@








#-reg NHEK