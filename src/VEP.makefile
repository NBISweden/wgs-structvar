###
# some usual VEP calls
#

# cmd: variant_effect_predictor.pl -i P1426_101.sv.sort.filt.vcf -cache -dir /pica/sw/apps/bioinfo/bcbio/20150205/milou/Cellar/vep/75_2014-06-12/lib/homo_sapiens/79_GRCh37/ -o stdout -regulatory --force_overwrite -gmaf

#VEP_HOME:=/sw/apps/bioinfo/bcbio/20150205/milou/Cellar/vep/75_2014-06-12
VEP_HOME:=/home/pallol/glob/apps/bioinfo/variant_effect_predictor
VEP:=$(VEP_HOME)/ensembl-tools-release-79/scripts/variant_effect_predictor/variant_effect_predictor.pl
#VEP_CACHE:=$(VEP_HOME)/lib/homo_sapiens/79_GRCh37
VEP_CACHE:=$(VEP_HOME)/.vep
# specifying vep cache with -dir_cache does not work for the bcbio installation


%.vcf.vep : %.vcf
	$(VEP) -i $< -cache --dir $(VEP_CACHE) -o $@.tmp --vcf --merged --regulatory --force_overwrite --sift b --polyphen b --symbol --numbers --biotype --total_length --canonical --ccds --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE --offline && mv $@.tmp $@

%.vcf.vep : %.vcf.gz
	zcat $< | $(VEP) -i $< -cache --dir $(VEP_CACHE) -o $@.tmp --vcf --merged --regulatory --force_overwrite --sift b --polyphen b --symbol --numbers --biotype --total_length --canonical --ccds --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE --offline && mv $@.tmp $@

%.vcf.vep : %.vcf.bgz
	zcat $< | $(VEP) -i $< -cache --dir $(VEP_CACHE) -o $@.tmp --vcf --merged --regulatory --force_overwrite --sift b --polyphen b --symbol --numbers --biotype --total_length --canonical --ccds --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE --offline && mv $@.tmp $@

%.vcf.vep.cadd : %.vcf.bgz
	zcat $< | $(VEP) -i $< -cache --dir $(VEP_CACHE) -o $@.tmp --plugin CADD,/sw/data/uppnex/ToolBox/whole_genome_SNVs_inclAnno.tsv.gz --vcf --merged --regulatory --force_overwrite --sift b --polyphen b --symbol --numbers --biotype --total_length --canonical --ccds --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE --offline && mv $@.tmp $@
