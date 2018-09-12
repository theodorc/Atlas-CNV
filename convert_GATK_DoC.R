
# --------------------------------------------------------------------------------
# R script that computes RPKM based on the GATK_DoC *.interval_summary file.
# and adds gene_exon, gender, and midpool
# Called by main atlas_cnv.pl v0, Sep 11, 2018. Ted Chiang
# Copyright 2016-2018, Baylor College of Medicine Human Genome Sequencing Center.
# All rights reserved.
# --------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	  stop("Missing GATK_DoC input file: *.interval_summary or wrong file format", call.=FALSE)
} 
gatk_interval_summary_file = args[1]
GATK_DoC_IS = read.table(file=gatk_interval_summary_file, header=T)
ggg = GATK_DoC_IS$average_coverage
names(ggg) = GATK_DoC_IS$Target
rrr=t(as.matrix(ggg))

panel_design = args[2]
PANEL = read.table(file=panel_design, header=T)

gender = args[3]
midpool = args[4]


# Reduce rrr matrix to only targets that are Call_CNV_Yes in panel_design.
rrr = rrr[1, which(PANEL$Call_CNV=='Y')] 
rrr = t(rrr)  # changes vector with names to matrix (w/o really transposing), same as: rrr = t(as.matrix(rrr))

# RPKM
# ----
# Get the total exon lengths  SUM of all exons{ (stop-start+1) }

exon_coors = strsplit(colnames(rrr), split=':')
exon_bp_coors = sapply(exon_coors, function(x){x[2]})
# start-stop+1 ...ie. exon length in bp.
exon_lengths = sapply(exon_bp_coors, function(x){ abs(eval(parse(text=x[])))+1 })
num_of_reads_in_exon = t(t(rrr)*exon_lengths) / 100
num_of_reads_in_sample_per_million = rowSums(num_of_reads_in_exon)/10^6
#write.table(num_of_reads_in_sample_per_million, file="Reads_in_sample_per_million.txt", sep="\t", quote=F)
ave_cov_10 = rrr * 10 
rpkm = ave_cov_10 / num_of_reads_in_sample_per_million
rpkm = t(rpkm)

#check if exon targets from GATK_Doc is identical to those in the panel design.
check = identical( colnames(rrr), as.vector(PANEL$Exon_Target[which(PANEL$Call_CNV=='Y')] ) )
if (check) { 
	rpkm = cbind(as.vector(PANEL$Gene_Exon[which(PANEL$Call_CNV=='Y')]), rownames(rpkm), rpkm, rep(gender, dim(rpkm)[1]), rep(midpool, dim(rpkm)[1])) 
}  else { 
	stop("convert_GATK_DoC.R: GATK_DoC Exon_Target coords do not match the panel design. Exit R and die.", call.=FALSE)
	}
rownames(rpkm) = NULL
colnames(rpkm) = c('Gene_Exon', 'Exon_Target', 'RPKM', 'Gender', 'Midpool')
write.table(rpkm, file=paste (gatk_interval_summary_file, "rpkm.txt", sep='.'), sep="\t", quote=F, row.names=F)

