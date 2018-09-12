# --------------------------------------------------------------------------------
# atlas_cnv.R v0. called by main atlas_cnv.pl v0, Sep 11, 2018. Ted Chiang
# Copyright 2016-2018, Baylor College of Medicine Human Genome Sequencing Center.
# All rights reserved.
# --------------------------------------------------------------------------------
library(reshape2)
library("optparse")
library(plotrix)
 
option_list = list(
  make_option(c("--rpkm_file"), type="character", default=NULL, 
              help="RPKM matrix data filename", metavar="character"),
  make_option(c("--panel_file"), type="character", default="HG19_CAREseqv2_eMERGE3_LD_panel_design_3565_exons.tab", 
              help="panel design filename. 3 tab-delimited columns: 'Exon_Target', 'Gene_Exon', 'Call_CNV'
                [default = %default]", metavar="character"),
  make_option(c("--threshold_del"), type="numeric", default="-0.6", 
              help="log2 ratio cutoff to call a CNV deletion. [default= %default]", metavar="numeric"),
  make_option(c("--threshold_dup"), type="numeric", default="0.4", 
              help="log2 ratio cutoff to call a CNV duplication. [default= %default]", metavar="numeric"),
  # ???????#make_option(c("--threshold_exonQC_Zscore"), type="numeric", default="0.2", 
  #            help="ExonQC cutoff to fail an exon target, ie. stddev > this value are flagged as failed. [default= %default]", metavar="numeric"),
  make_option(c("--threshold_sampleQC"), type="numeric", default="0.2", 
              help="SampleQC cutoff to fail a sample, ie. stddev > this value are flagged as failed. [default= %default]", metavar="numeric"),
  make_option(c("--threshold_sample_anova"), type="numeric", default="0.05", 
              help="In the coefficients of one way ANOVA of sample means, the Pr(>|t|) cutoff to flag a sample. [default= %default]", metavar="numeric")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$rpkm_file)) {
   print_help(opt_parser)
   stop("RPKM matrix datafile must be supplied (--rpkm_file).", call.=FALSE)
}
if (!file.exists(opt$panel_file)) {
   print_help(opt_parser)
   stop("Panel design file is not supplied, or default cannot be found. (--panel_file).", call.=FALSE)
}

# ccc: df of rpkms for all samples, all targets, ie. RPKM_matrix.IRC_MB-MID0075.CLEMRG_1pA (already 2066 columns, ie. only Call_CNV=Y)
rpkm_matrix_file = opt$rpkm_file
ccc = read.table(file=rpkm_matrix_file, header=T, check.names=F, row.names=1)

# Midpool summary result outputfile & midpool directory
midpool_summary_results = gsub('RPKM_matrix', 'atlas_cnv_summary', rpkm_matrix_file) 
midpool_dir = gsub('/.*', '/', rpkm_matrix_file)

# ppp: df of the panel design file with 3 columns: ppp$Exon_Target, ppp$Gene_Exon, ppp$Call_CNV, ppp$RefSeq
panel_design_file = opt$panel_file
ppp = read.table(file=panel_design_file, header=T, check.names=F)

# =============================================================================
# CNV default cutoffs:
# 2^(0.5)  = 1.414214....  so if RPKM ratio is >= 1.41 then, it's a DUP call.
# 2^(0.4692898) = 1.3844278

# 1-2^(-0.7) = 0.3844278

# 2^(-0.7) = 0.615572....  so if RPKM ratio is <= 0.61 then, it's a DEL call.  <--OLD
# 2^(-0.415) = 0.75....  so if RPKM ratio is <= 0.75 then, it's a DEL call.
# =============================================================================
threshold_del          = opt$threshold_del
threshold_dup          = opt$threshold_dup
#threshold_exonQC       = opt$threshold_exonQC
threshold_sampleQC     = opt$threshold_sampleQC
threshold_sample_anova = opt$threshold_sample_anova


# some in-house functions:
find_median_sample = function (x) {
	median_rpkm = mmm[x]
	targ_coor = names(mmm[x])   # "1:1220087-1220186"
	# ccc$"1:1220087-1220186"  works but not as a variable ccc$"targ_coor".  
	# ccc[targ_coor]     <-- column slice
	# ccc[[targ_coor]]   <-- column vector, same as ccc$"1:1220087-1220186" 
	median_rpkm = max(ccc[[targ_coor]][ which(ccc[[targ_coor]] <= median_rpkm) ]);  
	sample_idx = which(ccc[[targ_coor]] == median_rpkm)
	sample_idx = sample_idx[1]           # in case there are duplicate samples of 0 rpkm, then just pick one.
	sample_id  = rownames(ccc)[sample_idx]
	med = c(targ_coor, median_rpkm, sample_id, sample_idx)
	return (med)
}
my_row_sd = function(x) {
	aaa = DDD[x,][FFF[x,]]
	mysd = sd( aaa )
	return (mysd) 
}
my_row_count = function(x) {
	count = sum(FFF[x,])
	# cat ('n(for stddev cal)=', count, "\n")
	return (count) 
}
my_col_sd = function(x) {
	aaa = DDD[,x][FFF[,x]]
	mysd = sd( aaa )
	return (mysd) 
}
my_col_count = function(x) {
	count = sum(FFF[,x])
	# cat ('n(for stddev cal)=', count, "\n")
	return (count) 
}
my_col_sd_wo_outliers = function(x) {
	aaa = DDD_wo_outliers[,x][FFF_wo_outliers[,x]]
	mysd = sd( aaa )
	return (mysd) 
}
my_col_x_based_on_zscore_neg3_wo_outliers = function(x) {
	aaa = DDD_wo_outliers[,x][FFF_wo_outliers[,x]]
	mymean = mean( aaa )
	mysd = sd( aaa )
	#myx = (-4)*mysd + mymean
	myx = (-2.576)*mysd + mymean
	return (myx) 
}
my_col_x_based_on_zscore_pos3_wo_outliers = function(x) {
	aaa = DDD_wo_outliers[,x][FFF_wo_outliers[,x]]
	mymean = mean( aaa )
	mysd = sd( aaa )
	#myx = (4)*mysd + mymean
	myx = (2.576)*mysd + mymean
	return (myx) 
}
pass_fail = function(x) { 
	if (x==TRUE) { return ("Fail") }
	else         { return ("Pass") }
}
call_cnvs = function (x) {   # calling cnv per sample via idx of ccc df
	sample_id = rownames(fff)[x]
	# --------------------------------------------------------------------------------------------------------------------------
	# check to see if sample failed sample QC or ANOVA.  If so, name the file:  .cnv.FAILED_sampleQC or .cnv.FAILED_sampleANOVA
	# --------------------------------------------------------------------------------------------------------------------------
	anova_sample_id = paste('sample_id', sample_id, sep='')
	if       (is.element(sample_id, failed_samples)  &  !is.element(anova_sample_id, rownames(failed_samples_by_anova_pval_lt_5pct) ))  {  
		outfile = paste(midpool_dir, sample_id, '.cnv.FAILED_sampleQC', sep='')
	}
	else if  (!is.element(sample_id, failed_samples)  &  is.element(anova_sample_id, rownames(failed_samples_by_anova_pval_lt_5pct) )) {
		outfile = paste(midpool_dir, sample_id, '.cnv.FAILED_sampleANOVA', sep='')
	}
	else if  (is.element(sample_id, failed_samples)  &  is.element(anova_sample_id, rownames(failed_samples_by_anova_pval_lt_5pct) )) {
		outfile = paste(midpool_dir, sample_id, '.cnv.FAILED_sampleQC_and_sampleANOVA', sep='')
	}
	else {
		outfile = paste(midpool_dir, sample_id, '.cnv', sep='')
	}
	header = rbind(c('Exon_Target', 'Gene_Exon', 'cnv', 'log2R', 'rpkm', 'median_rpkm', 'Exon_Status', 'E_StDev', 'c_Score', 'RefSeq'))
	write.table(header, file=outfile, append=T, quote=F, row.names=F, col.names=F, sep="\t")
	# call dels
	cat ("Calling dels on:", sample_id)
	sample_dels=NULL
	#dels = which(fff[x,1:(ncol(fff)-2)]<=threshold_del)
	dels = which(fff[x,1:(ncol(fff)-2)]<=threshold_del & fff[x,1:(ncol(fff)-2)]<=threshold_del_soft)
	if (length(dels)==0) { cat (" has no cnv dels.\n") }
	else {
		cat (" has cnv dels.\n")
		Gene_Exon = ppp$Gene_Exon[ppp$Exon_Target %in% colnames(fff)[dels]]
		RefSeq = ppp$RefSeq[ppp$Exon_Target %in% colnames(fff)[dels]]
		cnv = rep('del', length(dels))
		log2R = as.data.frame( t(fff[x, dels]) )
		rpkm  = as.data.frame( t( ccc[x, colnames(ccc) %in% colnames(fff)[dels]] ) )
		median_rpkm = mmm$median_rpkm[mmm$targ_coor %in% colnames(fff)[dels]]
		ExonQC_PF = is.element(colnames(fff)[dels], rownames(failed_exons_by_ExonQC))
		ExonQC_PF = sapply(ExonQC_PF, pass_fail)
		ExonQC_value = as.numeric(fff["ExonQC", dels])
		Confidence_Score = round(log2R/ExonQC_value, 2)
		sample_dels = data.frame( as.character(Gene_Exon), cnv, log2R, rpkm, median_rpkm, ExonQC_PF, ExonQC_value, Confidence_Score, RefSeq, stringsAsFactors=FALSE )
		colnames(sample_dels) = c('Gene_Exon', 'cnv', 'log2R', 'rpkm', 'median_rpkm', 'Exon_Status', 'E_StDev', 'c_Score', 'RefSeq')
		rownames(sample_dels) = colnames(fff)[dels]
		write.table(sample_dels, file=outfile,  append=T, quote=F, row.names=T, col.names=F, sep="\t")
	}
	# call dups
	cat ("Calling dups on:", sample_id)
	sample_dups=NULL
	#dups = which(fff[x,1:(ncol(fff)-2)]>=threshold_dup)
	dups = which(fff[x,1:(ncol(fff)-2)]>=threshold_dup & fff[x,1:(ncol(fff)-2)]>=threshold_dup_soft)
	if (length(dups)==0) { cat (" has no cnv dups.\n") }
	else {
		cat (" has cnv dups.\n")
		Gene_Exon = ppp$Gene_Exon[ppp$Exon_Target %in% colnames(fff)[dups]]
		RefSeq = ppp$RefSeq[ppp$Exon_Target %in% colnames(fff)[dups]]
		cnv = rep('dup', length(dups))
		log2R = as.data.frame( t(fff[x, dups]) )
		rpkm  = as.data.frame( t( ccc[x, colnames(ccc) %in% colnames(fff)[dups]] ) )
		median_rpkm = mmm$median_rpkm[mmm$targ_coor %in% colnames(fff)[dups]]
		ExonQC_PF = is.element(colnames(fff)[dups], rownames(failed_exons_by_ExonQC))
		ExonQC_PF = sapply(ExonQC_PF, pass_fail)
		ExonQC_value = as.numeric(fff["ExonQC", dups])
		Confidence_Score = round(log2R/ExonQC_value, 2)
		sample_dups = data.frame( as.character(Gene_Exon), cnv, log2R, rpkm, median_rpkm, ExonQC_PF, ExonQC_value, Confidence_Score, RefSeq, stringsAsFactors=FALSE )
		colnames(sample_dups) = c('Gene_Exon', 'cnv', 'log2R', 'rpkm', 'median_rpkm', 'Exon_Status', 'E_StDev', 'c_Score', 'RefSeq')
		rownames(sample_dups) = colnames(fff)[dups]
		write.table(sample_dups, file=outfile,  append=T, quote=F, row.names=T, col.names=F, sep="\t")
	}
	# plots
	plot_cnvs(rbind(sample_dels, sample_dups), sample_id)
}
plot_cnvs = function(x,y) {   # x is df of sample_dels/dups, y is sample_id
	pdffile = paste(midpool_dir, y, '.pdf', sep='')
	pdf.options(bg='honeydew')
        pdf(file=pdffile);
	par(mfrow=c(2,1))
	sample_id = gsub('[_A-Z]', '', y)  
	if (length(nrow(x))==0) {
		plot.new()
		text (0.5, 1, paste('No called CNVs for: ', sample_id, sep=''))
		garbage = dev.off()
	} 
	else {
		# gene plots
		cnv_genes = rle( gsub('-.*', '', x$Gene_Exon) )   # rle is like unix:uniq,  rm the transcript & exon numbering of Gene_Exon:  MTHFR-001_11, ie. '-001_11'.
		#cnv_genes_2plot = cnv_genes$values[ which(cnv_genes$lengths >= 1) ]  # Do gene plots w/ 1+ exons
		cnv_genes_2plot = cnv_genes$values[ which(cnv_genes$lengths >= 2) ]  # Do gene plots w/ 2+ exons
		if (length(cnv_genes_2plot) > 0) {   # are there any genes to plot?
			 for (g in 1:length(cnv_genes_2plot)) {
				gene = cnv_genes_2plot[g]
				gene_pattern = paste('^', cnv_genes_2plot[g], '-', sep='')  #  bug: grep 'NF2' gets 'RNF207' 
				cat ('Gene plot for sample: ', y, 'using gene pattern: ', gene_pattern, "\n")
				gene_exon_all = ppp$Gene_Exon[ grep(gene_pattern, ppp$Gene_Exon)  ]  # ppp: df of the panel design file with 3 columns: ppp$Exon_Target, ppp$Gene_Exon, ppp$Call_CNV
				gene_targ_coor_all = ppp$Exon_Target[ grep(gene_pattern, ppp$Gene_Exon)  ]  # ppp: df of the panel design file with 3 columns: ppp$Exon_Target, ppp$Gene_Exon, ppp$Call_CNV
				gene_targ_coor_start = which(colnames(fff) == gene_targ_coor_all[1])
				gene_targ_coor_stop  = which(colnames(fff) == gene_targ_coor_all[length(gene_targ_coor_all)])
				lll = fff[y, gene_targ_coor_start:gene_targ_coor_stop]   # df: one row w/column slice of log2R for given gene targets.
				lll = 2^(unlist(lll))  # convert log2R back to ratio for barplot
				names(lll) = gene_exon_all
				cnv_gene_exons = x$Gene_Exon[ grep(gene_pattern, x$Gene_Exon) ]
				col_vector=rep('grey', length(gene_exon_all))
				col_vector[ gene_exon_all %in% cnv_gene_exons ] = 'red'
				mean_Confidence_Score = round( mean( x$c_Score[ grep(gene_pattern, x$Gene_Exon) ] ), 2)
				barplot(lll, col=col_vector, ylim=c(0,2), ylab='sample / median' , las=2, cex.names=0.5, axis.lty=1, main='Gene plot', space=0)
				mtext(gene, cex=0.7)
				legend("top", c(paste('mean c_Score', mean_Confidence_Score, sep=' = ')), col=c('black'), pch=15, cex=0.5, bg="transparent")
				legend("topright", c('c_Score: (-)loss, (+)gain', 'High: abs(c_Score) > 5', 'Med: abs(c_Score) = 3-5', 'Low: abs(c_Score) < 3'), col=c('black', 'grey', 'grey', 'grey'), pch=15, cex=0.5, bg="transparent")
				#legend("topright", c('dip: S/M = 0.85-1.15', 'del(het): S/M = 0.4-0.6', 'del(hom): S/M = 0-0.1', 'dup: S/M = >1.2': ), col=c('black', 'black', 'black', 'black'), pch=15, cex=0.5, bg="transparent")
				axis(1, at=match(cnv_gene_exons,names(lll))-0.5, labels=cnv_gene_exons, las=2, cex.axis=0.5, col.axis='red')
				abline(h=0.5)
				abline(h=1)

				table_rows = grep(gene_pattern, x$Gene_Exon)   # each plot.new() can only handle 12 rows per table, so make appropriate number of plots.
				if (length(table_rows) <= 12) {
					plot.new()
					p_start = 1
					p_stop = length(table_rows)
					addtable2plot(x='top', table=x[ table_rows[p_start:p_stop], 1:8], cex=0.75, bty='o', hlines=T, vlines=T, bg='honeydew', xpad=0.25, ypad=1)
				}
				else if (length(table_rows) > 12) {
					num_of_tables = floor(length(table_rows)/12)
					p_stop = 0
					for ( p in 1:num_of_tables ) {
						plot.new()
						p_start = p_stop + 1
						p_stop  = p*12
						addtable2plot(x='top', table=x[ table_rows[p_start:p_stop], 1:8], cex=0.75, bty='o', hlines=T, vlines=T, bg='honeydew', xpad=0.25, ypad=1)
					}
					if ((length(table_rows) %% 12) > 0) {    # a last table if needed.
						plot.new()  # last table
						p_start = p_stop + 1
						p_stop = length(table_rows)
						addtable2plot(x='top', table=x[ table_rows[p_start:p_stop], 1:8], cex=0.75, bty='o', hlines=T, vlines=T, bg='honeydew', xpad=0.25, ypad=1)
					}
			 	}
			 }
		}
		# exon plots
		for (i in 1:nrow(x)) {
			cnv_targ_coor = rownames(x)[i]
			cnv_gene_exon = x$Gene_Exon[i]
			cnv_exon_qc_value  = fff["ExonQC",cnv_targ_coor]
			cnv_exon_qc_status = x$Exon_Status[i]
			cnv_exon_confidence_score = x$c_Score[i]
			lll = fff[1:(nrow(fff)-2), cnv_targ_coor , drop=FALSE]    # df: column slice of log2R for given target.  Note:drop=FALSE keeps dimnames when subsetting to one dim.
			nnn = gsub('[_A-Z]', '', rownames(lll))                   # remove the '_HMWCYBCXX' and 'IDMB', for x-axis labels in the plot.
			j = match(cnv_targ_coor, mmm$targ_coor)
			median_rpkm       = mmm$median_rpkm[j]
			median_sample_id  = mmm$sample_id[j]
			median_sample_id  = gsub('[_A-Z]', '', median_sample_id)
			median_sample_idx = mmm$sample_idx[j]
			col_vector=rep('grey', dim(ccc)[1])
			threshold_del_soft_this_exon = threshold_del_soft[[cnv_gene_exon]]
			threshold_dup_soft_this_exon = threshold_dup_soft[[cnv_gene_exon]]
			col_vector[ which(lll<=threshold_del) ] = 'darkgoldenrod3' # hom/het dels
			col_vector[ which(lll>=threshold_dup) ] = 'darkgoldenrod3' # dups
			col_vector[ which(lll<=threshold_del & lll<=threshold_del_soft_this_exon) ] = 'darkgoldenrod3' # hom/het dels
			col_vector[ which(lll>=threshold_dup & lll>=threshold_dup_soft_this_exon) ] = 'darkgoldenrod3' # dups
			col_vector[ which(rownames(ccc)==y)] = 'red'      # this given sample
			col_vector[ median_sample_idx ] = 'blue'          # median sample used for reference
			lll=as.vector(t(2^(lll)))          # convert log2R back to ratio for barplot
			names(lll)=nnn                     # as.vector conversion causes rownames to be lost
			barplot(lll, col=col_vector, ylim=c(0,2), ylab='sample / median', las=2, cex.names=0.5, axis.lty=1, main='Exon plot', space=0)
			mtext(paste(cnv_gene_exon, cnv_targ_coor, sep=', '), cex=0.7)
			legend("top", c(paste('exon status', cnv_exon_qc_status, sep=': '), paste('E[StDev]', cnv_exon_qc_value, sep=' = '), paste('c_Score', cnv_exon_confidence_score, sep=' = ')), col=c('black', 'black', 'black'), pch=15, cex=0.5, bg="transparent")
			legend("topright", c('sample', 'median'), col=c('red', 'blue'), pch=15, cex=0.5, bg="transparent")
			axis(1, at=match(sample_id,names(lll))-0.5, labels=sample_id, las=2, cex.axis=0.5, col.axis='red')
			axis(1, at=match(median_sample_id,names(lll))-0.5, labels=median_sample_id, las=2, cex.axis=0.5, col.axis='blue')
			abline(h=0.5)
			abline(h=1)
		}
		garbage = dev.off() 
	}
}
remove_5pc_outliers = function (aaa) {
 for (i in seq_along(aaa)) {
        a = which( aaa[[i]] < ( -1.96*sd(aaa[[i]]) + mean(aaa[[i]]) ) )
        b = which( aaa[[i]] > (  1.96*sd(aaa[[i]]) + mean(aaa[[i]]) ) )
        aaa[a,i]=NA
        aaa[b,i]=NA
 }
 return(aaa)
}





# mmm: df of medians for each target: rownames(targ_coor) ==> colnames(mmm$median_rpkm, mmm$sample_id, mmm$sample_idx)
# ie. the artificial sample of reference medians.
# mmm = apply(ccc, 2, median)  # 2 for column, 1 for row.
ccc_wo_outliers = remove_5pc_outliers(ccc)
mmm = apply(ccc_wo_outliers, 2, median, na.rm=TRUE)
mmm = lapply(seq_along(mmm), find_median_sample)
mmm = data.frame(matrix(unlist(mmm), ncol=4, byrow=T), stringsAsFactors=FALSE) # suppresses data.frame()’s default behaviour which turns strings into factors.
colnames(mmm) = c('targ_coor', 'median_rpkm', 'sample_id', 'sample_idx')
mmm = data.frame(targ_coor=mmm$targ_coor, median_rpkm=as.integer(mmm$median_rpkm), sample_id=mmm$sample_id, sample_idx=as.integer(mmm$sample_idx), stringsAsFactors=FALSE)

# ddd: df of log2 ratios of: samples (rows) vs. targets
ddd = as.data.frame(t( log2( t(ccc)/mmm$median_rpkm) ), stringsAsFactors=FALSE)   # a matrix has to be coerced into a dataframe!
DDD = as.matrix(ddd)    # apply function requires matrix.
FFF = is.finite(DDD)    # T/F matrix to rm the Inf, NaN.
# ---------------------------------------------------------------------------------------------
# Same as above 3, but on ccc_wo_outliers. THIS IS FOR COMPUTING THE ExonQC_wo_outliers vector.
# ---------------------------------------------------------------------------------------------
ddd_wo_outliers = as.data.frame(t( log2( t(ccc_wo_outliers)/mmm$median_rpkm) ), stringsAsFactors=FALSE)   # a matrix has to be coerced into a dataframe!
DDD_wo_outliers = as.matrix(ddd_wo_outliers)    # apply function requires matrix.
FFF_wo_outliers = is.finite(DDD_wo_outliers)    # T/F matrix to rm the Inf, NaN.
# ---------------------------------------------------------------------------------------------
# Vector of exon stddev wo_outliers 
ExonQC             = sapply(seq_along(colnames(ddd)), my_col_sd)
ExonQC_wo_outliers = sapply(seq_along(colnames(ddd_wo_outliers)), my_col_sd_wo_outliers)
threshold_del_soft = sapply(seq_along(colnames(ddd_wo_outliers)), my_col_x_based_on_zscore_neg3_wo_outliers)    # The 'x' log2Ratio threshold (soft and data-derived) to call deletions (wo/outliers)
threshold_dup_soft = sapply(seq_along(colnames(ddd_wo_outliers)), my_col_x_based_on_zscore_pos3_wo_outliers)    # The 'x' log2Ratio threshold (soft and data-derived) to call deletions (wo/outliers)
names(ExonQC)             = ppp$Gene_Exon[ which(ppp$Call_CNV=='Y') ]
names(ExonQC_wo_outliers) = ppp$Gene_Exon[ which(ppp$Call_CNV=='Y') ]
names(threshold_del_soft) = ppp$Gene_Exon[ which(ppp$Call_CNV=='Y') ]
names(threshold_dup_soft) = ppp$Gene_Exon[ which(ppp$Call_CNV=='Y') ]
exon_qc_outfile             = gsub('RPKM_matrix', 'ExonQC', rpkm_matrix_file) 
exon_qc_wo_outliers_outfile = gsub('RPKM_matrix', 'ExonQC_wo_outliers', rpkm_matrix_file) 
threshold_del_soft_outfile = gsub('RPKM_matrix', 'Exon_threshold_del_soft_wo_outliers', rpkm_matrix_file) 
threshold_dup_soft_outfile = gsub('RPKM_matrix', 'Exon_threshold_dup_soft_wo_outliers', rpkm_matrix_file) 
suppressWarnings(write.table(ExonQC,             file=exon_qc_outfile,             append=T, quote=F, row.names=T, col.names=T, sep="\t"))
suppressWarnings(write.table(ExonQC_wo_outliers, file=exon_qc_wo_outliers_outfile, append=T, quote=F, row.names=T, col.names=T, sep="\t"))
suppressWarnings(write.table(threshold_del_soft, file=threshold_del_soft_outfile, append=T, quote=F, row.names=T, col.names=T, sep="\t"))
suppressWarnings(write.table(threshold_dup_soft, file=threshold_dup_soft_outfile, append=T, quote=F, row.names=T, col.names=T, sep="\t"))
# 95%, 97.5%, 99%, 99.95%, 99.99% => zscore: 1.645, 1.96, 2.576, 3.291, 4
#threshold_exonQC = 2.576 * sd(ExonQC_wo_outliers) + mean(ExonQC_wo_outliers)
threshold_exonQC = 3.291 * sd(ExonQC_wo_outliers) + mean(ExonQC_wo_outliers)
#threshold_exonQC = 100  # To disable ExonQC metric.
names(threshold_exonQC) = 'threshold_exonQC'
cat ("ExonQC threshold based on 99% of sd distribution is: ", threshold_exonQC, "\n")
cat ("CNV Exon threshold del/dup: ", threshold_del, ", ", threshold_dup, "\n")
#sample_size = round(0.1*length(ExonQC_wo_outliers_for_Call_CNV_Y))  # sample 10% of the data to compute a exonQC threshold.
#threshold_exonQC = NULL
#for (i in 1:1) {
#	hhh = sample(ExonQC_wo_outliers_for_Call_CNV_Y, size=sample_size, replace=F)
#	threshold_exonQC = c(threshold_exonQC, 1.96 * sd(hhh) + mean(hhh))   # the sd cutoff is 95% of this background sampling
#}
#threshold_exonQC = mean(threshold_exonQC)
#cat ("ExonQC [ave of 1 events, each sampling 10% of bkgrd exon SD thresholded at 95% is]:", threshold_exonQC, "\n")


# add SampleQC stddev of log2 ratios, and row counts of data w/o Inf, NaN.
SampleQC               = sapply(seq_along(rownames(ddd)), my_row_sd)
Sample_count_wo_InfNaN = sapply(seq_along(rownames(ddd)), my_row_count)
# add ExonQC of target stddev of log2 ratios, and col counts of data w/o Inf, NaN.
Exon_count_wo_InfNaN   = sapply(seq_along(colnames(ddd)), my_col_count)

# ==> as two columns at the right end of matrix.
fff = cbind(  ddd, data.frame( SampleQC = SampleQC, Sample_count_wo_InfNaN = Sample_count_wo_InfNaN  ) )
# ==> as two rows at bottom end of matrix.
eee = as.data.frame(rbind( c(ExonQC_wo_outliers, NA, NA), c(Exon_count_wo_InfNaN, NA, NA) ))
rownames(eee) = c('ExonQC', 'Exon_count_wo_InfNaN')
colnames(eee) = colnames(fff)
fff = rbind(fff, eee)
fff = round(fff, 2)
rm(DDD, FFF, eee)

# fff: df of ddd w/ removed panel targets not used for CNV calling ie. (rm N ==> ppp$Call_CNV)
# fff = ddd[, which(ppp$Call_CNV=='Y') ]
# fff = cbind( fff, SampleQC=ddd$SampleQC, Sample_count_wo_InfNaN=ddd$Sample_count_wo_InfNaN )

# rrr: df of sample rpkm means & stddev: rownames(sampleids) ==> colnames(rrr$rpkm_mean, rrr$rpkm_stddev)
rrr = data.frame(rpkm_mean = rowMeans(ccc), rpkm_stddev = apply(ccc, 1, sd))
rrr = round(rrr, 2)
sss = data.frame(rpkm_mean = rowMeans(ccc), rpkm_stddev = apply(ccc, 1, sd), SampleQC = SampleQC)
sss = round(sss, 2)
# write sample rpkm means & stddev results to outfile:
sample_mean_stddev_outfile = gsub('RPKM_matrix', 'Sample_RPKM-means-stddevs_log2-stddevs', rpkm_matrix_file) 
suppressWarnings(write.table(sss, file=sample_mean_stddev_outfile, append=T, quote=F, row.names=T, col.names=T, sep="\t"))


# c_score matrix and pvalue matrix
cscore_outfile = gsub('RPKM_matrix', 'Cscore_matrix', rpkm_matrix_file)
pval_outfile   = gsub('RPKM_matrix', 'Pval_matrix', rpkm_matrix_file) 
c_scores = as.matrix(  t( t(ddd) / ExonQC_wo_outliers)  )   # apply function requires matrix.
#pvals    = as.data.frame( apply(c_scores, 2, function(x) { 2*(1 - pnorm(abs(x))) } ),  stringsAsFactors=FALSE )  # a matrix has to be coerced into a dataframe!        two tailed pvalues
pvals   = as.data.frame( apply(c_scores, 2, function(x) { 1 - pnorm(abs(x)) } ),  stringsAsFactors=FALSE )  # a matrix has to be coerced into a dataframe!          one tail pvalues
c_scores = as.data.frame(c_scores)
c_scores = round(c_scores, 2)
c_scores = cbind(Samples = rownames(c_scores), c_scores)
pvals    = round(pvals, 2)
pvals    = cbind(Samples = rownames(pvals), pvals)
suppressWarnings(write.table(c_scores, file=cscore_outfile, append=T, quote=F, row.names=F, col.names=T, sep="\t"))
suppressWarnings(write.table(pvals, file=pval_outfile, append=T, quote=F, row.names=F, col.names=T, sep="\t"))


# Failed exons by Exon QC > computed exon sd without outliers  (row slice is a dataframe already)
fff_idx = which(fff["ExonQC",]>threshold_exonQC)
if (length(fff_idx)==1) {    # when row slice is only 1 datapoint which is not automatically a row slice.
	fff_idx = which(fff["ExonQC",]>threshold_exonQC)
	failed_exons_by_ExonQC = as.data.frame(fff["ExonQC", fff_idx])
	rownames(failed_exons_by_ExonQC) = colnames(fff)[fff_idx]
	colnames(failed_exons_by_ExonQC) = 'ExonQC'
} else {   # when row slice has 0 datapoints or 2 or more which would be a row slice.
	failed_exons_by_ExonQC = as.data.frame(  t( fff["ExonQC", which(fff["ExonQC",]>threshold_exonQC)] )  )
}
Exon_count_wo_InfNaN   = as.data.frame(  t( fff["Exon_count_wo_InfNaN", rownames(failed_exons_by_ExonQC)] )  )
failed_exons_by_ExonQC = data.frame( 
				Gene_Exon=as.character(ppp$Gene_Exon[ppp$Exon_Target %in% rownames(failed_exons_by_ExonQC)]), 
				failed_exons_by_ExonQC, 
				Exon_count_wo_InfNaN, 
				stringsAsFactors=FALSE
				)
# Failed samples by SampleQC > 0.2
failed_samples = rownames(fff)[which(fff$SampleQC>threshold_sampleQC)]
failed_samples_by_SampleQC = data.frame( 
				failed_samples = failed_samples,
				sd_SampleQC = fff$SampleQC[which(fff$SampleQC>threshold_sampleQC)], 
				Sample_count_wo_InfNaN = fff$Sample_count_wo_InfNaN[which(fff$SampleQC>threshold_sampleQC)],
				rpkm_mean = rrr$rpkm_mean[rownames(rrr) %in% failed_samples],
				stringsAsFactors=FALSE
				)
# sample rpkm stats
sample_rpkm_stats = data.frame(sample_stats = c(min(rrr$rpkm_mean), max(rrr$rpkm_mean), median(rrr$rpkm_mean), min(rrr$rpkm_stddev), max(rrr$rpkm_stddev), median(rrr$rpkm_stddev)))
rownames(sample_rpkm_stats) = c('min_rpkm', 'max_rpkm', 'median_rpkm', 'min_stddev', 'max_stddev', 'median_stddev')

# Failed samples by anova, use rrr$rpkm_mean
aaa = melt(t(ccc))   # from reshape2 library to collapse df to data for lm.
colnames(aaa) = c('targ_coor', 'sample_id', 'rpkm')
fstatistic = anova(lm(rpkm ~ sample_id, data=aaa))$F[1]
pvalue     = anova(lm(rpkm ~ sample_id, data=aaa))$P[1]
s = summary(lm(rpkm ~ sample_id, data=aaa))
anova = data.frame(fstatistic,pvalue)
failed_samples_by_anova_pval_lt_5pct = data.frame(Prob_gt_t = s$coefficients[ which(s$coefficients[,4]<threshold_sample_anova), 4])

# write midpool summary results of above 5 dfs:
suppressWarnings(write.table(threshold_exonQC,                     file=midpool_summary_results, append=T, quote=F, row.names=T, col.names=T, sep="\t"))
suppressWarnings(write.table(failed_exons_by_ExonQC,               file=midpool_summary_results, append=T, quote=F, row.names=T, col.names=T, sep="\t"))
suppressWarnings(write.table(failed_samples_by_SampleQC,           file=midpool_summary_results, append=T, quote=F, row.names=F, col.names=T, sep="\t"))
suppressWarnings(write.table(sample_rpkm_stats,                    file=midpool_summary_results, append=T, quote=F, row.names=T, col.names=T, sep="\t"))
suppressWarnings(write.table(anova,                                file=midpool_summary_results, append=T, quote=F, row.names=T, col.names=T, sep="\t"))
suppressWarnings(write.table(failed_samples_by_anova_pval_lt_5pct, file=midpool_summary_results, append=T, quote=F, row.names=T, col.names=T, sep="\t"))


# Calling CNVs on fff dataframe: call dels first since it creates the .cnv outputfile.
garbage = sapply(seq_along(rownames(ccc)),  call_cnvs)



# -------under development, think about this-----------#
# dels_targs = cbind(colnames(fff)[dels], rep(1, length(dels)))
# dels_targs = split( dels_targs, seq(nrow(dels_targs)) )
# x = c("7:6013030-6013173", 1)
# targ_zscore = lapply( dels_targs, get_targ_zscore)
# get_targ_zscore = function (x) {  # x is a list.
# 	i = x[1]
# 	j = as.integer(x[2])
# 	x1 = ccc[[i]][j]
# 	xbar =  mean(ccc[[i]])
# 	mu   =  sd(ccc[[i]])
# 	zscore = (x1 - xbar)/mu
# 	return(zscore)
# }
# targ_pval = # write function to get the pval of the targ_zscore above.

