


#install.packages('VGAM')
library('VGAM')
library(dplyr)
source('s_asb_meta_data.R')

library(DESeq2)
source('s_asb_meta_data.R')
source('test/t_deseq_norm.r')





###### ASB Agreement between DEseq2 and binomial ############################

#source("http://bioconductor.org/biocLite.R")
#biocLite("DEseq2")


f_compare_DEseq_and_binomial <- function(loc_cell, loc_tf, base_dir, db_version, dp_threshold, fdr, deseq_ver = '2' ){
  library(DESeq)
  data_replicates = f_get_epi_from_paird_database(cell = loc_cell, loc_cell, tf=loc_tf,  f_p("%s/%s/%s/",base_dir, loc_cell, db_version),dp_threshold, fdr, file_suffix = 'database' ,raw_data = TRUE, 
                                                  cell_filter = loc_cell, het_filter=c('het'), labs= lab_list ,target_lab =NA, rep = '.', strong_allele_exchange = F )
  
  
  
  data_merge = f_get_epi_from_paird_database(cell = loc_cell, loc_cell, tf=loc_tf,  f_p("%s/%s/%s/",base_dir, loc_cell, db_version),dp_threshold, fdr, file_suffix = 'database' ,raw_data = FALSE, 
                                             cell_filter = loc_cell, het_filter=c('het'), labs=lab_list ,target_lab =NA, rep = '.', strong_allele_exchange = F )
  
  data_replicates = data_replicates[rownames(data_merge),]
  
  
  
  data_merge$total_dp = data_merge$ref_tf_dp + data_merge$alt_tf_dp
  
  #hist(data_merge$total_dp, breaks = 200)
  
  
  best_list = f_select_highest_target_lab_and_rep(data_replicates, loc_tf, lab_list)
  best_lab = best_list$best[4]
  replicate_count = sum(best_list$stats[best_lab,c('rep1','rep2','rep3')] > 0)
  ASB_data_replicates = data_replicates[rownames(subset(data_merge, ASB == 'ASB')),]
  
  ref_group = sort(grep(f_p('ref_%s_dp_%s[0-9]$', loc_tf, best_lab), colnames(ASB_data_replicates), value = TRUE))
  alt_group = sort(grep(f_p('alt_%s_dp_%s[0-9]$', loc_tf, best_lab), colnames(ASB_data_replicates), value = TRUE))
  
  ref_rep = ref_group[1]
  i = 1
  counts = c(0,0,0)
  for(i in 1:length(ref_group) ){
    rep_cols = grep(f_p('.*_%s%s',best_lab, i), c(ref_group, alt_group), value = TRUE)
    if(length(rep_cols) == 0) next
    counts[i] = sum(data_replicates[,rep_cols])
  }
  
  
  
  ref1_sign_ratio = (ASB_data_replicates[ref_group[1]] > ASB_data_replicates[alt_group[1]])
  ref2_sign_ratio = (ASB_data_replicates[ref_group[2]] > ASB_data_replicates[alt_group[2]])
  ref2_ratio = (ASB_data_replicates[ref_group[2]] + 0.1)/(ASB_data_replicates[ref_group[2]] + ASB_data_replicates[alt_group[2]] + 0.2)
  ref1_ratio = (ASB_data_replicates[ref_group[1]] + 0.1)/(ASB_data_replicates[ref_group[1]] + ASB_data_replicates[alt_group[1]] + 0.2)
  
  rep_correlation = cor(ref1_ratio, ref2_ratio)
  rep_sign_correlation = cor(ref1_sign_ratio, ref2_sign_ratio)
  #qplot(ref1_ratio$ref_mef2a_dp_haib1, ref2_ratio$ref_mef2a_dp_haib2) + geom_point(position = 'jitter')
  
  
  pasillaCountTable = data_replicates[rownames(data_merge), c(ref_group, alt_group)]
  pasillaCountTable= pasillaCountTable[rowSums(pasillaCountTable) >= 0,]
  
  
  pasillaDesign = data.frame(
    row.names = colnames( pasillaCountTable ),
    condition = rep(c('ref','alt'), each = length(alt_group)),
    libType = rep('single-end', each = 2*length(alt_group)) )
  
  pairedSamples = pasillaDesign$libType == "single-end"
  countTable = pasillaCountTable[ , pairedSamples ]
  condition = factor(pasillaDesign$condition[ pairedSamples ])
  
  if (deseq_ver == '1'){
    cds = newCountDataSet( countTable, condition )
    
    cds = estimateSizeFactors( cds )
    sizeFactors( cds )
    result =try(cds <- estimateDispersions( cds))
    if (class(result) == "try-error"){
      cds = estimateDispersions( cds, fitType = 'local')
    }
    normalizedCounts <- t( t(counts(cds)) / sizeFactors(cds) )
    plotDispEsts( cds )
    res = nbinomTest( cds, "alt", "ref")
    head(res)
    plotMA(res)
  }else if (deseq_ver == '2'){
    
    dds <- DESeqDataSetFromMatrix(countData = countTable,
                                  colData = pasillaDesign,
                                  design = ~ condition)
    
    dds <- DESeq(dds,)
    normalizedCounts <- t( t(counts(dds)) / sizeFactors(dds) )
    
    res <- as.data.frame(results(dds,independentFiltering = F))
    res$pval = res$pvalue
    res$id = rownames(res)
    #View(res)
    
  }
  
  if(TRUE){
    
    replicate_counts = data_replicates[rownames(data_merge), c(ref_group, alt_group)]
    
    return_list <- f_deseq_norm(replicate_counts, rep_cols)
    normalizedCounts = return_list$data
    normalized_ASB=f_asb_calling(as.data.frame(round(normalizedCounts)), loc_tf = loc_tf, target_lab = best_lab, lab_cols = colnames(normalizedCounts),calling_method = 'binomial', fdr, dp_threshold = 0)
    scale_factor = max(return_list$scale_factor)
    
  }
  

  edgeR_flag = F
  if (edgeR_flag == FALSE){
    library(edgeR)
    #x <- read.delim("fileofcounts.txt",row.names="Symbol")
    x = countTable

    group <- factor(rep(1:2, each = length(ref_group)))
    y <- DGEList(counts=x,group=group)
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    #et <- exactTest(y)
    
    
    group = y$samples$group
    design = model.matrix(~group)
    rownames(design) = rownames(y$samples)
    
    et <- glmFit(y, design)   
    et2 <- glmLRT(et)
    edgeR_predictions=topTags(et2,adjust.method = 'BH',n = nrow(x),)$table
    head(edgeR_predictions)
    #str(edgeR_predictions)
  }

  
  edgeR_pred=rownames(edgeR_predictions)[(edgeR_predictions$FDR < fdr)]
  
  rownames(res) = res$id

  tiff(f_p("%s/asb_pvalue/%s_%s_pvlaue_dist.tiff", writing_directory, loc_cell, loc_tf))
  par(mfrow=c(3,1)) 
  hist(res$pval, breaks=20, col="skyblue", border="slateblue", main="DEseq")
  hist(as.numeric(edgeR_predictions$PValue), breaks=20, col="skyblue", border="slateblue", main="EdgeR")
  hist(data_merge$p_value, breaks=20, col="skyblue", border="slateblue", main="Binomial")
  dev.off()
  
  
  null_pvalues = rownames(subset(res, pval > 0.99))
  head(res)
  non_na_ids =intersect(intersect(subset(res, !is.na(pval))$id, rownames(data_merge)), rownames(edgeR_predictions))
  asb_ids = intersect(non_na_ids, rownames(subset(data_merge, ASB == 'ASB')))
  res$sig = res$padj < 0.05
  
  ggplot(res, aes(baseMean, log(pval))) + geom_point(size = 0.1)
  
  ggplot(data_merge, aes(total_dp, fill = ASB)) + geom_histogram() + xlim(c(8,100))
  
  
  res_subset = res[ complete.cases(res) & res$baseMean > 5,  ]
  res_subset$p_adj = p.adjust(res_subset$pval, method ='fdr')
  
  
  perf_data = data.frame(binomial = data_merge[non_na_ids, 'p_value'], binomial_norm = normalized_ASB[non_na_ids, 'p_value'], deseq = res[non_na_ids, 'pval'], edgeR = edgeR_predictions[non_na_ids, 'PValue'] )
  rownames(perf_data) = non_na_ids
  perf_data = f_sort_by_col(perf_data, index = 'deseq')
  ASB_data = perf_data[asb_ids,]
  ggplot(perf_data, aes(log10(binomial), log10(deseq)) ) + geom_point()
  #ggplot(perf_data, aes((binomial), (deseq)) ) + geom_point()
   
  correlation = cor(perf_data$binomial, perf_data$deseq, method = 'spearman')
  edgeR_correlation = cor(perf_data$binomial, perf_data$edgeR, method = 'spearman')
  edge_overlap = intersect(edgeR_pred, asb_ids)
  bioNorm_overlap = intersect(rownames(subset(normalized_ASB, ASB=='ASB')), asb_ids)
  bioNorm_cor = cor(perf_data$binomial, perf_data$binomial_norm, method = 'spearman')
  
  asb_deseq_correlation = cor(ASB_data$binomial, ASB_data$deseq, method = 'spearman')
  asb_edgeR_correlation = cor(ASB_data$binomial, ASB_data$edgeR, method = 'spearman')
  asb_bioNorm_correlation = cor(ASB_data$binomial, ASB_data$binomial_norm, method = 'spearman')
  
  cat('\n Rank cordrelation of two methods:', correlation, '\n')
  cat('Deseq:' ,sum(res$sig, na.rm = T), 'EdgeR', sum(edgeR_predictions$FDR <= fdr), 'Binomial', sum(data_merge$ASB == 'ASB'), 'BinomialNorm', sum(normalized_ASB$ASB == 'ASB'),'\n')
  return_list = list(stats = c(correlation, edgeR_correlation, bioNorm_cor, sum(res$sig, na.rm = T), sum(edgeR_predictions$FDR <= fdr) ,sum(data_merge$ASB == 'ASB'), sum(normalized_ASB$ASB == 'ASB'),
                               rep_correlation, rep_sign_correlation, sort(counts), asb_deseq_correlation, asb_edgeR_correlation, length(edge_overlap),length(bioNorm_overlap), scale_factor, replicate_count),
                     pvalue = data.frame(DEseq2 = res$pval, Binomial = data_merge$p_value, edgeR = edgeR_predictions[rownames(data_merge),'PValue'])
                     )
  
  return (return_list)
}




source('test/t_beta_binomial_call.r')
f_compare_betaBinomial_and_binomial <- function(loc_cell, loc_tf, base_dir, db_version, dp_threshold, fdr){
  
  
  data_merge = f_get_epi_from_paird_database(cell = loc_cell, loc_cell, tf=loc_tf,  f_p("%s/%s/%s/",base_dir, loc_cell, db_version),dp_threshold, fdr, file_suffix = 'database' ,raw_data = FALSE, 
                                             cell_filter = loc_cell, het_filter=c('het'), labs=lab_list ,target_lab =NA, rep = '.', strong_allele_exchange = F , wgs_rm_flag =F)
  
  
  cmp_list=f_beta_binomial_call_wgs(data_merge, dp_threshold)
  cat('Dispersion:', cmp_list$dispersion)
  dy = cmp_list$dy
  
  guest_data = data_merge
  minor_change  = (abs(guest_data$pwm_ref_score - guest_data$pwm_alt_score) < 0.01)
  subset_data = subset(guest_data, best_match == 'Best' & minor_change & ASB == 'nonASB')
  
  
  pwm_list=f_beta_binomial_call_wgs(subset_data, dp_threshold)
  
  
  sum(data_merge$ASB == 'ASB')
  overall_cor=cor(data_merge$p_value, dy,method = 'spearman')
  
  data_merge$betabinomial = dy 
  
  fold_change = ((data_merge$ref_tf_dp + 0.01)/(data_merge$alt_tf_dp + 0.01) <= 0.66 | (data_merge$ref_tf_dp + 0.01)/(data_merge$alt_tf_dp + 0.01) >= 1.5)
  
  qplot(log(betabinomial), log(p_value), data=data_merge) + xlim(c(-120, 0)) + ylim(c(-120,0))
  
  ASB_data = subset(data_merge, ASB == 'ASB')
  ASB_cor=cor(ASB_data$p_value, ASB_data$betabinomial ,method = 'spearman') 
  
  
  
  Binomial_call=which(data_merge$ASB == 'ASB')
  Beta_call=which(p.adjust(data_merge$betabinomial, method = 'BH' ) < fdr )
  Binomil_call2 = which(p.adjust(data_merge$p_value, method = 'BH' ) < fdr & fold_change )
  cat('Binomial:', length(Binomil_call2), 'Beta:', length(Beta_call), 'Binomial2:', length(Binomial_call), '\n')
  intersect_call=intersect(Binomial_call, Beta_call)
  return (c(nrow(ASB_data), length(Beta_call), length(intersect_call), overall_cor, ASB_cor, as.numeric(cmp_list$dispersion), as.numeric(pwm_list$dispersion), nrow(subset_data) ))
}


#mtrace(f_compare_DEseq_and_binomial)

f_add_gm_cells_in_datalist<-function(cell_data_list){
  gm_cell_list = c('gm12864','gm12872','gm12873','gm19238','gm19239','gm19240')#,'gm12891','gm12892')
  
  for (loc_cell in gm_cell_list){
    
    cell_data_list[[loc_cell]][['ctcf']] = 0
    
  }
  
  return (cell_data_list)
  
  
}



f_get_batch_stats_DEseq_binomial <- function(cell_data_list, base_dir, db_version, dp_threshold, fdr){
  
  loc_cell  = 'gm12878'
  cell_loc = 'gm12878'
  tf_list = names(cell_data_list[[loc_cell]])
  
  loc_tf = 'ctcf'
  
  library('DESeq2')
  cor_stats = data.frame()
  pvalues= data.frame()
  edgeR_flag = 'TRUE'
  for(loc_cell in names(cell_data_list)){
    loc_tf_list = names(cell_data_list[[loc_cell]])
    for(loc_tf in loc_tf_list){
      dp_threshold = 10
      cat('TF', loc_tf,'\n')
      #debug(f_compare_DEseq_and_binomial)
      #mtrace(f_compare_DEseq_and_binomial)
      #undebug(f_compare_DEseq_and_binomial)
      loc_correlation=f_compare_DEseq_and_binomial(loc_cell, loc_tf, base_dir, db_version, dp_threshold, fdr)
      cor_stats = rbind(cor_stats,c(loc_cell, loc_tf, loc_correlation$stats))
      pvalues = rbind(pvalues, loc_correlation$pvalue )
      #stop()
    }
    
  }
  
  names(cor_stats) = c('cell', 'tf', 'correlation', 'edgeR_correlation' ,'bioNorm_cor', 'deseq', 'edgeR', 'binomial','bioNorm',
                       'replicate_cor', 'rep_sign_cor', 'rep3', 'rep2', 'rep1', 'deseq_ASB_cor','edgeR_ASB_cor', 'edgeR_overlap'
                       ,'bioNormOverlap', 'scale_factor', 'replicate_count')
  mean(as.numeric(cor_stats$correlation))
  min(as.numeric(cor_stats$correlation))
  mean(as.numeric(cor_stats$replicate_cor))
  mean(as.numeric(cor_stats$rep_sign_cor))
  number_change = sum(as.numeric(cor_stats$binomial))/sum(as.numeric(cor_stats$deseq))
  
  
  sum(as.numeric(cor_stats$deseq) > 40)
  sum(as.numeric(cor_stats$binomial) >= 40)
  
  cat('Number of events','Binomial: ', sum(as.numeric(cor_stats$binomial)), 'edgeR', sum(as.numeric(cor_stats$edgeR)) ,
      'DEseq', sum(as.numeric(cor_stats$deseq)), 
      'edgeR overlap', sum(as.numeric(cor_stats$edgeR_overlap)) ,'\n')
  cat('Mean edgeR correlation:', mean(as.numeric(cor_stats$edgeR_correlation)), 'overlap ratio with Binomial:', sum(as.numeric(cor_stats$edgeR_overlap))/sum(as.numeric(cor_stats$edgeR)) )
  dim(pvalues)
  hist(pvalues$DEseq2)
  hist(pvalues$Binomial)
  hist(pvalues$edgeR)
  
  
  
  ###Replicate lib size distribution########
  cor_stats$rep_sign_cor = as.numeric(cor_stats$rep_sign_cor)
  cor_stats$correlation = as.numeric(cor_stats$correlation)
  cor_stats$size_ratio = as.numeric(cor_stats$rep1)/(as.numeric(cor_stats$rep1) + as.numeric(cor_stats$rep2) + as.numeric(cor_stats$rep3))
  write.table(cor_stats, file = f_p('%s/cor_stats.txt', writing_directory),quote = F, sep = '\t', row.names = F, col.names = T)

  
  cat('Range of replicate correlation in ASB events mean',mean(as.numeric(cor_stats$rep_sign_cor)), 'min', min(as.numeric(cor_stats$rep_sign_cor)))
  cat('DEseq and Binomial correlation:', mean(as.numeric(cor_stats$correlation)), 'Min', min(as.numeric(cor_stats$correlation)))
  
  hist(cor_stats$correlation)
  
  #hist(cor_stats)
  
  return (cor_stats)
  
}
f_get_batch_beta_and_binomial <- function(cell_data_list, base_dir, db_version, dp_threshold, fdr){
  
  loc_cell  = 'gm12872'
  cell_loc = 'gm12878'
  tf_list = names(cell_data_list[[loc_cell]])
  
  loc_tf = 'ctcf'
  
  library('DESeq2')
  cor_stats = data.frame()
  pvalues= data.frame()
  for(loc_cell in names(cell_data_list)){
    loc_tf_list = names(cell_data_list[[loc_cell]])
    for(loc_tf in loc_tf_list){
      dp_threshold = 10
      cat('\n', loc_cell, 'TF', loc_tf,'\n')
      #debug(f_compare_DEseq_and_binomial)
      #mtrace(f_compare_DEseq_and_binomial)
      #undebug(f_compare_DEseq_and_binomial)
      #mtrace(f_compare_betaBinomial_and_binomial)
      #mtrace(f_beta_binomial_call)
      loc_correlation=f_compare_betaBinomial_and_binomial(loc_cell, loc_tf, base_dir, db_version, dp_threshold, fdr)
      cor_stats = rbind(cor_stats,c(loc_cell, loc_tf, loc_correlation))
    }
    
  }
  
  names(cor_stats) = c('cell', 'tf', 'binomial' , 'Beta', 'overlap', 'overall_cor','ASB_cor', 'dispersion')
 
  
  cat('Number of events','Binomial: ', sum(as.numeric(cor_stats$binomial)), 'Beta', sum(as.numeric(cor_stats$Beta)),'Overlap', sum(as.numeric(cor_stats$overlap)),'\n')
  return (cor_stats)
  
}



f_main<- function(){
  lab_list = c('sydh','uw','uta', 'haib')
  load(file = './data/tmp/cell_data_list')
  
  cell_data_list = f_add_gm_cells_in_datalist(cell_data_list)
  names(cell_data_list)
  
  #cell_data_list = cell_data_list[c('gm12872')]
  
  #debug(f_compare_DEseq_and_binomial)
  #debug(f_compare_betaBinomial_and_binomial)
  #debug(f_beta_binomial_call)
  
  cell_data_list[c('gm12878','gm12872')]
  
  edgeR_stats=f_get_batch_stats_DEseq_binomial(cell_data_list, base_dir, db_version, dp_threshold, fdr)
  
  mean(as.numeric(edgeR_stats$bioNorm_cor))
  mean(as.numeric(edgeR_stats$bioNormOverlap)/as.numeric(edgeR_stats$binomial))
  
  #beta_stats = f_get_batch_beta_and_binomial(cell_data_list, base_dir, db_version, dp_threshold, fdr)
  
  #mean(as.numeric(edgeR_stats$edgeR_correlation))
  #mean(edgeR_stats$correlation)
  
  #write.table(beta_stats, file = f_p('%s/beta_stats.txt', writing_directory),quote = F, sep = '\t', row.names = F, col.names = T)
  
  return (edgeR_stats)
  
  
}

beta_stats = f_main()
edgeR_stats = beta_stats
head(edgeR_stats)
edgeR_stats$scale_factor = as.numeric(edgeR_stats$scale_factor)
edgeR_stats$bioNorm_overlapRatio = as.numeric(edgeR_stats$bioNormOverlap)/as.numeric(edgeR_stats$bioNorm)
mean(edgeR_stats$bioNorm_overlapRatio)
edgeR_stats$bioNorm = as.numeric(edgeR_stats$bioNorm)

sum(as.numeric(edgeR_stats$bioNormOverlap))


cor_test=cor(edgeR_stats$bioNorm_overlapRatio, edgeR_stats$scale_factor,method = 'spearman')

scale_plot<-ggplot(edgeR_stats, aes(x = bioNorm_overlapRatio, y = scale_factor )) + geom_point() + 
  xlab('Percentage of called ASB events in normalized approach \noverlapped with direct-sum approach') +
  ylab('Scale factor of larger replicate') + theme_Publication(base_size = 9)

#normalized_grid<-arrangeGrob(scale_plot, density_plot, widths = c(4,3), default.units = "in", ncol=2)
ggsave(filename = f_p('%s/scale_factor_and_overlap.tiff', writing_directory ),plot = scale_plot, width = 5, height = 5, units = 'in')


range(edgeR_stats$bioNorm_cor)
head(edgeR_stats)
edgeR_stats %>% filter(scale_factor > 2, bioNorm_overlapRatio < 0.9)


beta_stats[,3:8] = data.matrix(beta_stats[,3:8])
print(beta_stats %>% summarise(mean_cor = mean(overall_cor), mean_asb_cor = mean(ASB_cor)))








