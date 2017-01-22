
setwd(dir = 'E:/Projects/R/')
source('s_function.R')
source('./test/t_get_epi_from_database3.r')
source('s_file_function.R', chdir = T)
source('s_library.R', chdir = T)
source('s_ggplot2_theme.R', chdir = T)
library(gridExtra)
source('s_asb_meta_data.R')
source('./ASB_calling/t_beta_binomial_call.r')
library('futile.logger')

library(logging)
library(dplyr)


fig_directory = './ASB_calling/figures'
dir.create(fig_directory)
options(stringsAsFactors = FALSE)

############Functions#################
#binomial_pwm_data = all_data
method = 'binom'

f_method_metrics <- function(method, binomial_pwm_data){
  #dnase_loc = cor(binomial_pwm_data$dnase_ratio > 0.5, binomial_pwm_data[[f_p('ASB_%s', method)]] == 'ASB')
  
  
  binomial_pwm_data = subset(binomial_pwm_data, total_dp < 20)
  
  test_obj= cor.test(binomial_pwm_data$dnase_ratio, -log(binomial_pwm_data[[f_p('pvalue_%s', method)]]), method = 'spearman',
                                                      use = 'pairwise.complete.obs')

  dnase_loc = test_obj$estimate

  
  
  #Concordance, remove NA.
  #dnase_loc = sum((binomial_pwm_data$dnase_ratio > 0.6) == (binomial_pwm_data[[f_p('ASB_%s', method)]] == 'ASB'), na.rm = T)/sum(!is.na(binomial_pwm_data$dnase_ratio))
  
  if(all(binomial_pwm_data$class == 'class3') | length(unique(binomial_pwm_data[[f_p('ASB_%s', method)]])) ==1 ){
    loc_estimate = 0
  }else{
    loc_estimate=fisher.test(binomial_pwm_data$class == 'class1', binomial_pwm_data[[f_p('ASB_%s', method)]] == 'ASB')$estimate
  }
  
  #    if( length(unique(binomial_pwm_data$dnase_ratio > 0.5)) ==1 | length(unique(binomial_pwm_data[[f_p('ASB_%s', method)]])) ==1 ){
  #      dnase_loc = 0
  #    }else{
  #      test_obj = fisher.test(binomial_pwm_data$dnase_ratio > 0.5, binomial_pwm_data[[f_p('ASB_%s', method)]] == 'ASB')
  #      if (test_obj$p.value < 0.01){
  #        dnase_loc = test_obj$estimate
  #        
  #      }else{
  #        
  #        dnase_loc = 0
  #      }
  #      
  #    }
  return_vector = c(dnase_loc, test_obj$p.value, loc_estimate, sum(binomial_pwm_data[[f_p('ASB_%s', method)]] == 'ASB'))
  names(return_vector) = c(f_p('dnase_%s', method), f_p('estimate_%s', method))
  
  return ( return_vector)  
                  
}

f_get_significant_correlation_coeff <- function(input_data){
  target_methods = c('binom', 'likelihood','edgeR','betaPool', 'betaRep')
  target_cols = paste0('pvalue_', target_methods)
  
  #install.packages('psych')
  library("psych")
  #data(sat.act)
  a=corr.test(input_data[,target_cols])
  a$r[a$p > 0.05] = 0
  
  return (a$r)
}

f_normalize_by_scale_factor <- function(input_data, target_cols, scale_factor){
  
  return_data = input_data
  repf_group = grep('ref', target_cols)
  
  for(i in 1:(length(target_cols)/2)){
    rep_i_cols = grep(f_p('.*_dp_.*%s$',i), target_cols, value = T)
    return_data[,rep_i_cols] = input_data[,rep_i_cols]/scale_factor[i]
  }
  return (data=return_data[,target_cols])
}



###############Sync##########
#system('./winscp_sync.bat')

########Part 1, feature analysis#########

#PWM score
#tf_list = c('ctcf', 'znf143', 'pu1', 'ebf1', 'bhlhe40')


raw_stats = read.table(file = f_p('./result/s_database_count/%s.txt', guest_cell), sep ='\t', header = TRUE)
tf_list_pwm = c( subset(raw_stats, ASB > ASB_threshold & PWM == TRUE)$tf, 'batf', 'runx3' )
tf_list   = subset(raw_stats, ASB > ASB_threshold )$tf
print(paste(tf_list, collapse = '\', \''))


tf_data_list = list()

statistic_table = data.frame(tf=tf_list, ASB_pwm_pvalue = 0, nonASB_pwm_pvalue = 0, ASB_dnase_pvalue = 0, nonASB_dnase_pvalue = 0,
                             snps_access_number =0, dp_filter = 0, ASB_num = 0, loop_end_pvalue = 1, loop_region_pvalue = 1,
                             ASB_best = 0, peak_pwm_ratio = 0)


rownames(statistic_table) = tf_list


########1 load the data#####
#' ##Loading the ASB data
#' 
#' 
cell_list = c('gm12878','helas3')


#cell_list = c('gm12878')

library(edgeR)
source('./ASB_calling/s_asbcall.R')
source('./ASB_calling/s_replicates_betabinomial.R')


cell_loc = 'gm12878'
tf = 'ctcf'
db_version = 'asb_call'

read_data = FALSE
#Pax5 is from pax5c20 data. Change all the pax5c20 to pax5 in the colnames.
if(read_data == TRUE){
  
  cell_tf_list = list()
  for (cell_loc in cell_list){
    
    tf_data_list = list()
    raw_stats = read.table(file = f_p('./result/s_database_count/%s.%s.txt', cell_loc, db_version), sep ='\t', header = TRUE)
    tf_list_pwm = subset(raw_stats, ASB > ASB_threshold & PWM == TRUE)$tf
    #tf_list   = grep('Total', subset(raw_stats, ASB > ASB_threshold)$tf, value = T, invert = T) 
    tf_list = raw_stats$tf
    
    for (tf in tf_list){
      
      f_header_in_loop(header = f_p('%s-%s', cell_loc, tf),level = 4)
      
      #undebug(f_get_epi_from_paird_database)
      guest_data_raw = f_get_epi_from_paird_database(cell_loc, cell_loc, tf, f_p("%s/%s/%s/",base_dir, cell_loc, db_version),dp_threshold, fdr,
                                                     cell_filter = cell_loc, het_filter=c("het"), raw_data=TRUE, diffpeak_method = diffpeak_method)
      
      tf_data_list[[tf]] = guest_data_raw
    }
    
    cell_tf_list[[cell_loc]] = tf_data_list
  }
  
  save(cell_tf_list,file = './ASB_calling/cell_tf_list')
  
}else{
  load(file = './ASB_calling/cell_tf_list')
}

str(cell_tf_list, max.level = 2)

source('./test/t_deseq_norm.r')


###2 Merge all the het sites for DHS####
scale_factor=list()
for (cell_loc in cell_list){
  raw_stats = read.table(file = f_p('./result/s_database_count/%s.%s.txt', cell_loc, db_version), sep ='\t', header = TRUE)
  tf_list = raw_stats$tf
  
  dnase_data = data.frame()
  
  for (tf in tf_list){
    
    guest_data_raw = cell_tf_list[[cell_loc]][[tf]]
    
    tf_dnase_data = guest_data_raw[,grep('dnase_dp', colnames(guest_data_raw), value = T)]
    
    tf_dnase_data$loc = rownames(tf_dnase_data)
    rownames(tf_dnase_data) = NULL
    dnase_data = rbind(dnase_data, tf_dnase_data)
  }
  
  dim(dnase_data)
  if (FALSE){
    scale_factor[[cell_loc]] = c(1)
  }else{
    dnase_colnames=sort(grep('dnase_dp', colnames(dnase_data), value = T))
    dnase_data_rmdup = dnase_data[!duplicated(dnase_data$loc) & rowSums(dnase_data[,dnase_colnames] == 0) != 4,]
    scale_factor[[cell_loc]]=f_deseq_norm(input_data = dnase_data_rmdup[,dnase_colnames], grep('ref', dnase_colnames, value = T))$scale_factor 
  }

}



#####3 Data Analysis#####

cor_matrix = data.frame(matrix(0, nrow= 5, ncol = 5))
ASB_collection = data.frame() 
stats_collect = data.frame()

dim(stats_collect)

dim(ASB_collection)
tf ='ebf1'
cell_loc = 'gm12878'

for (cell_loc in cell_list){
    
    raw_stats = read.table(file = f_p('./result/s_database_count/%s.%s.txt', cell_loc, db_version), sep ='\t', header = TRUE)
    tf_list_pwm = subset(raw_stats, ASB > ASB_threshold & PWM == TRUE)$tf
    #tf_list   = grep('Total', subset(raw_stats, ASB > ASB_threshold)$tf, value = T, invert = T) 
    tf_list = raw_stats$tf
    
    for (tf in tf_list){
      
      f_header_in_loop(header = f_p('%s-%s', cell_loc, tf),level = 4)
      
      #undebug(f_get_epi_from_paird_database)
      guest_data_raw = cell_tf_list[[cell_loc]][[tf]]
      statistic_table[tf,'snps_access_number']=nrow(subset(guest_data_raw, het_type == 'het'))
      
      #debug(f_get_epi_from_paird_database)
      #undebug(f_get_epi_from_paird_database)
      loc_rep = '.'
      
      dnase_cols = grep('(ref|alt)_dnase_dp', colnames(guest_data_raw), value = T)
      
      guest_data_raw[,dnase_cols]=(f_normalize_by_scale_factor(input_data = guest_data_raw, target_cols = dnase_cols, scale_factor =scale_factor[[cell_loc]] ))
      #guest_data_raw$ref_dnase_dp_broad2 = NULL
      #guest_data_raw$alt_dnase_dp_broad2 = NULL
      
      binomial_data = f_get_epi_from_paird_database(cell_loc, cell_loc, tf, f_p("%s/%s/%s/",base_dir, cell_loc, db_version), dp_threshold = 10, fdr,
                                                 cell_filter = cell_loc, het_filter=c("het"), add_het = TRUE,
                                                 db_data = guest_data_raw, rep=loc_rep, labs= lab_list, target_lab = NA, diffpeak_method = diffpeak_method,
                                                 calling_method = 'binomial')
      
      
      table(binomial_data$ASB)
      likelihood_data = f_get_epi_from_paird_database(cell_loc, cell_loc, tf, f_p("%s/%s/%s/",base_dir, cell_loc, db_version), dp_threshold = 10, fdr,
                                                    cell_filter = cell_loc, het_filter=c("het"), add_het = TRUE,
                                                    db_data = guest_data_raw, rep=loc_rep, labs= lab_list, target_lab = NA, diffpeak_method = diffpeak_method,
                                                    calling_method = 'likelihood')
      
      
      
      edgeR_data = f_get_epi_from_paird_database(cell_loc, cell_loc, tf, f_p("%s/%s/%s/",base_dir, cell_loc, db_version), dp_threshold = 10, fdr,
                                                      cell_filter = cell_loc, het_filter=c("het"), add_het = TRUE,
                                                      db_data = guest_data_raw, rep=loc_rep, labs= lab_list, target_lab = NA, diffpeak_method = diffpeak_method,
                                                      calling_method = 'edgeR')
      
      binomial_data=f_add_tf_data_snp_class(toupper(tf), binomial_data)
      
      
      dim(binomial_data)
      dim(guest_data_raw)
      
      replicate_pvalue=f_replicates_betabinomial(binomial_data, guest_data_raw, tf,distance = -1)
      #replicateD35_pvalue=f_replicates_betabinomial(binomial_data, guest_data_raw, tf, distance = 35)
      
      betabinomial_pvalue = f_beta_binomial_call(binomial_data)
      
      dim(binomial_data)
      all_data = binomial_data
      all_data$tf_ratio = all_data$ref_tf_dp/(all_data$ref_tf_dp + all_data$alt_tf_dp)
      all_data$total_dp = all_data$ref_tf_dp + all_data$alt_tf_dp
      
      
      all_data$pvalue_binom = (binomial_data$p_value)
      all_data$pvalue_likelihood = (likelihood_data$p_value)
      all_data$pvalue_edgeR = (edgeR_data$p_value)
      all_data$pvalue_betaRep = replicate_pvalue
      all_data$pvalue_betaPool = betabinomial_pvalue$pvalue
      
      tf_correlation=(f_get_significant_correlation_coeff(all_data))
      colnames(cor_matrix) = colnames(tf_correlation)
      rownames(cor_matrix) = rownames(tf_correlation)
      cor_matrix = cor_matrix + tf_correlation
        
      all_data$ASB_edgeR = edgeR_data$ASB
      all_data$ASB_likelihood = likelihood_data$ASB
      all_data$ASB_binom = binomial_data$ASB

      #all_data$ASB_likelihood5e3 = f_fdr_ASB(all_data$pvalue_likelihood, all_data$tf_ratio, fdr = 0.005 )
      #all_data$ASB_likelihood1e2 = f_fdr_ASB(all_data$pvalue_likelihood, all_data$tf_ratio, fdr = 0.01 )
      #all_data$ASB_likelihood5e2 = f_fdr_ASB(all_data$pvalue_likelihood, all_data$tf_ratio, fdr = 0.05 )
      #all_data$ASB_likelihood1e1 = f_fdr_ASB(all_data$pvalue_likelihood, all_data$tf_ratio, fdr = 0.1 )
      all_data$ASB_betaPool = f_fdr_ASB(betabinomial_pvalue$pvalue, all_data$tf_ratio)
      all_data$ASB_betaRep = f_fdr_ASB(replicate_pvalue, all_data$tf_ratio)
      #all_data$ASB_betaRepD35 = f_fdr_ASB(replicateD35_pvalue, all_data$tf_ratio)
      
      ASB_cols=grep('ASB_', colnames(all_data), value = T)
      ASB_collection = rbind(ASB_collection, all_data[,ASB_cols])
      
      table(all_data$ASB_betaPool)
      table(all_data$ASB_betaRep)
      table(all_data$ASB_likelihood)

      

      f_select_highest_target_lab_and_rep(db_data = guest_data_raw, tf = 'dnase', 'broad')
      
      all_data$dnase_ratio=(all_data$ref_dnase_dp + 0.01)/(all_data$alt_dnase_dp + all_data$ref_dnase_dp + 0.02)
      #all_data$dnase_ratio[is.na(all_data$dnase_ratio)] = 0.5
      
      
      subset_all=all_data[!is.infinite(all_data$pvalue_likelihood),]
      
      methods = c('binom', 'edgeR', 'likelihood', 'betaPool', 'betaRep', 'likelihood5e3', 'likelihood1e2', 'likelihood5e2', 'likelihood1e1')

      methods = c('binom', 'edgeR', 'likelihood', 'betaPool', 'betaRep')
      
      #debug(f_method_metrics)
      #Correlation with the DHS
      tmp_vector = sapply(X = c('binom', 'edgeR', 'likelihood', 'betaPool', 'betaRep'), f_method_metrics, subset_all)
      
      stats_collect = rbind( stats_collect, c(tf, cell_loc, nrow(subset_all), as.vector(t(tmp_vector)) ))
      
    }
    
  }

source('./s_library.R')





########P-value correlation###################
#Figure 1A: Heatmap
approach_official_names = c('Binom+Pool', 'Binom+Rep', 'edgeR', 'Beta+Pool', 'Beta+Rep')
colnames(cor_matrix) = c('Binom+Pool', 'Binom+Rep', 'edgeR', 'Beta+Pool', 'Beta+Rep')
rownames(cor_matrix) = c('Binom+Pool', 'Binom+Rep', 'edgeR', 'Beta+Pool', 'Beta+Rep')

tiff(file = f_p("%s/correlation.tiff", fig_directory),res = 310, width = 7, height = 7, units='in')
par(0,2,0,0)
heatmap.2(as.matrix(cor_matrix/cor_matrix[1,1]),scale = 'none', trace = 'none',
          main = 'Correlation Coefficient', margins = c(9, 9), 
          col = colorRampPalette(brewer.pal(9, "OrRd"))(100))
dev.off()


f_read_tiff_as_raster <- function(img_file){
  library(tiff)
  img <- readPNG(img_file)
  rasterGrob(img)
}

cor_matrix_percentage=cor_matrix/cor_matrix[1,1]

average_cor_within_binom_family=(sum(cor_matrix_percentage[c(1,2,4,5),c(1,2,4,5)]) - 4)/12
average_cor_with_edgeR=(sum(cor_matrix_percentage[3,]) - 1)/4


#raster::merge(venn_fig, cor_heatmap, filename = f_p("%s/fig1.tiff", fig_directory) )


#Figure 1B: Venn Diagram

#install.packages('venneuler')
library(venneuler)

plot.VennDiagram2 <- function(x, col, col.fn = function(col) hcl(col * 360, 130, 60), alpha=0.3, main=NULL, edges=200, border=NA, col.txt=1, ...) {
  # calculate total extents
  xtp <- x$centers + x$diameters / 2
  xtm <- x$centers - x$diameters / 2
  xr <- range(c(xtp[,1], xtm[,1]))
  yr <- range(c(xtp[,2], xtm[,2]))
  # create canvas
  plot.new()
  plot.window(xr-0.1, yr, "", asp = 1)
  # adjust alpha for all colors if specified
  n <- length(x$diameters)
  if (missing(col)) col <- col.fn(x$colors)
  if (length(col) < n) col <- rep(col, length.out=n)
  if (!is.na(alpha)) {
    col <- col2rgb(col) / 255
    col <- rgb(col[1,], col[2,], col[3,], alpha)
  }
  # prepare circle coordinates
  s <- seq.int(edges) / edges * 2 * pi
  sx <- cos(s) / 2 # VD uses diameter, not radius
  sy <- sin(s) / 2
  if (!is.null(border)) border <- rep(border, length.out=n)
  # plot all circles
  for (i in seq.int(n))
    polygon(x$centers[i, 1] +  x$diameters[i] * sx, x$centers[i, 2] + x$diameters[i] * sy, col = col[i], border = border[i])
  # if col.txt is not NA, plot the circle text
  if (!all(is.na(col.txt))) text(x$centers, x$labels, col=col.txt)
  # finish with title
  title(main = main, ...)
  invisible(NULL)
}


edgeR_ratio = sum(rowSums(ASB_collection == 'ASB') == 5) / sum(ASB_collection$ASB_edgeR == 'ASB')

likelihood_ratio = sum(ASB_collection$ASB_likelihood == 'ASB') / sum(rowSums(ASB_collection == 'ASB') > 0)


cat('Total heterozygous sites:', nrow(ASB_collection), '\n')

vennObj = venneuler(ASB_collection[,c('ASB_betaPool','ASB_binom','ASB_edgeR','ASB_betaRep','ASB_likelihood')] == 'ASB')
vennObj$labels = rep('',5)

plot.VennDiagram2(vennObj)

ASB_total= colSums(ASB_collection == 'ASB')
vennObj$centers[,'x'] = vennObj$centers[,'x'] + 0.1



tiff(file = f_p("%s/venn.tiff", fig_directory),res = 310, width = 8, height = 7, units='in')
plot.VennDiagram2(vennObj)
par(mai = rep(0,4))
mov = 0.10
font_size = 1.5
text(x = 0.5548763 + mov, y = 0.5240646,  labels = f_p('edgeR\n(%s)', ASB_total['ASB_edgeR']), cex = font_size )

text(x = 0.15 + mov , y = 0.30,  labels = f_p('Binom+Rep\n(%s)', ASB_total['ASB_likelihood'] ), cex = font_size)
lines(c(0.23,0.33)+ mov, y = c(0.30,0.30))

text(x = 0.16 + mov , y = 0.68,  labels = f_p('Binom+Pool\n(%s)', ASB_total['ASB_binom'] ), cex = font_size)
lines(c(0.23,0.33) + mov, y = c(0.68,0.68))

text(x = 0.15 + mov, y = 0.55,  labels = f_p('Beta+Pool\n(%s)', ASB_total['ASB_betaPool'] ), cex = font_size)
lines(c(0.23,0.28) + mov, y = c(0.55,0.55))

text(x = 0.15 + mov, y = 0.42,  labels = f_p('Beta+Rep\n(%s)', ASB_total['ASB_betaRep'] ), cex = font_size)
lines(c(0.23,0.33)+ mov, y = c(0.42,0.42))
dev.off()





f_read_tiff_as_raster <- function(img_file){
  library(tiff)
  img <- readTIFF(img_file)
  rasterGrob(img)
}

cor_heatmap <- f_read_tiff_as_raster(f_p("%s/correlation.tiff", fig_directory))
venn_fig <- f_read_tiff_as_raster(f_p("%s/venn.tiff", fig_directory))


str(venn_fig)
str(cor_heatmap)

tiff(f_p('%s/figure1.tiff', fig_directory), width = 14, height = 7, units = 'in',res = 310)
grid.arrange(arrangeGrob(venn_fig, cor_heatmap, ncol =2))
grid.text(label = '(A)',x=unit(0.05, "npc"), y=unit(0.95, "npc"), gp=gpar(fontsize=16))
grid.text(label = '(B)',x=unit(0.46, "npc"), y=unit(0.95, "npc"), gp=gpar(fontsize=16))
dev.off()


######
f_method_metrics(method = 'betaRep', subset_all)


colnames(stats_collect) = c('tf', 'cell', 'total' , paste0('dnase_', methods), paste0('pvalue_', methods), paste0('estimate_', methods), paste0('ASB_', methods) )


table(subset_all$ASB_binom)


input_stats = stats_collect

f_change_coefficient_based_on_pvalue <- function(input_stats, method){
  fdr_col = p.adjust(input_stats[[f_p('pvalue_%s', method)]], method = 'BH')
  input_stats[fdr_col>0.05,f_p('dnase_%s', method)] = NA
  cat(sum(fdr_col <= 0.05), 'significant correlations in', method)
  return (input_stats)
}


stats_collect_bak = stats_collect
pvalue_list = rep(1, 5)
names(pvalue_list) = methods
edgeR_pvalue = pvalue_list
for (loc_method in methods){
  stats_collect = f_change_coefficient_based_on_pvalue(stats_collect, loc_method)
  pvalue_list[loc_method]=wilcox.test(x = as.numeric(stats_collect$dnase_binom), 
                                      y = as.numeric(stats_collect[,f_p('dnase_%s', loc_method)]),
                                      #mu = 0.01,
                                      alternative = 'greater',
                                      paired = T)$p.value
  
  edgeR_pvalue[loc_method]=wilcox.test(x = as.numeric(stats_collect$dnase_edgeR), 
                                      y = as.numeric(stats_collect[,f_p('dnase_%s', loc_method)]),
                                      #mu = 0.01,
                                      #alternative = 'greater',
                                      paired = T)$p.value
}



library(tidyr)

##Need add the p-value between each other.



all_data$dnase_total = all_data$ref_dnase_dp + all_data$alt_dnase_dp

subset_data <-all_data %>% filter(ref_tf_dp + alt_tf_dp < 20, ref_dnase_dp + alt_dnase_dp > 10)
dnase_dot_data<- subset_data %>% select(dnase_ratio, pvalue_binom, pvalue_edgeR) %>% gather(class, value, -dnase_ratio)

coefs = c(cor.test(subset_data$dnase_ratio, y = -log(subset_data$pvalue_edgeR))$estimate,
          cor.test(subset_data$dnase_ratio, y = -log(subset_data$pvalue_binom))$estimate
)
levels(dnase_dot_data$class) = approach_official_names[c(1,3)]
dot_annotation_data = data.frame( class = approach_official_names[c(3,1)],
                              coef = paste0('R: ',f_p('%.2f', coefs )))




cor_plot<-ggplot(dnase_dot_data, aes(dnase_ratio, -log10(value))) + geom_point(color = '#386cb0', size = 1)+ geom_smooth(color = 'black', lwd = 0.3) +
  facet_wrap(~class) + theme_Publication(12) + geom_text(data = dot_annotation_data, aes(0.20, 13, label = coef ), size = 4) +
  xlab('DHS allelic imbalance on the favored allele') + ylab('-log(p-value)') + theme_axis
cor_plot








colnames(stats_collect)
stats_collect_complete=stats_collect[complete.cases(stats_collect),]
dnase_cor <- stats_collect_complete %>% select(tf, cell, matches('dnase_.*')) %>% gather(class, value,dnase_binom:dnase_betaRep)


dnase_cor$binom = rep(stats_collect_complete$dnase_binom, times = 5)

ggplot(dnase_cor, aes(class, as.numeric(value))) + geom_boxplot()


other_methods = names(pvalue_list[2:5])
annotation_data = data.frame( class2 = approach_official_names[c(3,2,4,5)],
                              pvalues = paste0('P-value: ',f_p('%.1e', as.numeric(pvalue_list[other_methods]))))



dnase_cor$class2 = dnase_cor$class

levels(dnase_cor$class2) = approach_official_names[c(1,3,2,4,5)]

#####Second plot#######
pair_cmp_plot<-ggplot(subset(dnase_cor, class != 'dnase_binom'), aes(as.numeric(binom), as.numeric(value))) + 
geom_point(color = '#386cb0') + geom_abline() + facet_wrap(~class2) + theme_bw(base_size = 16) + 
  geom_text(data = annotation_data, aes(0.30, 0.55, label = pvalues ), size = 4) + xlab('DHS correlation of Binom+Pool') +
  ylab('DHS correlation of other approach') + theme_Publication(12) + theme_axis





tiff(f_p('%s/figure2.tiff', fig_directory), width = 7, height = 10, units = 'in',res = 310)
grid.arrange(arrangeGrob(cor_plot,pair_cmp_plot,nrow = 2, heights = c(1,2)))
grid.text(label = '(A)',x=unit(0.02, "npc"), y=unit(0.98, "npc"), gp=gpar(fontsize=12))
grid.text(label = '(B)',x=unit(0.02, "npc"), y=unit(0.66, "npc"), gp=gpar(fontsize=12))
dev.off()


ggplot(subset(dnase_cor, class != 'dnase_binom'), aes(as.numeric(binom), as.numeric(value))) + theme_Publication()

#fdr_data <- stats_collect %>% select(tf, cell, matches('dnase_likelihood')) %>% gather(class, value,dnase_likelihood5e3:dnase_likelihood1e1)
#ggplot(fdr_data, aes(class, as.numeric(value))) + geom_boxplot()

#pwm_data <- stats_collect %>% select(tf, cell, matches('estimate_')) %>% gather(class, value,estimate_likelihood5e3:estimate_likelihood1e1)
#ggplot(pwm_data, aes(class, as.numeric(value))) + geom_boxplot()


###PWM score#######

wilcox.test(x = as.numeric(stats_collect$estimate_binom), y = as.numeric(stats_collect$estimate_likelihood),paired = T)
wilcox.test(x = as.numeric(stats_collect$estimate_binom), y = as.numeric(stats_collect$estimate_edgeR),paired = T)
wilcox.test(x = as.numeric(stats_collect$estimate_binom), y = as.numeric(stats_collect$estimate_betaRep),paired = T)
wilcox.test(x = as.numeric(stats_collect$estimate_binom), y = as.numeric(stats_collect$estimate_betaPool),paired = T)
wilcox.test(x = as.numeric(stats_collect$estimate_binom), y = as.numeric(stats_collect$estimate_betaRepD35),paired = T)


ggplot(stats_collect, aes(as.numeric(estimate_binom), as.numeric(estimate_likelihood))) + geom_point() + geom_abline()
ggplot(stats_collect, aes(as.numeric(estimate_binom), as.numeric(estimate_edgeR))) + geom_point() + geom_abline()
ggplot(stats_collect, aes(as.numeric(estimate_binom), as.numeric(estimate_betaRep))) + geom_point() + geom_abline()
ggplot(stats_collect, aes(as.numeric(estimate_binom), as.numeric(estimate_betaPool))) + geom_point() + geom_abline()





my_ggplot_setting <- geom_point(size=2) + #, colour = '#00BFC4') + 
  theme(plot.title = element_text(size = 20), 
        axis.title=element_text(size=14,face="bold"), 
        strip.text.x = element_text(size = 14) )



raw_stats = read.table(file = f_p('./result/s_database_count/%s.%s.txt', guest_cell, db_version), sep ='\t', header = TRUE)
tf_list_pwm = subset(raw_stats,  PWM == TRUE)$tf
tf_list   = raw_stats$tf
tf_data_list = cell_data_list[[guest_cell]]
