#' Merge all the PWM data in the predicted TFBS.
#' Makde the classical plot.
source('../s_function.R', chdir = T)
source('../s_asb_manuscript_comotif_func.R', chdir = T)
load(file = '../data/tmp/cell_data_list')
merge_pwm_data = data.frame()

pwm_information_content = read.table('../data/server/snv/pwm_information_content.list',header = T)
rownames(pwm_information_content) = tolower(pwm_information_content$tf)

#str(cell_data_list)

for (loc_cell in names(cell_data_list)){
  
  tf_data_list = cell_data_list[[loc_cell]]
  
  for (loc_tf in names(tf_data_list)){
    
    tf_data = subset(tf_data_list[[loc_tf]], het_type == 'het')
    ##Add the information content
    pwm_len = pwm_information_content[loc_tf, 'length']
    
    if(!is.na(pwm_len) & !is.null(tf_data$pwm_ref_start)){
      
      information_content=as.numeric(unlist(str_split(string = pwm_information_content[loc_tf, 'ic'], pattern = ',')))
      forward_strand = tf_data$pwm_ref_strand == "+"
      dim(tf_data)
      tf_data$pwm_position = tf_data$pwm_ref_start
      tf_data$pwm_position[forward_strand]=pwm_len - tf_data$pwm_ref_start[forward_strand] -1 #1 based
      tf_data$information_content = information_content[tf_data$pwm_position + 1]
      
    }else{
      cat('Missing PWM of', loc_tf)
      tf_data$pwm_position = 0
      tf_data$information_content = 0
    }
    
    
    tf_data$pvalue_ratio_tf = tf_data[[paste0('pvalue_ratio_',loc_tf)]]
    
    tf_pwm_data = tf_data[,grep(pattern = '(.*pwm.*score|ASB|best_match|pwm_position|information_content|p_value|pvalue_ratio_tf$)', colnames(tf_data), value = T)] 
    
    if( !"pvalue_ratio_tf" %in% colnames(tf_pwm_data) ) tf_pwm_data$pvalue_ratio_tf = 0
    
    
    tf_pwm_data$tf = loc_tf
    tf_pwm_data$cell = loc_cell
    colnames(tf_pwm_data)
    tf_pwm_data=f_add_tf_data_snp_class(loc_tf, tf_pwm_data)
    
    colnames(merge_pwm_data)
    colnames(tf_pwm_data)
    
    if (all(colnames(merge_pwm_data) %in% colnames(tf_pwm_data))){
      
      merge_pwm_data = rbind(merge_pwm_data, tf_pwm_data)
    }else{
      
      flog.warn(f_p('Error of the data: %s',loc_tf))
    }
    
    
    
  }
  
  
}





################PWM plot and p-values###########################

dim(merge_pwm_data)
table(merge_pwm_data$class)
classI_data = subset(merge_pwm_data, class == 'class1')
ASB_data = subset(classI_data, ASB == 'ASB')
nonASB_data = subset(classI_data, ASB != 'ASB')

classI_data$ASB[classI_data$ASB != 'ASB'] = 'non-ASB'

dim(ASB_data)

p_value=binom.test(sum(ASB_data$pwm_ref_score > ASB_data$pwm_alt_score), n = nrow(ASB_data))$p.value
cat('Significance of PWM_ref > PWM_alt', p_value, '\n')

#p_value=binom.test(sum(ASB_data$pwm_ref_score > ASB_data$pwm_alt_score), n = nrow(ASB_data), p = sum(nonASB_data$pwm_ref_score > nonASB_data$pwm_alt_score)/nrow(nonASB_data))$p.value

library(ggplot2)
head(classI_data)
p1 <- ggplot(data = classI_data, aes(x = pwm_ref_score, y = pwm_alt_score, color=ASB)) + xlab('Motif score of favored allele') + 
  ylab('Motif score of unfavored allele') +
  geom_point(alpha = I(0.25), size = 1) + geom_abline(size = 0.3) + facet_grid(.~ASB) + coord_fixed()#+ geom_text(data = NULL, x = 0.75, y = 1, label = f_p('P-value = %.3e', p_value))

p1 + theme(legend.position = "none") 
#p1 + geom_density2d(colour="green")
#p1 + stat_binhex(bins = 50)

pp1 <-p1 + theme_Publication(base_size = 9) + scale_colour_Publication() + theme(legend.position = "none")
ggsave(filename = f_p('%s/pwm_classI_merged_TFs.tiff', writing_directory ),plot = pp1, width = 178, height = 101, units = 'mm')



#Save the figure for all the PWM data.
pwm_data = subset(x = merge_pwm_data, pwm_ref_score > 0)
pwm_data$ASB[pwm_data$ASB != 'ASB'] = 'non-ASB'
p1 <- ggplot(data = pwm_data, aes(x = pwm_ref_score, y = pwm_alt_score, color=ASB)) + xlab('Motif score of favored allele') + ylab('Motif score of unfavored allele') +
  geom_point(alpha = I(0.1), size = 1 ) + geom_abline() + facet_grid(.~ASB) + coord_fixed() + guides(colour = guide_legend(override.aes = list(alpha = 1))) + theme_set(theme_gray(base_size = 25))#+ geom_text(data = NULL, x = 0.75, y = 1, label = f_p('P-value = %.3e', p_value))

p1 + guides(color=FALSE)
p1 + theme_Publication() + scale_colour_Publication() + theme(legend.position = "none")
#p1 + theme(legend.position = "none") 
#p1 + geom_density2d(colour="green")
#p1 + stat_binhex(bins = 50)
ggsave(filename = f_p('%s/pwm_merged_TFs.tiff', writing_directory ),width = 178, height = 101, units = 'mm')


#hexbinplot(pwm_alt_score ~ pwm_ref_score | ASB, pwm_data)
#check ASB is enriched in classI sites.

#############PWM ratio plot vs TF binding imbalance and best pwm score########################
dim(classI_data)
colnames(classI_data)
(classI_data$pwm_ref_score - classI_data$pwm_alt_score) > 0

f_plot_pwm_change_against_pwm_rank <- function(pwm_data,ratio_type = 'p_value', window_size = 100){
  library(zoo)
  if(ratio_type == 'p_value'){
    pwm_data$pwm_change = pwm_data$pvalue_ratio
  }else{
    pwm_data$pwm_change = pwm_data$pwm_ref_score/pwm_data$pwm_alt_score
    
  }
  
  
  
  pwm_data_sorted = f_sort_by_col(pwm_data,index = 'p_value')
  head(pwm_data_sorted)
  window100 = rollmean(pwm_data$pwm_change,k = window_size)
  length(window100)
  head(window100)
  nrow(pwm_data)
  plot_data = pwm_data_sorted[(window_size/2):(length(window100) + window_size/2 -1 ),]
  plot_data$window = window100
  
  
  p1<- ggplot(plot_data, aes(x=-log10(p_value), y= window )) + geom_point() + #xlim(c(0,30)) +
    xlab('TF binding imbalance of two alleles( -log10(p-value))') +
    ylab('PWM ratio of favor over unfavor allele')
  
  ggsave(filename = f_p('%s/PWM_ratio(%s)_vs_TF_imbalance.png', writing_directory , ratio_type), width = 10, height = 10)
  
  print(p1)
  
  pwm_data_sorted = f_sort_by_col(pwm_data,index = 'pwm_ref_score')
  window100 = rollmean(pwm_data_sorted$pwm_ref_score/pwm_data_sorted$pwm_alt_score,k = window_size)
  plot_data = pwm_data_sorted[(window_size/2):(length(window100) + window_size/2 -1 ),]
  plot_data$window = window100
  p2 <- ggplot(plot_data, aes(x=pwm_ref_score, y= window )) + geom_point() + 
    xlab('Best PWM score of favored alleles') +
    ylab('PWM ratio of favor over unfavor allele')
  
  
  ggsave(filename = f_p('%s/PWM_ratio(%s)_vs_best_pwm.png',  writing_directory, ratio_type ), width = 10, height = 10)
  
  print(p2)
  
  
}

f_plot_pwm_change_against_pwm_rank(pwm_data, ratio_type = 'p_value')
f_plot_pwm_change_against_pwm_rank(pwm_data, ratio_type = 'pwm_relative', window_size = 500)
#install.packages('zoo')



####################PWM position and Infomation content##############################
#undebug(f_proportion)

#debug(f_position_enrichment_in_pwm)

f_position_enrichment_in_pwm <- function(classI_data, max_flag =F, fix_label = F){
  library(dplyr)
  merge_positions = data.frame()
  for(loc_cell in unique(classI_data$cell)){
    
    cell_data = subset(classI_data, cell == loc_cell)
    
    for(loc_tf in unique(cell_data$tf)){
      tf_data = subset(cell_data, tf == loc_tf)
      ASB_data = subset(tf_data, ASB == 'ASB')
      nonASB_data = subset(tf_data, ASB != 'ASB')
      if(nrow(ASB_data) == 0) next
      
      if(max_flag == TRUE){
        ASB_max = max(table(ASB_data$pwm_position))
        nonASB_max = max(table(nonASB_data$pwm_position))
      }else{
        ASB_max = nrow(ASB_data)
        nonASB_max = nrow(nonASB_data)
      }
      
      
      ASB_dist    <- ASB_data %>% group_by(cell, tf, pwm_position) %>% dplyr::summarise(asb_freq = length(cell)/ASB_max, count = ASB_max, information_content = mean(information_content))
      
      nonASB_dist <- nonASB_data %>% group_by(cell, tf, pwm_position) %>% dplyr::summarise(asb_freq = length(cell)/nonASB_max, count = nonASB_max, information_content = mean(information_content))
      
      merge_dist = merge(ASB_dist, nonASB_dist, by = c('cell','tf', 'pwm_position', 'information_content'),all = T)
      merge_dist[is.na(merge_dist)] = 0
      
      merge_dist$ASB_count  = nrow(ASB_data)
      merge_dist$nonASB_count = nrow(nonASB_data)
      merge_dist$freq_diff = merge_dist$asb_freq.x - merge_dist$asb_freq.y
      merge_positions = rbind(merge_positions, merge_dist)
      
    }
    
  }
  
  library(dplyr)
  
  merge_positions_raw = merge_positions
  
  
  merge_positions <- merge_positions_raw %>% filter(ASB_count >= 10, nonASB_count >= 10)
  
  merge_positions %>% filter(information_content > 1.75, freq_diff < -0.1, ASB_count >10, nonASB_count > 10)
  
  merge_positions %>% filter(information_content > 1, freq_diff > 0.10, ASB_count >10, nonASB_count > 10)
  
  merge_positions %>% filter(cell == 'helas3', tf== 'cebpb')
  merge_positions_raw %>% filter(cell == 'helas3', tf== 'max')
  
  
  merge_positions$label = paste0(toupper(merge_positions$tf),'(', (merge_positions$pwm_position + 1),')')
  
  #low_ic_high_diff = (merge_positions$information_content < 0.5) & (merge_positions$freq_diff > 0.15) | 
  #  (merge_positions$information_content > 1) & (merge_positions$freq_diff < -0.15) | 
  #  merge_positions$freq_diff > 0.15 | merge_positions$freq_diff < -0.2 |
  #  merge_positions$information_content > 1.9 & merge_positions$tf == 'cebpb'
  
  low_ic_high_diff = (merge_positions$information_content > 1.9 & (merge_positions$tf == 'cebpb')) #| merge_positions$freq_diff > 0.15
  
  #if (max_flag == T) low_ic_high_diff =  merge_positions$freq_diff > 0.7# make a false array
  
  merge_positions$label[!low_ic_high_diff] = '' 
  
  
  cor_obj=cor.test(x = merge_positions$freq_diff, merge_positions$information_content, method = 'spearman')
  lm_obj=summary(lm(freq_diff ~ information_content, merge_positions))
  pvalue = lm_obj$coefficients[2,4]
  r.squared = lm_obj$r.squared
  
  pvalue = cor_obj$p.value
  r.squared = cor_obj$estimate
  merge_positions %>% filter(information_content>1.9, freq_diff > 0.18 | freq_diff < 0.05)
  
  #install.packages("FField", type = "source")
  #install.packages("ggplot2")
  #install.packages("gridExtra")
  library(FField)
  label_data = merge_positions[low_ic_high_diff,]
  dim(label_data)
  x.fact <- 100 / max(label_data$information_content)
  y.fact <- 100 / max(label_data$freq_diff)
  
  # Repel points
  coords <-
    FFieldPtRep(coords = cbind(label_data$information_content * x.fact, 
                               label_data$freq_diff * y.fact),
                rep.fact = 40)
  
  # Convert back to plot coordinates
  if(fix_label == T){
    label_data$x.t <- coords$x / x.fact - 0.05
    label_data$y.t <- coords$y / y.fact
  }else{
    
    label_data$x.t = rep(1.949849, 4)
    label_data$y.t = c(0.07, 0.05, 0.18040293, 0.0) 
    #label_data$x.t = c(2.05, 2.06, 2.049849, 2.055963, 2.049450, 2.049450, 2.049450, 2.06, 2.048696)
    #label_data$y.t = c( 0.06, 0.045, 0.180402930, 0.014491922, -0.034926471, 0.082720588, 0.141544118, 0.030241566,-0.003881878)
  }

  

  pp<-ggplot(merge_positions, aes(x = information_content, y = freq_diff)) + geom_point(color = '#00BFC4') +
    geom_smooth() + geom_text(data = label_data, aes(x = x.t, y = y.t, label = label), hjust=1, vjust=0, size = 6) +
    geom_segment(data = label_data,xend = label_data$x.t,yend = label_data$y.t) +  ylim(c(-0.2,0.2)) +
    annotate("text", x = 0.1, y = 0.165, label = f_p('R = %.2f', r.squared),hjust = 0,  size = 6) +
    annotate("text", x = 0.1, y = 0.185, label = f_p('P = %.2e', pvalue),hjust = 0,  size = 6) +
    #ylab('Frequency(ASB) - Frequency(non-ASB)') + 
    ylab('Positional impact') + 
    xlab('Information content')
  
  print(pp)
  
  #pp + theme_Publication() 
  
  ggsave(filename = f_p('%s/information_positional_impact_correlation.png', writing_directory ))
  
  #f_ggsave(filename = , single = T, scale = 2)
  
  
  
  library(sp)
  gb <- ggplot_build(pp)
  
  # get the CI data
  p <- gb$data[[2]]
  
  # make a polygon out of it
  poly <- data.frame(
    x=c(p$x[1],    p$x,    p$x[length(p$x)],    rev(p$x)), 
    y=c(p$ymax[1], p$ymin, p$ymax[length(p$x)], rev(p$ymax))
  )
  
  # test for original values in said polygon and add that to orig data
  # so we can color by it
  merge_positions$in_ci <- point.in.polygon(merge_positions$information_content, merge_positions$freq_diff, poly$x, poly$y)
  
  
  return (list(merge_positions=merge_positions, pp))
  
  
}

#undebug(f_position_enrichment_in_pwm)
#mtrace(f_position_enrichment_in_pwm)
return_list = f_position_enrichment_in_pwm(classI_data,max_flag = F)
pp<-last_plot() + theme_Publication() + theme(axis.title = element_text(face = "bold",size = 14))
return_list[[2]]
pwm_infor_stats_total = return_list$merge_positions
sum(pwm_infor_stats_total$in_ci)/nrow(pwm_infor_stats_total)
table(pwm_infor_stats_total$tf)
length(unique(pwm_infor_stats_total$tf))

head(pwm_infor_stats_total)
cebpb = pwm_infor_stats_total %>% filter(tf == 'cebpb')
cebpb$sign = cebpb$freq_diff < 0

cebpb_plot <- ggplot(data = cebpb, aes(x = pwm_position + 1, y = freq_diff, fill = sign)) + geom_bar(stat = 'identity')  +
theme_Publication() + scale_fill_Publication() + ylab('Positional impact') + xlab('Motif Positions') + theme(legend.position = "none") +
  theme(axis.title = element_text(face = "bold",size = 14))

empyt_block <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + theme_white()


position_grid<-arrangeGrob(pp , arrangeGrob(cebpb_plot,empyt_block,heights=c(3,5),ncol=1) , widths = c(8,6), heights=c(7,7), default.units = "in", ncol=2)

print(position_grid)
ggsave(filename = f_p('%s/pwm_position_and_inforamtion_content.png', writing_directory), plot = position_grid,width = 14, height = 8, dpi = 300)




#View(pwm_infor_stats_total)
write.table(x = pwm_infor_stats_total, file = './data/oriol/pwm_infor_stats_total.txt', quote = F, sep = '\t', row.names = F)

table(pwm_infor_stats_total$tf)

#dimerize = set(["ap2gamma", "batf", "bhlhe40", "cebpb", "ebf1", "jund", "mafk", "max", "mef2a", "myc", "srf", "tcf12", "usf1", "usf2"])

dimer_tfs = c('cebpb', "ebf1", "mef2a",  'max', 'srf', 'tcf12','usf1')

ggplot(subset(pwm_infor_stats_total, tf %in% dimer_tfs), aes(x=pwm_position, y =freq_diff)) + geom_bar(stat = 'identity') +
  facet_wrap(~tf,ncol = 4) + theme_set(theme_gray(base_size = 20))


ggplot(pwm_infor_stats_total, aes(x=pwm_position, y =freq_diff)) + geom_bar(stat = 'identity') +
  facet_wrap(~tf,ncol = 4) + theme_set(theme_gray(base_size = 20))

length(unique(pwm_infor_stats_total$tf))

#pwm_infor_stats_max = f_position_enrichment_in_pwm(classI_data,max_flag = T)
#write.table(x = pwm_infor_stats_max, file = './data/oriol/pwm_infor_stats_max.txt', quote = F, sep = '\t', row.names = F)
#View(pwm_infor_stats_max)


############Dimmer stats test###############


dimer_tf_data = subset(pwm_infor_stats_total, tf %in% dimer_tfs)
#dimer_positions = data.frame(tf = c('cebpb', "ebf1", 'max', "mef2a",   'srf', 'tcf12','usf1'),
#                             row.names = c('cebpb', "ebf1", 'max', "mef2a", 'srf', 'tcf12','usf1'),
#                             central_position = c(5, 5, 5, 7, 7, 6, 5))

dimer_positions = data.frame(tf = c('cebpb', 'max', 'tcf12','usf1'),
                             row.names = c('cebpb', 'max', 'tcf12','usf1'),
                             central_position = c(5, 5, 6, 5))


loc_tf = 'cebpb'

dimer_tf_data$class = 'no_select'
dimer_tf_data$rank = 0
rownames(dimer_tf_data) = paste(dimer_tf_data$tf, dimer_tf_data$pwm_position, sep = '-')
for ( loc_tf in dimer_positions$tf){
  
  dimer_tf_data[dimer_tf_data$tf==loc_tf,'rank'] =  rank(dimer_tf_data[dimer_tf_data$tf==loc_tf,'freq_diff'])
  
  
  
  row_name = f_p('%s-%s', loc_tf, dimer_positions[loc_tf,'central_position'] )
  central_pos_ic = dimer_tf_data[row_name, 'information_content']
  similar_ic_positions = dimer_tf_data$tf == loc_tf & dimer_tf_data$information_content >= 0.9 *central_pos_ic & dimer_tf_data$information_content <= 1.1 *central_pos_ic
  dimer_tf_data$class[similar_ic_positions]= 'selected'
  dimer_tf_data[row_name, 'class'] = 'center'  
}

dimer_test = wilcox.test(x = unlist(subset(dimer_tf_data, class == 'selected', select = rank)), alternative = 'less',
            y = unlist(subset(dimer_tf_data, class == 'center', select = rank)),)


cat(f_p('Dimer test p-value %s', dimer_test$p.value))



###########Plot the information content against positional impact ############

#install.packages('latticeExtra')
library(latticeExtra)
pwm_infor_stats_total$pwm_position = as.numeric(pwm_infor_stats_total$pwm_position)


plot_data = arrange(pwm_infor_stats_total, tf, pwm_position)
plot_data$tf = factor(f_rename_words_for_figures(plot_data$tf,upper = T))
plot_data$cell = factor(f_rename_words_for_figures(plot_data$cell,upper = F))
plot_data$tf = paste0(plot_data$tf,"(",plot_data$cell,")")

temp <- xyplot(information_content ~ pwm_position | tf,xlab = 'Moitf Position', ylab = 'Information content',
               data =plot_data , type = "l", )
rain <- xyplot(freq_diff ~ pwm_position | tf, data = plot_data, type = "h",xlab = 'Motif Position', ylab = 'Positional impact')

par(cex = 2)


library(latticeExtra)
old_setting=trellis.par.get('superpose.line')
new_setting = old_setting
new_setting$col = c('red','#6495ED', rep("#6495ED",5))

trellis.par.set('superpose.line',new_setting)

trellis.par.get('superpose.line')

tiff(filename = f_p("%s/positional_inpact_imformation_content_all_tfs.tif", writing_directory),width = 10.5, height = 10.5,units = 'in', res = 300)
doubleYScale(rain, temp, add.ylab2 = TRUE, use.style =T,  text = c( "Positional impact", "Information content"),  columns = 3)
dev.off()
trellis.par.set('superpose.line',old_setting)







