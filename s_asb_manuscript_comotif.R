
source('./test/t_feature_heatmap.r')
source('s_cofactor_search.r')
f_convert_to_compact_table <- function(corelated_motif){
  compact_tf_table = data.frame(tf = unique(corelated_motif$tf), comotif = '',row.names =unique(corelated_motif$tf) )
  
  for (loc_tf in compact_tf_table$tf){
    loc_comotifs = corelated_motif[corelated_motif$tf == loc_tf,'homer_motif']
    #compact_tf_table[loc_tf,'comotif']=paste(unique(str_replace(string = loc_comotifs, pattern = '[.][a-z0-9]*[.]homer', replacement = '' )),collapse = ', ')
    compact_tf_table[loc_tf,'comotif']=paste(loc_comotifs, collapse = ', ')
  }
  
  return (compact_tf_table)
  
}


f_find_selected_lab <- function(tf_data){
  
  lab_cols = grep('lab_.*', colnames(tf_data), value = T)
  lab=names(which(apply(tf_data[lab_cols] , 2, function(x) all(x != '.') )))
  lab = str_replace(string = lab, pattern = 'lab_', replacement = '')
  return (lab)
}


f_show_top_comotifs <- function(loc_cell, cell_data_list){
  
  
  motif_table = read.table(f_p('./data/server/chipseq_snv/%s_homer_top5_motif.txt', loc_cell),sep ='\t',header = T)  
  
  tf_set = unique(motif_table$target_tf)
  
  lab_tf = tf_set[1]
  
  tf_comotif_stats = data.frame()
  
  lab_tf = 'uw-helas3-ctcf'
  for(lab_tf in tf_set){
    
    tf_data = subset(motif_table, target_tf == lab_tf )  
    
    
    class(do.call(rbind.data.frame, str_split(tf_data$motif_name,pattern = '/'))[1])
    
    tf_line = c(rep('',times = 7))
    for ( i in 1:nrow(tf_data) )
    {
      tf_line[i + 1] = paste(str_split(tf_data$motif_name[i], pattern = '/')[[1]][1], tf_data$target_seq_percent[i], sep = '-')
    }
    
    tf_line[1] = lab_tf
    
    tf_name = str_split(lab_tf, pattern = '-')[[1]][3]
    
    if (length(grep(tf_name, tf_line[-1], ignore.case = T)) == 0){
      tf_line[7] = 'NA'
    }else{
      tf_line[7] =  tf_line[grep(tf_name, tf_line[-1], ignore.case = T)[1] + 1]
    }
    
    
    tf_comotif_stats=rbind(tf_comotif_stats, tf_line)
    
  }
  
  dim(tf_comotif_stats)
  colnames(tf_comotif_stats) = c('tf', 'top1', 'top2', 'top3', 'top4', 'top5','itself')
  
  rownames(tf_comotif_stats) = tf_comotif_stats$tf
  
  head(tf_comotif_stats)
  
  tf_name = 'ctcf'
  
  selected_stats = data.frame()
  for (tf_name in names(cell_data_list[[loc_cell]])){
    
    tf_data = cell_data_list[[loc_cell]][[tf_name]]
    
    #Find the which lab the data is using
    lab_cols = grep('lab_.*', colnames(tf_data), value = T)
    
    lab=names(which(apply(tf_data[lab_cols] , 2, function(x) all(x != '.') )))
    lab = str_replace(string = lab, pattern = 'lab_', replacement = '')
    
    
    selected_stats = rbind(selected_stats, tf_comotif_stats[paste(lab,loc_cell,tf_name,sep = '-'),])
    
  }
  
  selected_stats$tf = NULL
  
  return (selected_stats)
  
}


writing_directory = 'E:/Projects/Manuscript/ASB Manuscript/data'


#undebug(f_homer_motif)
f_homer_motif <- function(tf_data_list, best_flag = FALSE){
  
  cor_statics = data.frame()
  
  
  for (tf in names(tf_data_list)){
    
    cat(tf)
    guest_data = tf_data_list[[tf]]
    
    het_data = subset(guest_data, het_type == 'het')
    if (best_flag == TRUE & 'pwm_ref_score' %in% colnames(het_data)){
      het_data = subset(het_data, best_match == 'Best')
    }
    
    #print(table(het_data$ASB, het_data$dnase_uw_dis))
    
    
    homer_motifs = grep('pvalue.*homer', colnames(het_data), value = T)
    
    
    threshold_list=f_sort_by_col(f_read_table3('./data/server/snv/pwm_threshold.list',quiet = TRUE), index = "tf")
    rownames(threshold_list) = threshold_list$tf
    threshold_list['PU1',] = threshold_list['SPI1',]
    
    tf_threshold=as.numeric(threshold_list[toupper(tf),'threshold'])
    
    
    het_data$class = 'class3'
    
    if (!is.na(tf_threshold )){    
      het_data$class[het_data$pwm_alt_score > tf_threshold | het_data$pwm_ref_score > tf_threshold] = 'class1'
    }else{
      het_data[,c("pwm_ref_score", "pwm_alt_score", "peak_pwm_ref_score", "peak_pwm_alt_score")] = 0
    }
    
    
    
    
    for (feature in homer_motifs){
      
      binding_sites = het_data[[str_replace(string = feature, pattern = 'pvalue_ratio', replacement = 'binding_pwm')]]
      
      motif_data = het_data[binding_sites & het_data$class != 'class1',]
      #motif_data = het_data[binding_sites ,]
      feature_ratio = as.numeric(motif_data[[feature]])
      ASB_index = (motif_data$ASB == 'ASB')
      
      tf_ratio = (motif_data[[f_p('ref_tf_dp')]] + 1) / (motif_data[[f_p('alt_tf_dp')]] + 1)
      #tf_ratio = (motif_data[[f_p('ref_tf_dp')]]) / (motif_data[[f_p('alt_tf_dp')]] + motif_data[[f_p('ref_tf_dp')]])
      
      tf_ratio_flip = tf_ratio
      tf_ratio_flip[motif_data$alt_flip] = 1/tf_ratio_flip[motif_data$alt_flip]
      #tf_ratio_flip[motif_data$alt_flip] = 1 - tf_ratio_flip[motif_data$alt_flip]
      feature_ratio_flip = feature_ratio
      feature_ratio_flip[motif_data$alt_flip] = -1 * feature_ratio_flip[motif_data$alt_flip]
      
      #if ( nrow(subset(motif_data, ASB == 'ASB')) >= 5 & sum(binding_sites &  het_data$class == 'class1' )/sum(binding_sites) < 0.8){
      if ( nrow(subset(motif_data, ASB == 'ASB')) >= 1 & sum(!ASB_index) > 0){
      #if ( nrow(motif_data) >= 5){
        stats=cor.test(x = as.numeric(feature_ratio), y = tf_ratio, method ='spearman')
        flip_stats=cor.test(x = as.numeric(feature_ratio_flip), y = tf_ratio_flip,method ='spearman')
        
        wilcox_stats = wilcox.test(x = as.numeric(feature_ratio[ASB_index]),y = as.numeric(feature_ratio[!ASB_index]))
        
        if( !is.na(stats$p.value)){
          
          simple_feature_name = str_replace(feature, pattern = '(pvalue_ratio_|[.]homer)', replacement = '')    
          cor_statics = rbind(cor_statics, c(tf, simple_feature_name,  stats$p.value, wilcox_stats$p.value, flip_stats$p.value  ) )
        }
      }
      
      
      
      
      #print(qplot(ref_dnase_dp, alt_dnase_dp, data=het_data, position ='jitter' ,xlab="Allele with stronger ChIP-seq signal",ylab="Allele with weaker ChIP-seq signal",
      #) + ggtitle(f_p('DNase ChIP-seq: %s-%s', guest_cell, tf)) + facet_wrap(~ASB) + geom_abline()  + my_ggplot_setting)
      
    }
    
    
  }
  
  
  colnames(cor_statics) = c('tf', 'homer_motif', 'pavlue', 'wilcox_p.value','flip_p.value')
  cor_statics$fdr = (p.adjust(cor_statics$flip_p.value,method = ))
  return (cor_statics)
  
}


# debug(f_get_snps_in_comotifs)
# mtrace(f_get_snps_in_comotifs)
# undebug(f_get_snps_in_comotifs)
f_get_snps_in_comotifs <- function(compact_tf_table, tf_data_list){
  
  statistic_table = data.frame(tf = names(tf_data_list), peak_pwm_ratio = 0,row.names = names(tf_data_list),
                               peak_pwm_ratio=0, asb_snp_in_comotif_overlap =0,
                               asb_comotif_exclusive = 0, asb_snp_in_tfbs =0,
                               asb_best_ratio = 0, nonASB_snp_in_comotif_overlap =0,
                               nonASB_comotif_exclusive = 0, nonASB_snp_in_tfbs =0,
                               nonASB_best_ratio = 0, ASB_num = 0
  )
  
  
  for (loc_tf in names(tf_data_list)){
    het_data = subset(tf_data_list[[loc_tf]], het_type == 'het')
    
    threshold_list=f_sort_by_col(f_read_table3('./data/server/snv/pwm_threshold.list',quiet = TRUE), index = "tf")
    rownames(threshold_list) = threshold_list$tf
    threshold_list['PU1',] = threshold_list['SPI1',]
    
    tf_threshold=threshold_list[toupper(loc_tf),'threshold']
    
    
    statistic_table[loc_tf, 'ASB_num'] = sum(het_data$ASB == 'ASB')
    statistic_table[loc_tf, 'nonASB_num'] = sum(het_data$ASB != 'ASB')
    if (!is.na(tf_threshold )){
      
      statistic_table[loc_tf, 'peak_pwm_ratio'] = sum(het_data$peak_pwm_alt_score > tf_threshold | het_data$peak_pwm_ref_score > tf_threshold)/nrow(het_data)
      snp_in_tfbs=(het_data$pwm_alt_score > tf_threshold | het_data$pwm_ref_score > tf_threshold)
      
    }else{
      statistic_table[loc_tf, 'peak_pwm_ratio'] = 0
      snp_in_tfbs = (het_data$chr == 'xxx') #create false array
      
    }
    
    
    
    ####Comotifs###########
      
    comotifs = compact_tf_table[loc_tf, 'comotif']
    
    if (!is.na(comotifs)){
      
      other_motifs = grep(loc_tf,unlist(str_split(comotifs, pattern = ', ')),value = T, invert = T)
      
      #select the less overlap motifs
      selected_motifs = c()
      for(loc_motif in other_motifs){
        
        loc_binding_col=grep(f_p('binding_pwm_(%s).*',loc_motif), colnames(het_data), value = T)
        
        if(length(loc_binding_col) > 1){
          flog.warn(f_p('Multiple loc_binidng:%s', paste(loc_binding_col)))
          loc_binding_col = loc_binding_col[1]
        }
        
        if(sum(het_data[[loc_binding_col]] &  snp_in_tfbs )/sum(het_data[[loc_binding_col]]) < 0.8){
          
          selected_motifs = c(selected_motifs, loc_motif)
          
        } 
        
        
      }
      
      
      tfs_or=paste(selected_motifs, collapse = '|')
      
      
      
      binding_cols=grep(f_p('binding_pwm_(%s).*',tfs_or), colnames(het_data), value = T)
      
      
      
      if (tfs_or != ""){
        snp_in_comotif_overlap = rowSums(het_data[binding_cols]) > 0 & snp_in_tfbs == TRUE
        snp_in_comotif_exclusive = rowSums(het_data[binding_cols]) > 0 & snp_in_tfbs == FALSE
        
        statistic_table[loc_tf, 'asb_snp_in_comotif_overlap'] = sum(snp_in_comotif_overlap[het_data$ASB == 'ASB'])/sum(het_data$ASB == 'ASB')
        statistic_table[loc_tf, 'asb_comotif_exclusive'] = sum(snp_in_comotif_exclusive[het_data$ASB == 'ASB'])/sum(het_data$ASB == 'ASB')
        statistic_table[loc_tf, 'nonASB_snp_in_comotif_overlap'] = sum(snp_in_comotif_overlap[het_data$ASB != 'ASB'])/sum(het_data$ASB != 'ASB')
        statistic_table[loc_tf, 'nonASB_comotif_exclusive'] = sum(snp_in_comotif_exclusive[het_data$ASB != 'ASB'])/sum(het_data$ASB != 'ASB')
      }
      
    }
    
    f_get_proportion_TFBS_in_peaks <- function(input_data, tf_threshold){
      
      proportion = sum(input_data$peak_pwm_alt_score > tf_threshold | input_data$peak_pwm_ref_score > tf_threshold)/nrow(input_data)
      return (proportion)
    }
    

    statistic_table[loc_tf, 'asb_snp_in_tfbs'] = sum(snp_in_tfbs[het_data$ASB == 'ASB'])/sum(het_data$ASB == 'ASB') - statistic_table[loc_tf, 'asb_snp_in_comotif_overlap']    
    statistic_table[loc_tf, 'asb_best_ratio'] = f_get_proportion_TFBS_in_peaks(subset(het_data, ASB == 'ASB'), tf_threshold)
    
    statistic_table[loc_tf, 'nonASB_snp_in_tfbs'] = sum(snp_in_tfbs[het_data$ASB != 'ASB'])/sum(het_data$ASB != 'ASB') - statistic_table[loc_tf, 'nonASB_snp_in_comotif_overlap']      
    statistic_table[loc_tf, 'nonASB_best_ratio'] = f_get_proportion_TFBS_in_peaks(subset(het_data, ASB != 'ASB'), tf_threshold)
    
  
  }
    
  
  
  return (statistic_table)
  
}

#mtrace(f_comotif_search)
mtrace.off()
#undebug(f_comotif_search)
#undebug(f_homer_motif)


f_add_tf_data_snp_class <- function(loc_tf, tf_pwm_data){
  
  #class I: SNP in the PWM
  #Class II: SNP in the flanking of PWM
  #Class III: now known PWM 
  
  threshold_list=f_sort_by_col(f_read_table3('./data/server/snv/pwm_threshold.list',quiet = TRUE), index = "tf")
  rownames(threshold_list) = threshold_list$tf
  threshold_list['PU1',] = threshold_list['SPI1',]
  
  tf_threshold=as.numeric(threshold_list[toupper(loc_tf),'threshold'])
  
  flog.info(f_p('%s: %s',toupper(loc_tf), tf_threshold))
  
  tf_pwm_data$class = 'class3'
  
  if (!is.na(tf_threshold )){
    
    tf_pwm_data$class[tf_pwm_data$peak_pwm_alt_score > tf_threshold | tf_pwm_data$peak_pwm_ref_score > tf_threshold] = 'class2'
    
    tf_pwm_data$class[tf_pwm_data$pwm_alt_score > tf_threshold | tf_pwm_data$pwm_ref_score > tf_threshold] = 'class1'
  }else{
    tf_pwm_data[,c("pwm_ref_score", "pwm_alt_score", "peak_pwm_ref_score", "peak_pwm_alt_score")] = 0
  }
  
  
  return (tf_pwm_data)
  
}

#debug(f_comotif_search)
f_comotif_search <- function(cell_data_list){
  #According to the tf data in cell_data_list. Find the correlated motifs with ASB events.
  #Then write a word table and txt table to the writing directory
  
  cell_list = names(cell_data_list)
  compact_tf_stats = data.frame()
  #txt_table = data.frame()
  merge_corelated_motif = data.frame()
  for (loc_cell in cell_list){
    tf_data_list = cell_data_list[[loc_cell]]
    comotif_df = f_show_top_comotifs(loc_cell, cell_data_list)
    f_table_to_word(comotif_df, file_prefix = f_p('comotif_table_%s',loc_cell), data_dir = writing_directory)
    corelated_motif = f_homer_motif( tf_data_list )
    corelated_motif$cell = loc_cell
    merge_corelated_motif = rbind(merge_corelated_motif, corelated_motif)
  }
  
  colnames(merge_corelated_motif)
  merge_corelated_motif$fdr = p.adjust(merge_corelated_motif$flip_p.value)
  for(loc_cell in cell_list ){
    
    tf_data_list = cell_data_list[[loc_cell]]
    #comotif_df = f_show_top_comotifs(loc_cell, cell_data_list)
    
    
    
    compact_tf_table = f_convert_to_compact_table(subset(merge_corelated_motif, fdr < 0.05 & cell == loc_cell  ))
       
    new_statistic_table = f_get_snps_in_comotifs(compact_tf_table, tf_data_list)
    tmp_table = f_sort_by_col(subset(new_statistic_table, ASB_num > 0), index = 'peak_pwm_ratio', decr_flag = T)
    
    selected_columns = c(grep('.*asb_(snp|comotif)', colnames(tmp_table),value = T, ignore.case = T), 'peak_pwm_ratio')
    
    plot_data = data.frame(Ratio = unlist(tmp_table[selected_columns]), 
                           TF = factor(rep(tmp_table$tf, times = 7), levels = tmp_table$tf ),
                           name = rep(selected_columns, each = nrow(tmp_table))
    )
    
    
    
    plot_data$group = str_match(string = plot_data$name, pattern = '(^asb|nonASB|peak)')[,1]
    plot_data$color = str_match(string = plot_data$name, pattern = '(^peak|tfbs|comotif_overlap|comotif_exclusive)')[,1]
    plot_data$color[plot_data$color == 'peak'] = 'tfbs'
    plot_data$group[plot_data$group == 'peak'] = 'Peak TFBS'
    plot_data$group[plot_data$group == 'asb'] = 'ASB'
    plot_data$group = factor(plot_data$group, levels = c('Peak TFBS', 'ASB', 'nonASB'))
    
    
    plot_data$color = factor(plot_data$color, levels = c( 'comotif_exclusive',  'comotif_overlap', 'tfbs'), ordered = T) 
    
    
    ggplot(plot_data, aes(x =group, y=Ratio ,fill = color, order = -as.numeric(color))) + 
      geom_bar(stat="identity",colour="black" ,position = 'stack' ) + ggtitle('Proportion of SNPs in putative TFBS ') +
      facet_wrap(~TF, nrow = 1) + theme(axis.text.x = element_text(angle = 90, vjust = -0.5))
    
    #f_table_to_word(input_data = compact_tf_table, file_prefix = f_p('correlated_comotif_%s', loc_cell), data_dir = writing_directory)
    ggsave(filename = f_p('%s/snp_distributiion_in_tf_motif_and_comotif_%s.png', writing_directory, loc_cell ))
    
    compact_tf_table$cell = loc_cell
    
    compact_tf_stats = rbind(compact_tf_stats, compact_tf_table)
  }
  compact_tf_stats_org = compact_tf_stats
  compact_tf_stats$tf = f_rename_tfs(compact_tf_stats$tf)
  compact_tf_stats$comotif = f_rename_homer_motifs(compact_tf_stats$comotif)
  compact_tf_stats$cell = f_rename_words_for_figures(compact_tf_stats$cell)
  
  compact_tf_stats = rbind(compact_tf_stats, c(nrow(compact_tf_stats), length(unlist(str_split(compact_tf_stats$comotif,pattern = ', '))), 2))
  
  colnames(compact_tf_stats) = c('TF','Comotif','Cell')
  f_table_to_word(input_data = compact_tf_stats, file_prefix = 'correlated_comotif', data_dir = writing_directory)
  
  colnames(compact_tf_stats_org) = c('tf','features','cell')
  write.table(x = compact_tf_stats_org, file = f_p('%s/analytical_results/comotif.txt', writing_directory), quote = F, sep = '\t', row.names = F)
  
  return (compact_tf_stats_org)
  
}




f_rename_tfs <- function(input_vector){
  
  input_vector = toupper(input_vector)
  input_vector = str_replace(input_vector, pattern = 'gabp', replacement = 'gabpa')
  return (input_vector)
}

f_rename_homer_motifs <- function(input_vector){
  input_vector = str_replace_all(input_vector, pattern = '[.]homer', replacement = '')
  input_vector = str_replace(input_vector, pattern = 'nf.e2', replacement = 'nf-e2')
  input_vector = toupper(input_vector)
  return (input_vector)
  
}


#undebug(f_comotif_search2)
f_comotif_search2 <- function(cell_data_list){
  #Similar to the f_comotif_search1
  #Divide the TFs into four classes according to PWM and comotif.
  #Get the enrichment measure of ASB events in each class using fisher test.
  
  cell_list = names(cell_data_list)
  merge_table = data.frame()
  
  for(loc_cell in cell_list ){
    tf_data_list = cell_data_list[[loc_cell]]
    comotif_df = f_show_top_comotifs(loc_cell, cell_data_list)
    corelated_motif = f_homer_motif( tf_data_list )
    compact_tf_table = f_convert_to_compact_table(subset(corelated_motif, fdr <= 0.05 ) )
    new_statistic_table = f_get_snps_in_comotifs(compact_tf_table, tf_data_list)
    
    new_statistic_table = new_statistic_table[complete.cases(new_statistic_table),]
      
    tmp_table = f_sort_by_col(subset(new_statistic_table, ASB_num > 0), index = 'peak_pwm_ratio', decr_flag = T)
    
    selected_columns = c(grep('.*asb_(snp|comotif)', colnames(tmp_table),value = T, ignore.case = T), 'peak_pwm_ratio')
    
    
    tmp_table$cell = loc_cell
    merge_table = rbind(merge_table, tmp_table)
    
  }
    
  
  
  
  
  tmp_table = merge_table
  
  asb_cols = grep('^asb.*', colnames(tmp_table), value = T)
  nonASB_cols =grep('nonASB_num' ,grep('nonASB.*', colnames(tmp_table), value = T), value = T, invert = T)
  tmp_table[,asb_cols] = tmp_table[,asb_cols]*tmp_table$ASB_num
  tmp_table[,nonASB_cols] = tmp_table[,nonASB_cols]*tmp_table$nonASB_num
  tmp_table[,'peak_pwm_ratio'] = tmp_table[,'peak_pwm_ratio']*(tmp_table$nonASB_num + tmp_table$ASB_num)

  #Calculate the ASB events.
  class_I_snv=sum(tmp_table[,'asb_snp_in_tfbs'])/sum(tmp_table$ASB_num) #SNVs in the predicted TFBS
  class_II_snv=sum(tmp_table[,'asb_best_ratio'])/sum(tmp_table$ASB_num) - class_I_snv # SNVs outside the TFBS, but with a TFBS peak
  class_III_snv = 1 -  class_I_snv - class_II_snv #SNVs with no motif peak or TF without motif
  cat('\n', f_p('ASB SNV class I: %.4f, II: %.4f, III:%.4f',class_I_snv, class_II_snv,class_III_snv  ),'\n')
  
  #Calculate the nonASB events.
  class_I_snv=sum(tmp_table[,'nonASB_snp_in_tfbs'])/sum(tmp_table$nonASB_num)
  class_II_snv=sum(tmp_table[,'nonASB_best_ratio'])/sum(tmp_table$nonASB_num) - class_I_snv
  class_III_snv = 1 -  class_I_snv - class_II_snv
  cat(f_p('nonASB SNV class I: %.4f, II: %.4f, III:%.4f',class_I_snv, class_II_snv,class_III_snv  ),'\n')
  
  
  #Calculate addtional SNVs cover ASB events.
  colnames(tmp_table)
  tfs_with_comotif = subset(tmp_table, asb_comotif_exclusive > 0)
  additional_portion = sum(tfs_with_comotif[['asb_comotif_exclusive']])/sum(tfs_with_comotif[['ASB_num']])
  additional_portion_overall = sum(tfs_with_comotif[['asb_comotif_exclusive']])/sum(tmp_table[['ASB_num']])
  cat(f_p('Additional proportion explained by comotif in ASB events: %.4f', additional_portion  ),'\n')
  cat(f_p('Additional proportion explained by comotif in ASB events overall: %.4f', additional_portion_overall  ),'\n')
  
  #This is the TF class based on overlap with comotif.
  tmp_table$class = 'Class4'
  tmp_table$class[ tmp_table$peak_pwm_ratio > 0 & rowSums(tmp_table[,c('asb_snp_in_comotif_overlap', 'asb_comotif_exclusive')]) == 0 ] = 'TFs with PWM only'
  tmp_table$class[ tmp_table$peak_pwm_ratio > 0 & rowSums(tmp_table[,c('asb_snp_in_comotif_overlap', 'asb_comotif_exclusive')]) > 0 ] = 'TFs with PWM and comotif'
  tmp_table$class[ tmp_table$peak_pwm_ratio == 0 & rowSums(tmp_table[,c('asb_snp_in_comotif_overlap', 'asb_comotif_exclusive')]) > 0 ] = 'TFs with comotif only'
  
  loc_class = 'Class1'
  
  
  collect_stats = data.frame()
  for(loc_class in unique(tmp_table$class)){
    
    class_data = subset(tmp_table, class == loc_class)
    numeric_cols=grep('asb|pwm',colnames(class_data), value = T, ignore.case = T)
    class_sum = colSums(class_data[,numeric_cols])
    
    collect_stats = rbind(collect_stats , c(as.numeric(class_sum)))
    
    colnames(collect_stats) = c(names(class_sum))
  }
  
  four_classes = unique(tmp_table$class)
  
  #Fisher test
  for (i in 1:nrow(collect_stats)){
    
    ASB_motif=rowSums(collect_stats[i, asb_cols[1:3]])
    ASB_nonMotif = collect_stats[i, 'ASB_num'] - ASB_motif
    
    nonASB_motif = rowSums(collect_stats[i, nonASB_cols[1:3]])
    nonASB_nonMotif = collect_stats[i, 'nonASB_num'] - nonASB_motif
    
    pvalue=fisher.test(matrix(c(ASB_motif, nonASB_motif, ASB_nonMotif, nonASB_nonMotif),nrow=2,ncol=2),alternative="greater")$p.value
    
    print(f_p('Class%s: %s enrichment with ASB and nonASB:  %s', i,  four_classes[i], pvalue ))
    print(matrix(c(ASB_motif, nonASB_motif, ASB_nonMotif, nonASB_nonMotif),nrow=2,ncol=2))
    cat(ASB_motif/(ASB_motif + ASB_nonMotif),'ASB in motif while', nonASB_motif/(nonASB_motif + nonASB_nonMotif),'for nonASB')
    
  }
  
  collect_stats$Class = unique(tmp_table$class)
  collect_stats[,asb_cols] = collect_stats[,asb_cols]/collect_stats$ASB_num
  collect_stats[,nonASB_cols] = collect_stats[,nonASB_cols]/collect_stats$nonASB_num
  collect_stats[,'peak_pwm_ratio'] = collect_stats[,'peak_pwm_ratio']/(collect_stats$ASB_num + collect_stats$nonASB_num)
  #str(collect_stats)
  
  
  plot_data = data.frame(Ratio = unlist(collect_stats[selected_columns]), 
                           Class = factor(rep(collect_stats$Class, times = length(selected_columns) ), levels = collect_stats$Class ),
                           name = rep(selected_columns, each = nrow(collect_stats))
    )
    
    
    plot_data$group = str_match(string = plot_data$name, pattern = '(^asb|nonASB|peak)')[,1]
    plot_data$color = str_match(string = plot_data$name, pattern = '(^peak|tfbs|comotif_overlap|comotif_exclusive)')[,1]
    plot_data$color[plot_data$color == 'peak'] = 'tfbs'
    plot_data$group[plot_data$group == 'peak'] = 'Peak TFBS'
    plot_data$group[plot_data$group == 'asb'] = 'ASB'
    plot_data$group = factor(plot_data$group, levels = c('Peak TFBS', 'ASB', 'nonASB'))
    
    #plot_data = f_sort_by_col(plot_data, index = 'Position')
    
    plot_data$color = factor(plot_data$color, levels = c( 'comotif_exclusive',  'comotif_overlap', 'tfbs'), ordered = T) 
    plot_data$Position = as.character(plot_data$color) 
    plot_data$Position = str_replace(plot_data$Position, pattern = 'comotif_overlap|tfbs', 'Primary PWM') 
    plot_data$Position = str_replace(plot_data$Position, pattern = 'comotif_exclusive', 'Comotif') 
    
    plot_data$Position = as.factor(plot_data$Position)
    plot_data$Class
  
    ggplot(subset(plot_data, Class != 'Class4' & name != "peak_pwm_ratio"), aes(x =group, y=Ratio ,fill = color)) + 
      geom_bar(stat="identity",colour="black" ,position = 'stack' ) + ggtitle('Proportion of SNPs in putative TFBS ') +
      facet_wrap(~Class, nrow = 1) + theme(axis.text.x = element_text(angle = 90, vjust = -0.5)) + ylab('Fraction of SNPs in TFBS')
    
  sub_data = subset(plot_data, Class != 'Class4' & name != "peak_pwm_ratio")
  sub_data['asb_snp_in_comotif_overlap2','Ratio'] = sub_data['asb_snp_in_comotif_overlap2','Ratio'] + sub_data['asb_snp_in_tfbs2','Ratio']
  sub_data['nonASB_snp_in_comotif_overlap2','Ratio'] = sub_data['nonASB_snp_in_comotif_overlap2','Ratio'] + sub_data['nonASB_snp_in_tfbs2','Ratio']
  sub_data2 = sub_data[!rownames(sub_data) %in%c('asb_snp_in_tfbs2','nonASB_snp_in_tfbs2'),]
  f_sort_by_col(sub_data,'Class')
  
  levels(sub_data2$Class)
  
  sub_data2$Class = factor(sub_data2$Class, levels = c("TFs with comotif only","TFs with PWM and comotif", "TFs with PWM only"))
  sub_data2$Class = str_replace(sub_data2$Class, pattern = "TFs with", replacement = "ChIP-seq with\n")
  #sub_data2$Position = str_replace(sub_data2$Position, pattern = "Primary PWM", replacement = "Primary motif")
  
  pp <- ggplot(sub_data2, aes(x =group, y=Ratio ,fill = Position,  order = -as.numeric(Position))) + xlab('') +
    geom_bar(stat="identity",colour="black" ,position = 'stack' ) + # gtitle('Proportion of SNPs in putative TFBS ') +
    facet_wrap(~Class, nrow = 1) + theme(axis.text.x = element_text(angle = 90, vjust = -0.5)) + ylab('Fraction of SNVs in motifs')
    
  
  pp + theme_Publication() + labs(fill='') + scale_fill_discrete(labels = c("Comotif   ", "Primary motif")) + theme(legend.key.size= unit(0.6, "cm"))
  
  
  f_table_to_word(input_data = compact_tf_table, file_prefix = f_p('correlated_comotif_%s', loc_cell), data_dir = writing_directory)
  ggsave(filename = f_p('%s/snp_distributiion_2_cells.tiff', writing_directory, loc_cell ), dpi = 310, width = 7, height = 7, units = 'in')
  
  
  #print(f_sort_by_col(sub_data2,'Class') %>% filter(Ratio > 0))
  
  
  #Additional number.
  
    return (plot_data)
      
}



#####Analysis##########

#undebug(f_comotif_search2)
#mtrace(f_comotif_search2)
#' write the overall enrichment of TF PWM and comotifs
plot_data = f_comotif_search2(cell_data_list)

#' write the details of each TF
#debug(f_comotif_search)
#undebug(f_homer_motif)
compact_stats = f_comotif_search(cell_data_list)


#undebug(f_plot_comotif_correlation_with_tf_binding)
#mtrace(f_plot_comotif_correlation_with_tf_binding)
mtrace.off()

f_plot_comotif_correlation_with_tf_binding <- function(cell_data_list, compact_stats){
  #Do the dot correlation plot. 
  #check the comotoif_correlation_dot_plot.png in the data dir.
  
  
  plot_data = data.frame()
  
  for (i in 1:nrow(compact_stats)){

    loc_cell = compact_stats[i, 'cell']
    loc_tf   = compact_stats[i, 'tf']
    
    tf_data = cell_data_list[[loc_cell]][[loc_tf]]
    
    if (is.null(tf_data)) next
    
    tf_data = f_add_tf_data_snp_class(loc_tf, tf_data)
    
    table(tf_data$class)
    
    correlated_tfs = unlist(str_split(compact_stats[i,'features'], ', '))
    
    
    
    #change back to original state
    tf_data[tf_data$alt_flip == TRUE, c('ref_tf_dp','alt_tf_dp')] = tf_data[tf_data$alt_flip == TRUE, c('alt_tf_dp', 'ref_tf_dp')]
    
    for (comotif in correlated_tfs){
      
      if(comotif %in% loc_tf) next
      
      
      target_comotif = grep(f_p('pvalue_ratio_%s.*',comotif), colnames(tf_data), value = T)[1]
      binding_comotif = grep(f_p('binding_pwm_%s',comotif), colnames(tf_data), value = T)[1]
      tf_data[tf_data$alt_flip == TRUE, target_comotif] = -1 * tf_data[tf_data$alt_flip == TRUE, target_comotif]
      
      comotif_binding=subset(tf_data, class != 'class1' & tf_data[[binding_comotif]] == TRUE )
      #comotif_binding=subset(tf_data,  tf_data[[binding_comotif]] == TRUE )
      comotif_data=data.frame(
                              tf_ratio = (comotif_binding$ref_tf_dp + 1)/(comotif_binding$alt_tf_dp + 2 + comotif_binding$ref_tf_dp), 
                              #tf_ratio = (comotif_binding$ref_tf_dp + 1)/(comotif_binding$alt_tf_dp + 1), 
        
                              pwm_ratio = comotif_binding[[target_comotif]],
                              ASB = comotif_binding$ASB,
                              tf = f_p('%s-%s(%s)', loc_tf, str_split(comotif,'[.]')[[1]][1], loc_cell))
      
      plot_data = rbind(plot_data, comotif_data)
    }
    
  }
  
  plot_data$tf = f_rename_homer_motifs(plot_data$tf)
  plot_data$tf = f_rename_words_for_figures(plot_data$tf)
  
  #plot_data = subset(plot_data, tf != '')
  
  head(plot_data)
  dim(plot_data)
  plot_data$ASB[plot_data$ASB!= 'ASB'] = 'non-ASB' 
  plot_data$ASB = as.factor(plot_data$ASB)
  print(ggplot(data = plot_data, aes(x = pwm_ratio, y = tf_ratio, color = ASB)) + geom_point(size =1, alpha = 0.6) + facet_wrap(~tf,ncol = 4) +  xlim(-3, 3) +
    #geom_smooth(aes(group=tf), method="lm", fullrange=F, se = F,color = 'black', linetype = 'dashed') 
    xlab("Comotif alteration") + ylab('TF binding imbalance of two alleles'))
  
  ggsave(filename = f_p('%s/comotoif_correlation_dot_plot.png', writing_directory))
  
  
  
  #print(ggplot(data = subset(plot_data, tf == 'smc3-ctcf'), aes(x = pwm_ratio, y = tf_ratio, color = ASB)) + geom_point() +
  #  geom_smooth(aes(group = factor(1)), method="lm", se = F))
  
  #ggsave(filename = f_p('%s/comotoif_correlation_dot_plot-znf143-ctcf.png', writing_directory))
  
  
}


#debug(f_plot_comotif_correlation_with_tf_binding)
f_plot_comotif_correlation_with_tf_binding(cell_data_list, compact_stats)

#str(cell_data_list)
pp<-last_plot()
final_plot<-facetAdjust( pp + geom_point(size =1, alpha = 0.6)+ theme_Publication(base_size = 9) + scale_colour_Publication() + theme(strip.text = element_text(face="plain")) +
              theme(legend.position = c(0.85,0.15),legend.key.size= unit(0.5, "cm"),legend.direction = "vertical") + 
              scale_fill_discrete(labels = c("ASB", "non-ASB")) + labs(color='')
)

#ggsave(filename = f_p('%s/comotoif_correlation_dot_plot.pdf', writing_directory), plot = final_plot , width = 14, height = 14, units = 'in')
ggsave(filename = f_p('%s/comotoif_correlation_dot_plot.tiff', writing_directory), plot = final_plot, width = 7, height = 7, units = 'in')












