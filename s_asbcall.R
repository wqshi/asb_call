########My Beta-Binomial distribution approach to call ASB events#####
###By adding the mapping bias, this is the new stuff.
##But the data_replicates in missing, need to fix
  source('s_asb_meta_data.R')
  library(dplyr)
  library(VGAM)
  loc_cell = 'gm12878'
  loc_tf = 'ctcf'
   

  data_replicates = f_get_epi_from_paird_database(cell = loc_cell, loc_cell, tf=loc_tf,  f_p("%s/%s/%s/",base_dir, loc_cell, db_version),dp_threshold, 
                                                  fdr, file_suffix = 'database' ,raw_data = TRUE, cell_filter = loc_cell, het_filter=c('het'), 
                                                  labs= lab_list ,target_lab =NA, rep = '.', strong_allele_exchange = F )
  
  data_merge = f_get_epi_from_paird_database(loc_cell, loc_cell, loc_tf, f_p("%s/%s/%s/",base_dir, loc_cell, db_version), dp_threshold = 10, fdr, 
                                           cell_filter = loc_cell, het_filter=c("het"), add_het = TRUE,
                                           db_data = NULL, rep='.', labs= lab_list, target_lab = NA, diffpeak_method = diffpeak_method )

  
  table(data_merge$ASB)
  
  #undebug(f_read_table2)
  
  rep = 1
  lab = 'sydh'
  
  
  grep('ref_', colnames(data_replicates), value = T)
  reps = grep(f_p('%s%s', lab, rep),colnames((data_replicates)), value = TRUE)
  disperse_para = 0.01
  
  data_replicates$mapping_bias = with(data_replicates, ref_simulate_dp_broad1/(ref_simulate_dp_broad1 + alt_simulate_dp_broad1))
  data_replicates$total_dp = rowSums(data_replicates[,reps])
  
  data_replicates_subsets = data_replicates %>% filter(!is.na(mapping_bias), mapping_bias < 1, mapping_bias > 0, total_dp >= 10)
  data_replicates_subsets = f_rename_by_genome_cordinate(data_replicates_subsets)
  
  
  dim(data_replicates_subsets)
  head(data_replicates_subsets)
  
  data_replicates_subsets$A

  nonASB =  subset(data_merge, ASB == 'nonASB')
  rownames(nonASB)
  
  nonASB_replicates = data_replicates_subsets[intersect(rownames(data_replicates_subsets), rownames(nonASB)),]
  
  ref_allele = rowSums(nonASB_replicates[,grep('ref',reps, value = T )])
  total_allele = nonASB_replicates$total_dp
  
  mapping_bias_subset = with(nonASB_replicates, ref_simulate_dp_broad1/(ref_simulate_dp_broad1 + alt_simulate_dp_broad1))
  
  ll <- function(disperse_para){
    dy <- dbetabinom.ab(ref_allele, size = total_allele, shape1 = mapping_bias_subset*(1/disperse_para^2 - 1), shape2 = (1 - mapping_bias_subset)*(1/disperse_para^2 - 1), log = T)
    dy[is.infinite(dy)] = min(dy[!is.infinite(dy)]) - 1
    -sum(dy)
  }
  
  calculate_prop <- function(disperse_para, loc_ref, loc_total, mapping_bias = 0.5){
    dy <- dbetabinom.ab(loc_ref, size = loc_total, shape1 = mapping_bias*(1/disperse_para^2 - 1), shape2 = (1 - mapping_bias)*(1/disperse_para^2 - 1))
    #p_value = apply(cbind(dy, 1-dy), 1, min) * 2
    
    #p_value
    dy
  }
  
  
  x = seq(0.01, 0.99, by = 0.01)
  library(ggplot2)
  lapply(x, ll)
  qplot(x, unlist(lapply(x, ll)))
  
  library(stats4)
  mle_obj = mle(ll, start = list(disperse_para = 0.0001), method = 'Brent' ,lower = c(0.0001), upper = c(0.99))
  #mle(ll2, start = list(disperse_para = 0.0001), method = 'Brent' ,lower = c(0.0001), upper = c(0.99))
  
  
  loc_total = 9
  ref_dp = 7
  loc_total2 = 8 
  ref_dp2 = 6
  mapping_bias = 0.5
  dy =calculate_prop(disperse_para = mle_obj@coef, loc_ref = 0:loc_total, loc_total = loc_total, mapping_bias = 0.5)
  dy_pvalue = 1 - cumsum(dy)
  #replicate2
  dy2 =calculate_prop(disperse_para = mle_obj@coef, loc_ref = 0:loc_total2, loc_total = loc_total2, mapping_bias = 0.5)
  dy2_pvalue = 1 - cumsum(dy2)
  
  cur_pvalue =dy_pvalue[ref_dp + 1] * dy2_pvalue[ref_dp2 + 1]
  
  joint_density = as.matrix(data.frame(dy)) %*% as.matrix( t(data.frame(dy2)) )
  fn = as.vector(joint_density)
  
  joint_pvalue = as.matrix(data.frame(dy_pvalue)) %*% as.matrix( t(data.frame(dy2_pvalue)) )
  
  
  joint_selection = joint_pvalue
  joint_selection[,]=1
  joint_selection[,1:ref_dp]=0
  joint_selection[1:ref_dp2,]=0
  
  sum(joint_density[joint_selection])/sum(joint_density)
  sum(joint_density[joint_pvalue < cur_pvalue])/sum(joint_density)
  binom.test(ref_dp + ref_dp2, loc_total + loc_total2, alternative = 'greater')
  
  
  hist(log10(fn))
  
  
  
  heatmap.2(log10(joint_density),scale = 'none',trace = 'none',Rowv = FALSE, Colv = FALSE)
  
  
  calculate_prop(disperse_para = 0.18, loc_ref = 1, loc_total = 19) * calculate_prop(disperse_para = 0.18, loc_ref = 1, loc_total = 19)
  calculate_prop(disperse_para = 0.18, loc_ref = 2, loc_total = 38) 
  
  
  pvalue = binom.test(10, 30)
  pvalue2 = binom.test(10, 30)
  pvalue$p.value^2


  ###The likelyhood ratio test:




f_cmp_likelihood_and_binomial <- function(ref_dp, ref_dp2, loc_total, loc_total2, input_prob = 0.5){
  
  expected_prop = (ref_dp + ref_dp2) /(loc_total+loc_total2)
  
  null_prob = f_joint_prob(ref_dp, loc_total, ref_dp2, loc_total2, input_prob)
  alt_prob=f_joint_prob(ref_dp, loc_total, ref_dp2, loc_total2, expected_prop)
  
  test_static=-2 * log(null_prob/alt_prob)
  
  flog.info('
            Likeli-ratio test:
            Null prob: %s
            Alt prob: %s
            Alt prob: %s
            test static: %s', null_prob, alt_prob, expected_prop, test_static)
  
  print(1 - pchisq(test_static, df=1, ncp = 0))
  
  print(binom.test(ref_dp + ref_dp2,n = loc_total + loc_total2))
  
}

f_joint_prob <- function(ref_dp, loc_total, ref_dp2, loc_total2, input_prob){
  dbinom(x = ref_dp, size = loc_total, prob = input_prob) * dbinom(x = ref_dp2, size = loc_total2, prob = input_prob)    
}

f_likelihood_call_ASB <- function(ref_dp, ref_dp2, loc_total, loc_total2, null_prob = 0.5){
  expected_prop = (ref_dp + ref_dp2) /(loc_total+loc_total2)
  null_prob = f_joint_prob(ref_dp, loc_total, ref_dp2, loc_total2, null_prob)
  alt_prob=f_joint_prob(ref_dp, loc_total, ref_dp2, loc_total2, expected_prop)
  test_static=-2 * log(null_prob/alt_prob)  
  return (1 - pchisq(test_static, df=1, ncp = 0))
}


f_likelihood_call_ASB(c(1,3), c(2,4), c(3,3), c(4,4), 0.5)


f_cmp_likelihood_and_binomial(4,5,8,10)
f_cmp_likelihood_and_binomial(1,11,3,12)

#The test would be the same if the sum of the ref counts are the same.
#Ref1, ref2, total1, total2
f_cmp_likelihood_and_binomial(2,10,3,12)

#Is this one fair?
f_cmp_likelihood_and_binomial(0,12,3,12)
f_cmp_likelihood_and_binomial(3,9,3,12)




loc_total = 9
ref_dp = 8
loc_total2 = 8 
ref_dp2 = 8

library(ACSWR)
expected_prop = (ref_dp + ref_dp2)/ (loc_total + loc_total2)
MPbinomial(Hp = 0.5, Kp = expected_prop, alpha = 0.05, n = (loc_total + loc_total2 ))
binom.test(x = ref_dp + ref_dp2 , n = loc_total + loc_total2 )


##

#install.packages('lmtest')
library(lmtest)







