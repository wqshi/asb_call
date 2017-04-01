# asb_call
Compare five approaches to call ASB events, and evaluate them based on the DHS imbalance ratio

#Run five approaches on 37 ASB datasets, and calculate the DHS correlation for each method.

s_cmp_asb_call.R

##Beta-binomial(pool replicates) approach to call ASB events, and test cases.

t_beta_binomial_call.r


##Beta-binomial(consider replicates) approach to call ASB events, and test cases.

s_replicates_betabinomial.R

##Binomial, Binomial with replicates, and edgeR approach to call ASB events, and test cases.

t_asb_calling.r

##Compare the overdispersion parameter in two conditions: (1) TFBS alteration group (2) no TFBS alteration group.

s_overdispersion.R


