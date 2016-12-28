

# Frist approach: Summing the two rounds together
binom.test(head1 + head2, exp1 + exp2)
#p-value = 0.04904


#######Sencond approach based on joint probability of two rounds##########
#' #Generating the probability distribution of each round
dy =  dbinom(x = 0:exp1, exp1,prob = 0.5)
dy2 = dbinom(x = 0:exp2, exp2, prob = 0.5)

#' #Assuming the two rounds are different
observed_prob =dy[head1 + 1] * dy2_pvalue[head2 + 1]

#' #Calculate the joint probability of two rounds
joint_density = as.matrix(data.frame(dy)) %*% as.matrix( t(data.frame(dy2)) )

#Calculate  the events with lower probabiliy than the observed
sum(joint_density[joint_density < observed_prob ])
#0.23
##########This approach is wrong as it doesn't assume the probabilities of two reps are the same ######################



