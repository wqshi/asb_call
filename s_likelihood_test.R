
install.packages('epicalc')
library(epicalc)
model0 <- glm(case ~ induced + spontaneous, family=binomial, data=infert)
model1 <- glm(case ~ induced, family=binomial, data=infert)
lrtest (model0, model1)