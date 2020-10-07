### R code from vignette source 'HDMT.Rnw'

###################################################
### code chunk number 1: loadLibrary
###################################################

library(HDMT)


###################################################
### code chunk number 2: data
###################################################
data(snp_input)
dim(snp_input)
data(exercise_input)
dim(exercise_input)


###################################################
### code chunk number 3: input_pvalues
###################################################
input_pvalues <- snp_input

input_pvalues=input_pvalues[sample(1:nrow(input_pvalues),
                  size=ceiling(nrow(input_pvalues)/10)),]


###################################################
### code chunk number 4: nullprop
###################################################
nullprop <- null_estimation(input_pvalues,lambda=0.5)
nullprop


###################################################
### code chunk number 5: pnull1
###################################################
pnull1<-adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                        nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=1)


###################################################
### code chunk number 6: plot_snp
###################################################
pmax <- apply(input_pvalues,1,max)
correct_qqplot(pmax, pnull1)


###################################################
### code chunk number 7: fdr
###################################################
efdr <- fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)

plot(pmax[order(pmax)],efdr[order(pmax)],type="l",ylim=c(0,1),
     xlab="p-max",ylab="Estimated FDR")


###################################################
### code chunk number 8: input_pvalues
###################################################
input_pvalues <- exercise_input
#To save time, we use 10% of rows
input_pvalues=input_pvalues[sample(1:nrow(input_pvalues),
                  size=ceiling(nrow(input_pvalues)/10)),]


###################################################
### code chunk number 9: nullprop
###################################################
nullprop <- null_estimation(input_pvalues,lambda=0.5)


###################################################
### code chunk number 10: pnull
###################################################
pnull<-adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)


###################################################
### code chunk number 11: pnull1
###################################################
pnull1<-adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=1)


###################################################
### code chunk number 12: plot_excercise1
###################################################
pmax <- apply(input_pvalues,1,max)
correct_qqplot(pmax, pnull)


###################################################
### code chunk number 13: plot_excercise2
###################################################
correct_qqplot(pmax, pnull1)


###################################################
### code chunk number 14: sessionInfo
###################################################
sessionInfo()


