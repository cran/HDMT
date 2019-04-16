### R code from vignette source 'HDMT.Rnw'

###################################################
### code chunk number 1: loadLibrary
###################################################

library(HDMT)


###################################################
### code chunk number 2: data
###################################################
data(exampleinputdata)
dim(snp_input)
#We only included 10% of the excercise data from the paper 
#due to storage space limit.
dim(exercise_input)
#Each matrix contains two columns of p-values for candidate mediators.
#Column 1 is the p-value of testing if an exposure  is associated with 
#the mediator (alpha!=0). 
#Column 2 is the p-value of testing if a mediator is associated with 
#the outcome adjusted for the exposure (beta!=0) 



###################################################
### code chunk number 3: input_pvalues
###################################################
input_pvalues <- snp_input
#To save time for the illustration, we use 10% of rows; to reproduce the 
# figure in the paper, please don't run the following line
input_pvalues=input_pvalues[sample(1:nrow(input_pvalues),
                  size=ceiling(nrow(input_pvalues)/10)),]


###################################################
### code chunk number 4: nullprop
###################################################
nullprop <- nullestimation(input_pvalues,lambda=0.5)


###################################################
### code chunk number 5: pnull1
###################################################
pnull1<-adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
nullprop$alpha1,nullprop$alpha2,input_pvalues,method=1)


###################################################
### code chunk number 6: fdr
###################################################
fdr <- fdrest(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
nullprop$alpha1,nullprop$alpha2,input_pvalues,method=0)


###################################################
### code chunk number 7: plot_snp
###################################################
pmax <- apply(input_pvalues,1,max)
correct_qqplot(pmax, pnull1)


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
nullprop <- nullestimation(input_pvalues,lambda=0.5)


###################################################
### code chunk number 10: pnull
###################################################
pnull<-adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
nullprop$alpha1,nullprop$alpha2,input_pvalues,method=0)


###################################################
### code chunk number 11: pnull1
###################################################
pnull1<-adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
nullprop$alpha1,nullprop$alpha2,input_pvalues,method=1)


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


