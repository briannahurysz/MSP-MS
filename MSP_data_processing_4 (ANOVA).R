#statistic test, anova
#enter directory location
directory <- "~/Desktop/R Scripts"
#specify mass spec samples, replicates should be denoted with the same numbers
sample <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
#Optioanl, whether to fix random numbers
set.seed(123)
numsample <- length(sample)
library(data.table)
library(qvalue)
library(fitdistrplus)
library(truncnorm)
library(gtools)
#library(outliers)
#set working directory, all files should be in this directory
setwd(directory)

samplestr <- table(sample)
timept <- NULL
fccolnames <- NULL
for (i in 1:length(samplestr)){
  timept <- c(timept,as.numeric(samplestr[names(samplestr)==i]))
}

ls <- NULL
anova_pvalue <- NULL
cleavetab <- read.csv("protein-peptides_prepared.csv", sep = ",", header=TRUE, row.names = NULL, stringsAsFactors = FALSE)
  for (i in 1:nrow(cleavetab)){
    ls1 <- c(log2(as.numeric(cleavetab[i,11:(10 + numsample)])),sample)
    df <- data.frame(matrix(ls1,nrow = numsample,byrow=FALSE))
    colnames(df) <- c("intensity", "group")
    res.aov <- aov(as.numeric(df$intensity) ~ as.factor(df$group), data = df)
    pvalue <- as.numeric(unlist(summary(res.aov))[9])
    anova_pvalue <- c(anova_pvalue, pvalue)
  }

data8 <- cbind(cleavetab,anova_pvalue)


#calculate q-value, corrected p-value, permutation (default) or Benjamini-Hochber
#for information of qvalue plots, see https://www.bioconductor.org/packages/devel/bioc/vignettes/qvalue/inst/doc/qvalue.pdf
anova.q <- qvalue(anova_pvalue, fdr.level = 0.05, pfdr = FALSE, lfdr.out = TRUE, lambda = seq(0, 0.95, 0.05), pi0.method = "bootstrap", na.rm = TRUE)
names(anova.q)
summary(anova.q)
anova_qvalue <- anova.q$qvalues
localFDR <- anova.q$lfdr
pi0 <- anova.q$pi0
pdf(file = "qvalue.pdf")
layout(mat=matrix(c(1, 1, 1),nrow=1,ncol=1,byrow=T))
hist(anova_pvalue, nclass = 20)
hist(anova.q)
plot(anova.q)
dev.off()

data9 <- cbind(data8,anova_qvalue,anova.q$significant)



for (i in 1:(length(samplestr)-1)){
  fccolnames <- c(fccolnames, paste("foldchange", i, sep = " "))
}


fc <- NULL
for (i in 1:nrow(data9)){
  counter1 = 10
  for (j in 1:(length(timept)-1)){
    blank <- mean(as.numeric(data9[i, 11:(10 + timept[1])]),na.rm = TRUE)
    counter1 = counter1 + timept[j]
    counter2 = counter1 + timept[j+1]
    fc <- c(fc, mean(as.numeric(data9[i,(counter1+1):counter2]),na.rm = TRUE)/blank)
  }
}
fold_change <- as.matrix.data.frame(matrix(fc,ncol=(length(samplestr)-1),byrow=TRUE))

data10 <- as.matrix(cbind(data9,fold_change))
colnames(data10) <- c(colnames(data9)[1:(ncol(data9)-3)], "Anova p-value", "q-value", "Significant", fccolnames)
write.csv(data10, file = "protein-peptides_results.csv", append = FALSE,  sep = ",", row.names = FALSE)
