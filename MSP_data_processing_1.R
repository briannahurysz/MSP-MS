#enter directory location
directory <- "~/Desktop/R Scripts"
#file "protein-peptides.csv" exported from PEAKS label-free quantification
lfq_filename <- "protein-peptides.csv"
#file "protein-peptides.csv" exported from PEAKS identification
id_filename <- "protein-peptides-FDR.csv"
#specify mass spec samples, replicates should be denoted with the same numbers
sample <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)



##Optional, whether to fix random number generation
#set.seed(12345)
numsample <- length(sample)
library(data.table)
library(qvalue)
library(fitdistrplus)
library(truncnorm)
library(gtools)

#set working directory, all files should be in this directory
setwd(directory)

#function "FDRcorrection", find peptides < 1%FDR
FDRcorrection <- function(lfq, id){
        df <- read.csv(lfq, sep = ",", header=TRUE, row.names = NULL, stringsAsFactors = FALSE)
        refdata <- read.csv(id, sep = ",", header=TRUE, row.names = NULL, stringsAsFactors = FALSE)
        rels <- c(refdata[,4])
        ls <- NULL
        for (i in 1:nrow(df)){
                ls <- c(ls, df[i,],which(rels == df[i,4])[1])
        }
        data <- as.matrix.data.frame(matrix(ls,ncol=(ncol(df)+1),byrow=TRUE))
        cnames <- colnames(df)
        colnames(data) <- c(cnames, "matchwhich")
        data<-data[complete.cases(data[ , ncol(data)]),]
}

#function "matfornorm", filter out Quality <= 0.3 and prepare the file for Normalyzer
matfornorm <- function (sample,df){
        data <- NULL
        ls <- NULL
        cnames <- NULL
        for(i in 1:nrow(df)){
                if (df[i,7] > 0.3) {
                        ls <- c(ls,df[i,1:9],df[i,11:(10+numsample)])
                }
        }
        data <- as.matrix.data.frame(matrix(ls,ncol=(9+numsample),byrow=TRUE))
        cnames <- c(colnames(df)[1:9],colnames(df)[11:(10+numsample)])
        row0 <- c(0,0,0,0,0,0,0,0,0,sample)
        data <- rbind(row0,cnames,gsub("-", 0, matrix(gsub("-", 0, data[1:nrow(data),1:ncol(data)]),ncol = ncol(data))))
        colnames(data) <- cnames
        return (data)
}

data1 <- FDRcorrection(lfq_filename, id_filename)
data2 <- matfornorm(sample,data1)
write.csv(data1, file = "protein-peptides_FDRcorrected.csv", append = FALSE,  sep = ",", row.names = FALSE)
write.table(data2, file = "protein-peptides_for_Normalyzer.txt", append = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
#go to http://normalyzer.immunoprot.lth.se/normalize.php, upload "protein-peptides_for_Normalyzer.txt" for normalization

