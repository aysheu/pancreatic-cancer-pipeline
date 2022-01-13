library(edgeR) #https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

#Below script shows how the pre-processing of individual studies was conducted prior to merging all datasets.

#Set working directory
setwd("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE119794")

#Create a dataframe for the samples to be analysed
samples<-read.table("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE119794/GSE119794_design.csv",sep = ",", header=T) 

dat<-read.table("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE119794/GSE119794, PRJNA490335_Count data_for EdgeR.txt",sep=c("\t"),header=T)
counts<-dat[-1] # Remove the first column of dat which includes the Ensembl gene IDs
rownames(counts)<-dat[,1] # Assign the first column of dat as rownames to counts table

#Create counts object
counts<-DGEList(counts=counts,group=samples$Group) #Specify counts object and grouping variable

#Filter low expression genes
keep <- filterByExpr(counts) #ID low expression genes
counts <- counts[keep, , keep.lib.sizes=FALSE] # Remove low expression genes. The library size for each sample is recalculated to be the sum of the counts left in the rows after removing low expressed genes

#Generate CPM values
cpm_data<-cpm(counts,log=T)
write.table(cpm_data, "GSE119794_cpm_data_filtered_geneIDversions removed.csv",sep=",")