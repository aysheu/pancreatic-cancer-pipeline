# The object 'data' is a matrix genes x samples
# The object 'gene_list' is a list of genes (eg DEGs) which form the nodes of the hypernetwork
# The object 'rf_variable' is a variable separating the samples used for random forest

setwd("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS")

library(gplots) # heatmap.2 function required
library(dplyr)
library(entropy)

# ROBUSTNESS ANALYSIS

# Random gene sampling

# Data input - source is normalised and scaled (z-scored) transcriptomic data
data<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/GSE93326_matched plus additional stroma_scaledlogTMM_data.csv", sep = ",", header=T)

# Create a data matrix 
matrix_data <-data.matrix(data)
head(matrix_data)

# Create a vector of Ensembl IDs from the row names
rownames<-rownames(data)
head(rownames)
rownames<-as.data.frame(rownames) # dataframe of transcript IDs, protein names etc.

# Read the list of quality genes that passed the low expression threshold 
gene_names<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS/Gene names of all quality genes used in hypernetwork.csv",sep = ",", header = T)

# Identify ECM components that are present in the dataset so these can be used to extract ECM components in the galois
Matrisome<-read.table("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS/Matrisome_geneset.txt", sep = ",", header=T) # list ecm genes
Integrins<-dplyr::filter(gene_names, grepl("ITG",Gene.name))
Keratins<-dplyr::filter(gene_names, grepl("KRT",Gene.name))
ECMs<-gene_names[match(Matrisome$NABA_MATRISOME,gene_names$Gene.name),]
ECMs<-na.omit(ECMs)
ECMs_with_ITG_KRT<-rbind(ECMs,Integrins,Keratins)

# Random sampling for MCTs
iter<-1000 # Run 1000 iterations
N<-7 # Set the N number to 7 as this matches the number of DE-MCTs in central cluster that associated with ECMs
set.seed(1) # Set starting point for randomisation. This ensures all results, figures etc. are the same and reproducible

# Create lists
list_res<-list()
ECMlist<-list()
ECMnames<-list()
Galoislist<-list()

for (i in 1:iter){
  random_gene_rows <- sample(nrow(rownames), size = 7, replace = FALSE)
  genes <- rownames[random_gene_rows,]
  nrow(genes)

  #This one correlates the randomly selected genes with the transcriptome
  corr_data<-cor(t(matrix_data[na.omit(match(genes,rownames(matrix_data))),]), 
                t(matrix_data[-na.omit(match(genes,rownames(matrix_data))),]))
  corr_data<-corr_data[,which(!is.na(colSums(corr_data)))]
                 
  binary<-corr_data
  stdev<-sd(binary)
  binary<-abs(binary)
  binary[which(binary<stdev)]<-0
  binary[which(binary!=0)]<-1
  
  hyp<-binary%*%t(binary)
  
  hm<-heatmap.2(hyp,cexRow=1, cexCol=1, margins=c(10,10), srtCol=45, trace="none")
  
  dend<-as.hclust(hm$rowDendrogram)
  k<-2
  ct<- cutree(dend, k)
  
  res<-data.frame(matrix(NA,nrow=1,ncol=4))
  colnames(res)<-c("Iteration","Sum","Mean","Entropy")
  
  res$Iteration<-i
  res$Sum<-sum(hyp)
  res$Mean<-mean(hyp)
  res$Entropy<-entropy(hyp)/log(N^2) # normalised entropy
  list_res[[i]]<-res
  print(i) 
  
  Cluster1<-names(ct[which(ct!=1)])
  Cluster2<-names(ct[which(ct==1)])
 
  # Use the central cluster gene list to subset the binarized correlation matrix 
  hyp1 <- hyp[match(Cluster1, rownames(hyp)),match(Cluster1, rownames(hyp))]
  hyp2 <-hyp[match(Cluster2, rownames(hyp)),match(Cluster2, rownames(hyp))]
  
  # Calculate which is the more connected one
  mean(hyp1)
  mean(hyp2)
  
  if(mean(hyp1) > mean(hyp2)) {galoisgenes<-Cluster1} else {galoisgenes<-Cluster2}
 
  # Use the central cluster gene list to subset the binarized correlation matrix 
  galois <- binary[match(galoisgenes, rownames(binary)),]
  
  # Identify which of the columns of this matrix now have a 1 in every row (i.e. which of the non-DE genes are correlated with every one of the central cluster genes?)
  galois2<- galois[, which(colSums(galois) > nrow(galois)*0.9)]
  
  # Match colnames of galois2 with whole matrisome list
  nonDE_genes<-colnames(galois2)
  galois_ECM<-ECMs_with_ITG_KRT$Gene.name[match(nonDE_genes,ECMs_with_ITG_KRT$Gene.stable.ID)] 
  galois3<- na.omit(galois_ECM)
   
  # Calculate the size and the name of the list
  ECMlist[[i]]<-NROW(galois3)
  ECMnames[[i]]<-galois3
  Galoislist[[i]]<-galois2
}

# View the results including sum, mean, and entropy in all iterations
res<-as.data.frame(do.call(rbind,list_res))

ECMlist<-as.data.frame(ECMlist)
write.csv(ECMlist, "ECMlist from iteration for MCTs.csv")

saveRDS(ECMnames, file = "ECMnames_final.Rds")
saveRDS(ECMlist, file = "ECMist_final.Rds")
saveRDS(res, file = "ECMres_final.Rds")

# Random sampling for ECMs

setwd("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS")

library(gplots) # heatmap.2 function required
library(entropy)

# Data input - source is normalised and scaled (z-scored) transcriptomic data
data<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/GSE93326_matched plus additional stroma_scaledlogTMM_data.csv", sep = ",", header=T)

# Create a data matrix 
matrix_data <-data.matrix(data)
head(matrix_data)

# Create a vector of Ensembl IDs from the row names
rownames<-rownames(data)
head(rownames)
rownames<-as.data.frame(rownames) # dataframe of transcript IDs, protein names etc.

iter<-1000 # Run 1000 iterations
N<-376 # Set the N number to 376 as this matches the number of DE-ECMs in non-central cluster that associated with MCTs
set.seed(1) # set starting point for randomisation. This ensures all results, figures etc. are the same and reproducible

# Create lists
list_res<-list()
MCTlist<-list()
MCTnames<-list()

# See list of MCTs
MCTs<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/MCTlist.csv", sep=",", header=T)

for (i in 1:iter){
  random_gene_rows <- sample(nrow(rownames), size = 376, replace = FALSE)
  genes <- rownames[random_gene_rows,]
  nrow(genes)
  
  #This one correlates the randomly selected genes with the transcriptome
  corr_data<-cor(t(matrix_data[na.omit(match(genes,rownames(matrix_data))),]), 
                 t(matrix_data[-na.omit(match(genes,rownames(matrix_data))),]))
  corr_data<-corr_data[,which(!is.na(colSums(corr_data)))]
  
    binary<-corr_data
  stdev<-sd(binary)
  binary<-abs(binary)
  binary[which(binary<stdev)]<-0
  binary[which(binary!=0)]<-1
  
  hyp<-binary%*%t(binary)
  
  hm<-heatmap.2(hyp,cexRow=1, cexCol=1, margins=c(10,10), srtCol=45, trace="none")
  
  dend<-as.hclust(hm$rowDendrogram)
  k<-2
  ct<- cutree(dend, k)
  
  res<-data.frame(matrix(NA,nrow=1,ncol=4))
  colnames(res)<-c("Iteration","Sum","Mean","Entropy")
  
  res$Iteration<-i
  res$Sum<-sum(hyp)
  res$Mean<-mean(hyp)
  res$Entropy<-entropy(hyp)/log(N^2) # normalised entropy
  list_res[[i]]<-res
  print(i) 
  
  Cluster1<-names(ct[which(ct!=1)])
  Cluster2<-names(ct[which(ct==1)])

  
  # Use the cluster1 and cluster 2 gene list to subset the hypernetwork matrix 
  hyp1 <- hyp[match(Cluster1, rownames(hyp)),match(Cluster1, rownames(hyp))]
   hyp2 <-hyp[match(Cluster2, rownames(hyp)),match(Cluster2, rownames(hyp))]
  
  # Calculate which is the more connected one
  mean(hyp1)
  mean(hyp2)
  
  if(mean(hyp1) > mean(hyp2)) {galoisgenes<-Cluster1} else {galoisgenes<-Cluster2}
 
  # Use the central cluster gene list to subset the binarized correlation matrix 
  galois <- binary[match(galoisgenes, rownames(binary)),]
  
  # Identify which of the columns of this matrix now have a 1 in 90% of the rows (i.e. which of the non-DE genes are correlated with 90% of the central cluster genes?
  galois2<- galois[, which(colSums(galois) > nrow(galois)*0.9)]
  
  # Match colnames of galois2 with whole matrisome list
  nonDE_genes<-colnames(galois2)
  
  # Match the gene names of the "nonDE_gene" dataset with the genes in the MCT table
    galois_MCT<-MCTs$GeneName[match(nonDE_genes,MCTs$ENSG)]
    galois3<- na.omit(galois_MCT)
 
  # calculate the size of the list
  MCTlist[[i]]<-NROW(galois3)
  MCTnames[[i]]<-galois3
  
  }

MCTlist2<-as.data.frame(MCTlist)
write.csv(MCTlist2, "MCTlist_10000iter.csv")

res<-as.data.frame(do.call(rbind,list_res))

saveRDS(MCTnames, file = "MCTnames_10000iter.Rds")
saveRDS(MCTlist, file = "MCTlist_10000iter.Rds")
saveRDS(res, file = "res_10000iter.Rds")
save(nonDE_genes, res, galois_MCT, galois3, MCTlist, MCTnames, file = "ECMiteration_10000iter.Rdata")


# Direct path analysis

# Data input - source is normalised and scaled (z-scored) transcriptomic data
data<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/GSE93326_matched plus additional stroma_scaledlogTMM_data.csv", sep = ",", header=T)

# Create a vector of Ensembl IDs from the row names
rownames<-rownames(data)
head(rownames)

# Create a data matrix
matrix_data <-data.matrix(data)
head(matrix_data)

## Plot matrix_data
plot(matrix_data)

# Read the differentially expressed genes 
DEGs<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/Genes with 5% FDR-adjusted p-values_matched plus additional stroma.csv", sep = ",", header=T)

# Assign the first column of DEGs as row names
rownames(DEGs)<-DEGs$X

# Create a vector of row names from DEGs and convert it to a data frame
DEGs_rowname<-rownames(DEGs)
DEGs_rowname<-as.data.frame(DEGs_rowname)
DEGs2<-cbind(DEGs_rowname,DEGs)
rownames(DEGs2)<-1:nrow(DEGs2)


# Filter out the ECMs of interest from DEGs2

# First identify integrins and keratins from the list as they are not part of the matrisome geneset
Integrins<-dplyr::filter(DEGs2, grepl("ITG", Gene_name))
Keratins<-dplyr::filter(DEGs2, grepl("KRT", Gene_name))

# List of matrisome and matrisome associated genes from Molecular Signature Database (MSigDB)
# Ref paper:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5013529/
# MSigDB ref for matrisome genes:https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=NABA_MATRISOME

# Read matrisome geneset: Ensemble of genes encoding extracellular matrix and extracellular matrix-associated proteins
# The list includes members mapped to 1026 genes
Matrisome<-read.table("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS/Matrisome_geneset.txt", sep = ",", header=T)

# Change the column name in Matrisome dataset to be able to combine with Integrins and Keratins
names(Matrisome)[1]<- "Gene_name"

# Create a new column called Matrisome.match in DEGs2 to match the differentially-expressed genes with matrisome geneset
DEGs2['Matrisome.match'] <- NA
DEGs2$Matrisome.match <- Matrisome$Gene_name[match(DEGs2$Gene_name, Matrisome$Gene_name)]

# Drop NAs from the list
library(tidyr)
DEGs_matrisome<-DEGs2 %>% drop_na()
DEGs_matrisome<-DEGs_matrisome[,1:7] #to match the column numbers of integrins and keratins dataset to combine

# Combine DE-matrisome gene set with Integrins and Keratins
All_ECM_genes<- rbind(Integrins, Keratins, DEGs_matrisome)

# Filter out the MCTs from DEGs2 
MCTs<-dplyr::filter(DEGs2, grepl("SLC16", Gene_name))
MCTs<-MCTs[,-8]

# Combine ECMs and MCTs
MCT_ECM<-rbind(MCTs,All_ECM_genes)

# Assign the names of DE-matrisome genes to "gene_list" to be used when creating a correlation matrix downstream
gene_list<-MCT_ECM$DEGs_rowname

# Generate a correlation matrix 
# Correlation matrix based on complete dataset (epithelium and stroma). Correlating MCT_ECM against each other
cor_matrix<-cor(t(matrix_data[na.omit(match(gene_list,rownames(matrix_data))),]), # omit NA matches for instances where genes not found in matrix
                t(matrix_data[na.omit(match(gene_list,rownames(matrix_data))),])) 

# Examine the distribution of correlation matrix values for normality with a histogram
hist(cor_matrix, main = paste(""), xlab = "Correlation r-values", col = "darkgrey", cex.lab=1.5, cex.axis=1.5) #ref: https://stackoverflow.com/questions/4241798/how-to-increase-font-size-in-a-plot-in-r and https://www.datamentor.io/r-programming/histogram/
hist(cor_matrix[1:9,10:511])


# Calculate standard deviation
sd(cor_matrix)

# Generate absolute values for all correlation r-values. This matrix will be binarized
binary<-abs(cor_matrix) 

# Set r-value thresold to SD for binarisation
thresh<-sd(cor_matrix) # set r-value threshold as SD for binarisation
binary[which(binary>thresh)]<-1 # set any values greater than the threshold to 1
binary[which(binary!=1)]<-0 # set any values which aren't 1s as 0s 

# Hypernetwork generation and visualisation
hyp<-binary%*%t(binary) # Generate the hypernetwork by multiplying the binary matrix by the transpose of itself. ('%*%' is the operator for matrix multiplication ). This represents the adjacency matrix

# Generate heatmap (source for margin correction: https://stackoverflow.com/questions/21427863/how-to-add-more-margin-to-a-heatmap-2-plot-with-the-png-device)
# Further info on heatmap.2 function: https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2

#hm<-heatmap.2(hyp,trace="none") # generate a heatmap from the hypernetwork matrix, save the heatmap object as it contains dendrograms.Exclude the trace which aims to separates rows and columns
hm<-heatmap.2(hyp, cexRow=1, cexCol=1, margins=c(10,10), srtCol=45, trace="none")

# Customise heatmap axes labels
hm<-heatmap.2(hyp, cexRow=1, cexCol=1, margins=c(2,2), trace="none", xlab = "Differentially expressed MCT and ECM genes", ylab = "Differentially expressed MCT and ECM genes", labRow = FALSE, labCol = FALSE)

# Conduct direct path analysis

# S= (G - I + D((G-I)G))G^-1

# S = direct network (?)
# G = Hypernetwork ===>This equals to the hypernetwork/adjacency matrix
# I = Identity Matrix (diagonal ones) #poss also diagonal of hypernetwork?
# D = diagonal only

G<-cor_matrix # G is the correlation matrix
I<-matrix(0,nrow=nrow(G),ncol=ncol(G));diag(I)<-diag(G) #I is the identity matrix but the diagonals are set to the diagonals of G rather than 1s


GIG<-(G-I)%*%G
DGIG<-matrix(0,nrow=nrow(G),ncol=ncol(G));diag(DGIG)<-diag(GIG)

inv_G<-MASS::ginv(G)

S = (G - I + DGIG)%*%inv_G
colnames(S)<-rownames(S)

write.table(S,"MCT_ECM_Silenced_Hypernetwork_based on cor matrix.txt",sep="\t",col.names=NA)
write.table(G,"MCT_ECM_Comparison_Hypernetwork based on cor matrix.txt",sep="\t",col.names=NA)

hist(S)
hist(S, main = paste(""), xlab = "Correlation r-values in S", col = "darkgrey", cex.lab=1.5, cex.axis=1.5) #ref: https://stackoverflow.com/questions/4241798/how-to-increase-font-size-in-a-plot-in-r and https://www.datamentor.io/r-programming/histogram/

hist(S[1:9,10:511])

# Run statistics
zscore <- (0.100120846621097 - mean(S))/sd(S) # The input value here is the relationship directness score for a pair of genes
pval<-2*pnorm(zscore, mean = 0, sd = 1, lower.tail = FALSE)

zscore_bottom <- (-0.184793232644418 - mean(S))/sd(S) # The input value here is the relationship directness score for a pair of genes
pval_bottom<-2*pnorm(zscore_bottom, mean = 0, sd = 1, lower.tail = T)


# Customise heatmap axes labels
hm_sil<-heatmap.2(S, cexRow=1, cexCol=1, margins=c(2,2), trace="none", xlab = "Differentially expressed MCT and ECM genes", ylab = "Differentially expressed MCT and ECM genes", labRow = FALSE, labCol = FALSE)

# Generate directly comparable heatmaps
hm_hyp<-heatmap.2(G,trace="none")

G_re<-G[rev(hm_hyp$rowInd),hm_hyp$rowInd]
heatmap.2(G_re,trace="none",Rowv = F,Colv = F)

S_re<-S[rev(hm_hyp$rowInd),hm_hyp$rowInd]
heatmap.2(S_re,trace="none",Rowv = F,Colv = F)


# Save a list of the genes in hyp/S matrix to convert ENSG IDs to gene names in Ensembl's online Biomart tool
MCT_ECM_ENSG<-rownames(hyp)
write.table(MCT_ECM_ENSG,"Ensembl IDs of all DE-MCT and ECM genes.csv", quote = FALSE, row.names = FALSE)


##RANDOM FOREST ANALYSIS

#install.packages("ROCR")
#install.packages("pROC")
#install.packages("verification")

# RUN RANDOM FOREST BASED ON CLUSTERED MCTs AND ECMs COMBINED

setwd("F:/Latif Lab back up/AYSHE/PDAC/Birnbaum et al._Validation study")

## Read sample list from the design table
Group<-read.csv("./GSE93326_combined 57 stroma_design.csv", sep = ",", header=T)
Group$SampleID<-gsub("/",".",Group$SampleID)
samples_Birnbaum<-read.table("./Birnbaum_sample design.csv",sep = ",", header=T)

gsub("/",".",Group$SampleID)

# Data input - source is normalised and scaled (z-scored) transcriptomic data 
data_primary<-read.csv("./GSE93326_matched plus additional stroma_scaledlogTMM_data.csv", sep = ",", header=T)
data_birnbaum<-read.table("./Birnbaum_matched epithelium stroma_scaledlogTMM_data.csv", sep=",", header=T)

# Read clustered MCTs from the primary dataset
#MCT_cluster<-read.csv("./CCG_MCTs_stroma_using_r-value_threshold_of_0.2.csv", sep=",", header=T)
#MCT_cluster<-MCT_cluster[,-1]

# Read clustered ECMs from the primary dataset
#ECM_cluster<-read.csv("./NCG_ECM_stroma.csv", sep=",", header=T)
#ECM_cluster<-ECM_cluster[,-1]
#list<-c(MCT_cluster,ECM_cluster) #concatinate both list
#write.csv(list, "MCT_ECM combined.csv")

MCT_ECM_cluster<-read.csv("./MCT_ECM_names_combined.csv", sep=",", header=T)
View(head(data_primary))


# Birnbaum dataset

rf_data_birnbaum<-data_birnbaum[na.omit(match(MCT_ECM_cluster$Gene.name,rownames(data_birnbaum))),]
rf_data_birnbaum_t<-t(rf_data_birnbaum)
rf_data_birnbaum_t<-apply(rf_data_birnbaum_t,2,as.numeric)
rf_data_birnbaum_t<-as.data.frame(rf_data_birnbaum_t)
rownames(rf_data_birnbaum_t)<-samples_Birnbaum$SampleID
samples_Birnbaum$Group<-as.factor(samples_Birnbaum$Group)

colnames(rf_data_birnbaum_t)<-MCT_ECM_cluster$Gene.stable.ID[na.omit(match(colnames(rf_data_birnbaum_t),MCT_ECM_cluster$Gene.name))]

# Create RF datasets

rf_data_primary<-data_primary[na.omit(match(MCT_ECM_cluster$Gene.stable.ID,rownames(data_primary))),]
rf_data_primary_t<-t(rf_data_primary)
rf_data_primary_t<-apply(rf_data_primary_t,2,as.numeric)
rf_data_primary_t<-as.data.frame(rf_data_primary_t)
rownames(rf_data_primary_t)<-Group$SampleID
Group$Group<-as.factor(Group$Group)

# Match gene list

overlap<-intersect(colnames(rf_data_birnbaum_t), colnames(rf_data_primary_t))

rf_data_primary_t<-rf_data_primary_t[,na.omit(match(overlap,colnames(rf_data_primary_t)))]


# Run smote

library(smotefamily)

set.seed(222)
rf_smote_primary<-SMOTE(rf_data_primary_t,Group$Group, K=5, dup_size = 1)

rf_primary<-rf_smote_primary[1][["data"]]
rf_primary$class<-rf_primary$class

#BiocManager::install("randomForest")
library(randomForest)
library(datasets)
library(caret)

# Read the data
data<-rf_primary
str(data)

# Set the Group as a factor
data$class <- as.factor(data$class)

# Data partition

set.seed(222)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1,]
test <- data[ind==2,]

# Run random forest for the training set

rf_primary_train <- randomForest(class~., data=train, proximity=TRUE) 
print(rf_primary_train)


# Prediction & Confusion Matrix - train data

p1_train <- predict(rf_primary_train, train)
confusionMatrix(p1_train, train$ class)

# Prediction & Confusion Matrix - test data

p2_test <- predict(rf_primary_train, test)
confusionMatrix(p2_test, test$ class)

# Error rate of Random Forest

plot(rf_primary_train)

# Tune mtry

t <- tuneRF(train[,-316], train[,316],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 150,
            trace = TRUE,
            improve = 0.05)

# No. of nodes for the trees

hist(treesize(rf_primary_train),
     main = "No. of Nodes for the Trees",
     col = "green")

# Variable Importance
varImpPlot(rf_primary_train,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
importance(rf_primary_train)


# Validation

p_validation<- predict(rf_primary_train, rf_data_birnbaum_t)
confusionMatrix(p_validation,samples_Birnbaum$Group)

# OOB ROC of training model

rf<-rf_primary_train

# The `pROC' package implements various AUC functions.
# Calculate the Area Under the Curve (AUC).
pROC::roc(rf$y, as.numeric(rf$predicted))

# Calculate the AUC Confidence Interval.
pROC::ci.auc(rf$y, as.numeric(rf$predicted))

# List the importance of the variables.
rn <- round(randomForest::importance(rf), 2)
head(rn[order(rn[,1], decreasing=TRUE),])

library(verification)
aucc <- verification::roc.area(as.integer(as.factor(data[ind==1, ncol(data)]))-1,
                               rf$votes[,2])$A
verification::roc.plot(as.integer(as.factor(data[ind==1, ncol(data)]))-1,
                       rf$votes[,2], main="")
legend("bottomright", bty="n",
       sprintf("Area Under the Curve (AUC) = %1.3f", aucc))
title(main="OOB ROC Curve Random Forest A0",
      sub=paste("Rattle", format(Sys.time(), "%Y-%b-%d %H:%M:%S"), Sys.info()["user"]))


# ROC of validation model

newdata<-rf_data_birnbaum_t
newdata$class<-samples_Birnbaum$Group

library(ROCR)

# ROC Curve: requires the ggplot2 package.
library(ggplot2, quietly=TRUE)

# Generate an ROC Curve for the rf model on A0 [validate].

pr <- predict(rf, newdata=newdata,
              type    = "prob")[,2]

# Remove observations with missing target.

no.miss   <- na.omit(newdata$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr[-miss.list], no.miss)
} else
{
  pred <- prediction(pr, no.miss)
}

pe <- performance(pred, "tpr", "fpr")
au <- performance(pred, "auc")@y.values[[1]]
pd <- data.frame(fpr=unlist(pe@x.values), tpr=unlist(pe@y.values))
p <- ggplot(pd, aes(x=fpr, y=tpr))
p <- p + geom_line(colour="red")
p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")
p <- p + ggtitle("ROC Curve Random Forest Validation Tissue")
p <- p + theme(plot.title=element_text(size=10))
p <- p + geom_line(data=data.frame(), aes(x=c(0,1), y=c(0,1)), colour="grey")
p <- p + annotate("text", x=0.50, y=0.00, hjust=0, vjust=0, size=5,
                  label=paste("AUC =", round(au, 2)))
print(p)

# Calculate the area under the curve for the plot.


# Remove observations with missing target.

no.miss   <- na.omit(newdata$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr[-miss.list], no.miss)
} else
{
  pred <- prediction(pr, no.miss)
}
performance(pred, "auc")



# ROC of test model

test_data<-test

library(ROCR)

# ROC Curve: requires the ggplot2 package.
library(ggplot2, quietly=TRUE)

# Generate an ROC Curve for the rf model on A0 [validate].

pr_test <- predict(rf, newdata=test_data,
                   type    = "prob")[,2]

# Remove observations with missing target.

no.miss   <- na.omit(test_data$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr_test[-miss.list], no.miss)
} else
{
  pred <- prediction(pr_test, no.miss)
}

pe <- performance(pred, "tpr", "fpr")
au <- performance(pred, "auc")@y.values[[1]]
pd <- data.frame(fpr=unlist(pe@x.values), tpr=unlist(pe@y.values))
p <- ggplot(pd, aes(x=fpr, y=tpr))
p <- p + geom_line(colour="red")
p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")
p <- p + ggtitle("ROC Curve Random Forest Test Set")
p <- p + theme(plot.title=element_text(size=10))
p <- p + geom_line(data=data.frame(), aes(x=c(0,1), y=c(0,1)), colour="grey")
p <- p + annotate("text", x=0.50, y=0.00, hjust=0, vjust=0, size=5,
                  label=paste("AUC =", round(au, 2)))
print(p)



## RUN RANDOM FOREST BASED ON CLUSTERED MCTs ALONE

# Subset MCTs from the MCT_ECM_cluster
MCTs<-dplyr::filter(MCT_ECM_cluster, grepl("SLC16", Gene.name))

# Birnbaum dataset

rf_data_birnbaum<-data_birnbaum[na.omit(match(MCTs$Gene.name,rownames(data_birnbaum))),]
rf_data_birnbaum_t<-t(rf_data_birnbaum)
rf_data_birnbaum_t<-apply(rf_data_birnbaum_t,2,as.numeric)
rf_data_birnbaum_t<-as.data.frame(rf_data_birnbaum_t)
rownames(rf_data_birnbaum_t)<-samples_Birnbaum$SampleID
samples_Birnbaum$Group<-as.factor(samples_Birnbaum$Group)

colnames(rf_data_birnbaum_t)<-MCTs$Gene.stable.ID[na.omit(match(colnames(rf_data_birnbaum_t),MCTs$Gene.name))]

# Create RF datasets

rf_data_primary<-data_primary[na.omit(match(MCTs$Gene.stable.ID,rownames(data_primary))),]
rf_data_primary_t<-t(rf_data_primary)
rf_data_primary_t<-apply(rf_data_primary_t,2,as.numeric)
rf_data_primary_t<-as.data.frame(rf_data_primary_t)
rownames(rf_data_primary_t)<-Group$SampleID
Group$Group<-as.factor(Group$Group)

# Match gene list

overlap<-intersect(colnames(rf_data_birnbaum_t), colnames(rf_data_primary_t))

rf_data_primary_t<-rf_data_primary_t[,na.omit(match(overlap,colnames(rf_data_primary_t)))]


# Run smote

library(smotefamily)

set.seed(222)
rf_smote_primary<-SMOTE(rf_data_primary_t,Group$Group, K=5, dup_size = 1)

rf_primary<-rf_smote_primary[1][["data"]]
rf_primary$class<-rf_primary$class

#BiocManager::install("randomForest")
library(randomForest)
library(datasets)
library(caret)

# Read the data
data<-rf_primary
str(data)

# Set the Group as a factor
data$class <- as.factor(data$class)

# Data partition

set.seed(222)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1,]
test <- data[ind==2,]

# Run random forest for the training set

rf_primary_train <- randomForest(class~., data=train, proximity=TRUE) 
print(rf_primary_train)


# Prediction & Confusion Matrix - train data

p1_train <- predict(rf_primary_train, train)
confusionMatrix(p1_train, train$ class)

# Prediction & Confusion Matrix - test data

p2_test <- predict(rf_primary_train, test)
confusionMatrix(p2_test, test$ class)

# Error rate of Random Forest

plot(rf_primary_train)

# Tune mtry

t <- tuneRF(train[,-7], train[,7],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 150,
            trace = TRUE,
            improve = 0.05)

# No. of nodes for the trees

hist(treesize(rf_primary_train),
     main = "No. of Nodes for the Trees",
     col = "green")

# Variable Importance
varImpPlot(rf_primary_train,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
importance(rf_primary_train)


# Validation

p_validation<- predict(rf_primary_train, rf_data_birnbaum_t)
confusionMatrix(p_validation,samples_Birnbaum$Group)


# OOB ROC of training model

rf<-rf_primary_train

# The `pROC' package implements various AUC functions.
# Calculate the Area Under the Curve (AUC).
pROC::roc(rf$y, as.numeric(rf$predicted))

# Calculate the AUC Confidence Interval.
pROC::ci.auc(rf$y, as.numeric(rf$predicted))

# List the importance of the variables.
rn <- round(randomForest::importance(rf), 2)
head(rn[order(rn[,1], decreasing=TRUE),])

library(verification)
aucc <- verification::roc.area(as.integer(as.factor(data[ind==1, ncol(data)]))-1,
                               rf$votes[,2])$A
verification::roc.plot(as.integer(as.factor(data[ind==1, ncol(data)]))-1,
                       rf$votes[,2], main="")
legend("bottomright", bty="n",
       sprintf("Area Under the Curve (AUC) = %1.3f", aucc))
title(main="OOB ROC Curve Random Forest A0",
      sub=paste("Rattle", format(Sys.time(), "%Y-%b-%d %H:%M:%S"), Sys.info()["user"]))


# ROC of validation model

newdata<-rf_data_birnbaum_t
newdata$class<-samples_Birnbaum$Group

library(ROCR)

# ROC Curve: requires the ggplot2 package.
library(ggplot2, quietly=TRUE)

# Generate an ROC Curve for the rf model on A0 [validate].

pr <- predict(rf, newdata=newdata,
              type    = "prob")[,2]

# Remove observations with missing target.

no.miss   <- na.omit(newdata$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr[-miss.list], no.miss)
} else
{
  pred <- prediction(pr, no.miss)
}

pe <- performance(pred, "tpr", "fpr")
au <- performance(pred, "auc")@y.values[[1]]
pd <- data.frame(fpr=unlist(pe@x.values), tpr=unlist(pe@y.values))
p <- ggplot(pd, aes(x=fpr, y=tpr))
p <- p + geom_line(colour="red")
p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")
p <- p + ggtitle("ROC Curve Random Forest Validation Tissue")
p <- p + theme(plot.title=element_text(size=10))
p <- p + geom_line(data=data.frame(), aes(x=c(0,1), y=c(0,1)), colour="grey")
p <- p + annotate("text", x=0.50, y=0.00, hjust=0, vjust=0, size=5,
                  label=paste("AUC =", round(au, 2)))
print(p)

# Calculate the area under the curve for the plot.


# Remove observations with missing target.

no.miss   <- na.omit(newdata$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr[-miss.list], no.miss)
} else
{
  pred <- prediction(pr, no.miss)
}
performance(pred, "auc")


# ROC of test model

test_data<-test

library(ROCR)

# ROC Curve: requires the ggplot2 package.
library(ggplot2, quietly=TRUE)

# Generate an ROC Curve for the rf model on A0 [validate].

pr_test <- predict(rf, newdata=test_data,
                   type    = "prob")[,2]

# Remove observations with missing target.

no.miss   <- na.omit(test_data$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr_test[-miss.list], no.miss)
} else
{
  pred <- prediction(pr_test, no.miss)
}

pe <- performance(pred, "tpr", "fpr")
au <- performance(pred, "auc")@y.values[[1]]
pd <- data.frame(fpr=unlist(pe@x.values), tpr=unlist(pe@y.values))
p <- ggplot(pd, aes(x=fpr, y=tpr))
p <- p + geom_line(colour="red")
p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")
p <- p + ggtitle("ROC Curve Random Forest Test Set")
p <- p + theme(plot.title=element_text(size=10))
p <- p + geom_line(data=data.frame(), aes(x=c(0,1), y=c(0,1)), colour="grey")
p <- p + annotate("text", x=0.50, y=0.00, hjust=0, vjust=0, size=5,
                  label=paste("AUC =", round(au, 2)))
print(p)



## RUN RANDOM FOREST BASED ON CLUSTERED ECMs ALONE

# Subset ECMs from the MCT_ECM_cluster
ECMs<-MCT_ECM_cluster[-c(85,95,119,210,298,324,376),]

# Birnbaum dataset

rf_data_birnbaum<-data_birnbaum[na.omit(match(ECMs$Gene.name,rownames(data_birnbaum))),]
rf_data_birnbaum_t<-t(rf_data_birnbaum)
rf_data_birnbaum_t<-apply(rf_data_birnbaum_t,2,as.numeric)
rf_data_birnbaum_t<-as.data.frame(rf_data_birnbaum_t)
rownames(rf_data_birnbaum_t)<-samples_Birnbaum$SampleID
samples_Birnbaum$Group<-as.factor(samples_Birnbaum$Group)

colnames(rf_data_birnbaum_t)<-ECMs$Gene.stable.ID[na.omit(match(colnames(rf_data_birnbaum_t),ECMs$Gene.name))]

# Create RF datasets

rf_data_primary<-data_primary[na.omit(match(ECMs$Gene.stable.ID,rownames(data_primary))),]
rf_data_primary_t<-t(rf_data_primary)
rf_data_primary_t<-apply(rf_data_primary_t,2,as.numeric)
rf_data_primary_t<-as.data.frame(rf_data_primary_t)
rownames(rf_data_primary_t)<-Group$SampleID
Group$Group<-as.factor(Group$Group)

# Match gene list

overlap<-intersect(colnames(rf_data_birnbaum_t), colnames(rf_data_primary_t))

rf_data_primary_t<-rf_data_primary_t[,na.omit(match(overlap,colnames(rf_data_primary_t)))]


# Run smote

library(smotefamily)

set.seed(222)
rf_smote_primary<-SMOTE(rf_data_primary_t,Group$Group, K=5, dup_size = 1)

rf_primary<-rf_smote_primary[1][["data"]]
rf_primary$class<-rf_primary$class

#BiocManager::install("randomForest")
library(randomForest)
library(datasets)
library(caret)

# Read the data
data<-rf_primary
str(data)

# Set the Group as a factor
data$class <- as.factor(data$class)

# Data partition

set.seed(222)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1,]
test <- data[ind==2,]

# Run random forest for the training set

rf_primary_train <- randomForest(class~., data=train, proximity=TRUE) 
print(rf_primary_train)


# Prediction & Confusion Matrix - train data

p1_train <- predict(rf_primary_train, train)
confusionMatrix(p1_train, train$ class)

# Prediction & Confusion Matrix - test data

p2_test <- predict(rf_primary_train, test)
confusionMatrix(p2_test, test$ class)

# Error rate of Random Forest

plot(rf_primary_train)

# Tune mtry

t <- tuneRF(train[,-310], train[,310],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 150,
            trace = TRUE,
            improve = 0.05)

# No. of nodes for the trees

hist(treesize(rf_primary_train),
     main = "No. of Nodes for the Trees",
     col = "green")

# Variable Importance
varImpPlot(rf_primary_train,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
importance(rf_primary_train)


# Validation

p_validation<- predict(rf_primary_train, rf_data_birnbaum_t)
confusionMatrix(p_validation,samples_Birnbaum$Group)


# OOB ROC of training model

rf<-rf_primary_train

# The `pROC' package implements various AUC functions.
# Calculate the Area Under the Curve (AUC).
pROC::roc(rf$y, as.numeric(rf$predicted))

# Calculate the AUC Confidence Interval.
pROC::ci.auc(rf$y, as.numeric(rf$predicted))

# List the importance of the variables.
rn <- round(randomForest::importance(rf), 2)
head(rn[order(rn[,1], decreasing=TRUE),])

library(verification)
aucc <- verification::roc.area(as.integer(as.factor(data[ind==1, ncol(data)]))-1,
                               rf$votes[,2])$A
verification::roc.plot(as.integer(as.factor(data[ind==1, ncol(data)]))-1,
                       rf$votes[,2], main="")
legend("bottomright", bty="n",
       sprintf("Area Under the Curve (AUC) = %1.3f", aucc))
title(main="OOB ROC Curve Random Forest A0",
      sub=paste("Rattle", format(Sys.time(), "%Y-%b-%d %H:%M:%S"), Sys.info()["user"]))


# ROC of validation model

newdata<-rf_data_birnbaum_t
newdata$class<-samples_Birnbaum$Group

library(ROCR)

# ROC Curve: requires the ggplot2 package.
library(ggplot2, quietly=TRUE)

# Generate an ROC Curve for the rf model on A0 [validate].

pr <- predict(rf, newdata=newdata,
              type    = "prob")[,2]

# Remove observations with missing target.

no.miss   <- na.omit(newdata$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr[-miss.list], no.miss)
} else
{
  pred <- prediction(pr, no.miss)
}

pe <- performance(pred, "tpr", "fpr")
au <- performance(pred, "auc")@y.values[[1]]
pd <- data.frame(fpr=unlist(pe@x.values), tpr=unlist(pe@y.values))
p <- ggplot(pd, aes(x=fpr, y=tpr))
p <- p + geom_line(colour="red")
p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")
p <- p + ggtitle("ROC Curve Random Forest Validation Tissue")
p <- p + theme(plot.title=element_text(size=10))
p <- p + geom_line(data=data.frame(), aes(x=c(0,1), y=c(0,1)), colour="grey")
p <- p + annotate("text", x=0.50, y=0.00, hjust=0, vjust=0, size=5,
                  label=paste("AUC =", round(au, 2)))
print(p)

# Calculate the area under the curve for the plot.


# Remove observations with missing target.

no.miss   <- na.omit(newdata$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr[-miss.list], no.miss)
} else
{
  pred <- prediction(pr, no.miss)
}
performance(pred, "auc")


# ROC of test model

test_data<-test

library(ROCR)

# ROC Curve: requires the ggplot2 package.
library(ggplot2, quietly=TRUE)

# Generate an ROC Curve for the rf model on A0 [validate].

pr_test <- predict(rf, newdata=test_data,
                   type    = "prob")[,2]

# Remove observations with missing target.

no.miss   <- na.omit(test_data$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr_test[-miss.list], no.miss)
} else
{
  pred <- prediction(pr_test, no.miss)
}

pe <- performance(pred, "tpr", "fpr")
au <- performance(pred, "auc")@y.values[[1]]
pd <- data.frame(fpr=unlist(pe@x.values), tpr=unlist(pe@y.values))
p <- ggplot(pd, aes(x=fpr, y=tpr))
p <- p + geom_line(colour="red")
p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")
p <- p + ggtitle("ROC Curve Random Forest Test Set")
p <- p + theme(plot.title=element_text(size=10))
p <- p + geom_line(data=data.frame(), aes(x=c(0,1), y=c(0,1)), colour="grey")
p <- p + annotate("text", x=0.50, y=0.00, hjust=0, vjust=0, size=5,
                  label=paste("AUC =", round(au, 2)))
print(p)




## RUN RANDOM FOREST BASED ON MCTs NOT IN THE CENTRAL CLUSTER 

# Read list of MCTs not in the central cluster of the hypernetwork
MCTs_NC<-read.csv("./NCG_MCTs_stroma_using_r-value_threshold_of_0.2.csv", sep = ",", header=T)

MCT_NCC<-read.csv("./NCG_MCTs_gene names.csv", sep=",", header=T)


# Birnbaum dataset

rf_data_birnbaum<-data_birnbaum[na.omit(match(MCT_NCC$Gene.name,rownames(data_birnbaum))),]
rf_data_birnbaum_t<-t(rf_data_birnbaum)
rf_data_birnbaum_t<-apply(rf_data_birnbaum_t,2,as.numeric)
rf_data_birnbaum_t<-as.data.frame(rf_data_birnbaum_t)
rownames(rf_data_birnbaum_t)<-samples_Birnbaum$SampleID
samples_Birnbaum$Group<-as.factor(samples_Birnbaum$Group)

colnames(rf_data_birnbaum_t)<-MCT_NCC$Gene.stable.ID[na.omit(match(colnames(rf_data_birnbaum_t),MCT_NCC$Gene.name))]

# Create RF datasets

rf_data_primary<-data_primary[na.omit(match(MCT_NCC$Gene.stable.ID,rownames(data_primary))),]
rf_data_primary_t<-t(rf_data_primary)
rf_data_primary_t<-apply(rf_data_primary_t,2,as.numeric)
rf_data_primary_t<-as.data.frame(rf_data_primary_t)
rownames(rf_data_primary_t)<-Group$SampleID
Group$Group<-as.factor(Group$Group)

# Match gene list

overlap<-intersect(colnames(rf_data_birnbaum_t), colnames(rf_data_primary_t))

rf_data_primary_t<-rf_data_primary_t[,na.omit(match(overlap,colnames(rf_data_primary_t)))]


# Run smote

library(smotefamily)

set.seed(222)
rf_smote_primary<-SMOTE(rf_data_primary_t,Group$Group, K=5, dup_size = 1)

rf_primary<-rf_smote_primary[1][["data"]]
rf_primary$class<-rf_primary$class

#BiocManager::install("randomForest")
library(randomForest)
library(datasets)
library(caret)

# Read the data
data<-rf_primary
str(data)

# Set the Group as a factor
data$class <- as.factor(data$class)

# Data partition

set.seed(222)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1,]
test <- data[ind==2,]

# Run random forest for the training set

rf_primary_train <- randomForest(class~., data=train, proximity=TRUE) 
print(rf_primary_train)


# Prediction & Confusion Matrix - train data

p1_train <- predict(rf_primary_train, train)
confusionMatrix(p1_train, train$ class)

# Prediction & Confusion Matrix - test data

p2_test <- predict(rf_primary_train, test)
confusionMatrix(p2_test, test$ class)

# Error rate of Random Forest

plot(rf_primary_train)

# Tune mtry

t <- tuneRF(train[,-7], train[,7],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 150,
            trace = TRUE,
            improve = 0.05)

# No. of nodes for the trees

hist(treesize(rf_primary_train),
     main = "No. of Nodes for the Trees",
     col = "green")

# Variable Importance
varImpPlot(rf_primary_train,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
importance(rf_primary_train)


# Validation

p_validation<- predict(rf_primary_train, rf_data_birnbaum_t)
confusionMatrix(p_validation,samples_Birnbaum$Group)


## RUN RANDOM FOREST BASED ON ECMs and MCTs NOT IN THE CENTRAL CLUSTER 

# Read list of ECMs and MCTs not in the central cluster of the hypernetwork
ECM_NCC<-read.csv("./CCG_ECMs_NOT in central cluster_gene names.csv", sep=",", header=T)
MCT_NCC<-read.csv("./NCG_MCTs_gene names.csv", sep=",", header=T)
MCT_ECM_NCC<-rbind(MCT_NCC,ECM_NCC)

# Birnbaum dataset

rf_data_birnbaum<-data_birnbaum[na.omit(match(MCT_ECM_NCC$Gene.name,rownames(data_birnbaum))),]
rf_data_birnbaum_t<-t(rf_data_birnbaum)
rf_data_birnbaum_t<-apply(rf_data_birnbaum_t,2,as.numeric)
rf_data_birnbaum_t<-as.data.frame(rf_data_birnbaum_t)
rownames(rf_data_birnbaum_t)<-samples_Birnbaum$SampleID
samples_Birnbaum$Group<-as.factor(samples_Birnbaum$Group)

colnames(rf_data_birnbaum_t)<-MCT_ECM_NCC$Gene.stable.ID[na.omit(match(colnames(rf_data_birnbaum_t),MCT_ECM_NCC$Gene.name))]

# Create RF datasets

rf_data_primary<-data_primary[na.omit(match(MCT_ECM_NCC$Gene.stable.ID,rownames(data_primary))),]
rf_data_primary_t<-t(rf_data_primary)
rf_data_primary_t<-apply(rf_data_primary_t,2,as.numeric)
rf_data_primary_t<-as.data.frame(rf_data_primary_t)
rownames(rf_data_primary_t)<-Group$SampleID
Group$Group<-as.factor(Group$Group)

# Match gene list

overlap<-intersect(colnames(rf_data_birnbaum_t), colnames(rf_data_primary_t))

rf_data_primary_t<-rf_data_primary_t[,na.omit(match(overlap,colnames(rf_data_primary_t)))]


# Run smote

library(smotefamily)

set.seed(222)
rf_smote_primary<-SMOTE(rf_data_primary_t,Group$Group, K=5, dup_size = 1)

rf_primary<-rf_smote_primary[1][["data"]]
rf_primary$class<-rf_primary$class

#BiocManager::install("randomForest")
library(randomForest)
library(datasets)
library(caret)

# Read the data
data<-rf_primary
str(data)

# Set the Group as a factor
data$class <- as.factor(data$class)

# Data partition

set.seed(222)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1,]
test <- data[ind==2,]

# Run random forest for the training set

rf_primary_train <- randomForest(class~., data=train, proximity=TRUE) 
print(rf_primary_train)


# Prediction & Confusion Matrix - train data

p1_train <- predict(rf_primary_train, train)
confusionMatrix(p1_train, train$ class)

# Prediction & Confusion Matrix - test data

p2_test <- predict(rf_primary_train, test)
confusionMatrix(p2_test, test$ class)

# Error rate of Random Forest

plot(rf_primary_train)

# Tune mtry

t <- tuneRF(train[,-96], train[,96],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 150,
            trace = TRUE,
            improve = 0.05)

# No. of nodes for the trees

hist(treesize(rf_primary_train),
     main = "No. of Nodes for the Trees",
     col = "green")

# Variable Importance
varImpPlot(rf_primary_train,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
importance(rf_primary_train)


# Validation

p_validation<- predict(rf_primary_train, rf_data_birnbaum_t)
confusionMatrix(p_validation,samples_Birnbaum$Group)


# OOB ROC of training model

rf<-rf_primary_train

# The `pROC' package implements various AUC functions.
# Calculate the Area Under the Curve (AUC).
pROC::roc(rf$y, as.numeric(rf$predicted))

# Calculate the AUC Confidence Interval.
pROC::ci.auc(rf$y, as.numeric(rf$predicted))

# List the importance of the variables.
rn <- round(randomForest::importance(rf), 2)
head(rn[order(rn[,1], decreasing=TRUE),])

library(verification)
aucc <- verification::roc.area(as.integer(as.factor(data[ind==1, ncol(data)]))-1,
                               rf$votes[,2])$A
verification::roc.plot(as.integer(as.factor(data[ind==1, ncol(data)]))-1,
                       rf$votes[,2], main="")
legend("bottomright", bty="n",
       sprintf("Area Under the Curve (AUC) = %1.3f", aucc))
title(main="OOB ROC Curve Random Forest A0",
      sub=paste("Rattle", format(Sys.time(), "%Y-%b-%d %H:%M:%S"), Sys.info()["user"]))


# ROC of validation model

newdata<-rf_data_birnbaum_t
newdata$class<-samples_Birnbaum$Group

library(ROCR)

# ROC Curve: requires the ggplot2 package.
library(ggplot2, quietly=TRUE)

# Generate an ROC Curve for the rf model on A0 [validate].

pr <- predict(rf, newdata=newdata,
              type    = "prob")[,2]

# Remove observations with missing target.

no.miss   <- na.omit(newdata$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr[-miss.list], no.miss)
} else
{
  pred <- prediction(pr, no.miss)
}

pe <- performance(pred, "tpr", "fpr")
au <- performance(pred, "auc")@y.values[[1]]
pd <- data.frame(fpr=unlist(pe@x.values), tpr=unlist(pe@y.values))
p <- ggplot(pd, aes(x=fpr, y=tpr))
p <- p + geom_line(colour="red")
p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")
p <- p + ggtitle("ROC Curve Random Forest Validation Tissue")
p <- p + theme(plot.title=element_text(size=10))
p <- p + geom_line(data=data.frame(), aes(x=c(0,1), y=c(0,1)), colour="grey")
p <- p + annotate("text", x=0.50, y=0.00, hjust=0, vjust=0, size=5,
                  label=paste("AUC =", round(au, 2)))
print(p)

# Calculate the area under the curve for the plot.


# Remove observations with missing target.

no.miss   <- na.omit(newdata$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr[-miss.list], no.miss)
} else
{
  pred <- prediction(pr, no.miss)
}
performance(pred, "auc")


# ROC of test model

test_data<-test

library(ROCR)

# ROC Curve: requires the ggplot2 package.
library(ggplot2, quietly=TRUE)

# Generate an ROC Curve for the rf model on A0 [validate].

pr_test <- predict(rf, newdata=test_data,
                   type    = "prob")[,2]

# Remove observations with missing target.

no.miss   <- na.omit(test_data$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr_test[-miss.list], no.miss)
} else
{
  pred <- prediction(pr_test, no.miss)
}

pe <- performance(pred, "tpr", "fpr")
au <- performance(pred, "auc")@y.values[[1]]
pd <- data.frame(fpr=unlist(pe@x.values), tpr=unlist(pe@y.values))
p <- ggplot(pd, aes(x=fpr, y=tpr))
p <- p + geom_line(colour="red")
p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")
p <- p + ggtitle("ROC Curve Random Forest Test Set")
p <- p + theme(plot.title=element_text(size=10))
p <- p + geom_line(data=data.frame(), aes(x=c(0,1), y=c(0,1)), colour="grey")
p <- p + annotate("text", x=0.50, y=0.00, hjust=0, vjust=0, size=5,
                  label=paste("AUC =", round(au, 2)))
print(p)




## RUN RANDOM FOREST BASED ON NON-DEGs IN THE PRIMARY DATASET

# Read list of all genes that passed low expression threshold in primary dataset
All_genes_primary<-read.csv("./Unadjusted and adjusted p-values for ALL genes_matched plus additional stroma.csv", sep=",", header=T)

# Read list of DEGs with FDR<0.05 in primary dataset
DEGs_primary<-read.csv("./Genes with 5% FDR-adjusted p-values_matched plus additional stroma.csv", sep=",", header=T)
rownames(DEGs_primary)<-DEGs_primary$X
DEGs_primary<-DEGs_primary[,c(2:6)]

# Subset non-DEGs from all genes
nonDEGs<-subset(All_genes_primary, FDR > 0.05)
nonDEGs$Gene.stable.ID<-rownames(nonDEGs)


# Birnbaum dataset

rf_data_birnbaum<-data_birnbaum[na.omit(match(nonDEGs$Gene_name,rownames(data_birnbaum))),]
rf_data_birnbaum_t<-t(rf_data_birnbaum)
rf_data_birnbaum_t<-apply(rf_data_birnbaum_t,2,as.numeric)
rf_data_birnbaum_t<-as.data.frame(rf_data_birnbaum_t)
rownames(rf_data_birnbaum_t)<-samples_Birnbaum$SampleID
samples_Birnbaum$Group<-as.factor(samples_Birnbaum$Group)

colnames(rf_data_birnbaum_t)<-nonDEGs$Gene.stable.ID[na.omit(match(colnames(rf_data_birnbaum_t),nonDEGs$Gene_name))]


# Create RF datasets

rf_data_primary<-data_primary[na.omit(match(nonDEGs$Gene.stable.ID,rownames(data_primary))),]
rf_data_primary_t<-t(rf_data_primary)
rf_data_primary_t<-apply(rf_data_primary_t,2,as.numeric)
rf_data_primary_t<-as.data.frame(rf_data_primary_t)
rownames(rf_data_primary_t)<-Group$SampleID
Group$Group<-as.factor(Group$Group)

# Match gene list

overlap<-intersect(colnames(rf_data_birnbaum_t), colnames(rf_data_primary_t))

rf_data_primary_t<-rf_data_primary_t[,na.omit(match(overlap,colnames(rf_data_primary_t)))]

# Randomly select 315 nonDEGs
set.seed(222)
rf_data_primary_t<-rf_data_primary_t[,sample(ncol(rf_data_primary_t), 315)]

# Run smote

library(smotefamily)

set.seed(222)
rf_smote_primary<-SMOTE(rf_data_primary_t,Group$Group, K=5, dup_size = 1)

rf_primary<-rf_smote_primary[1][["data"]]
rf_primary$class<-rf_primary$class

#BiocManager::install("randomForest")
library(randomForest)
library(datasets)
library(caret)

# Read the data
data<-rf_primary
str(data)

# Set the Group as a factor
data$class <- as.factor(data$class)

# Data partition

set.seed(222)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1,]
test <- data[ind==2,]

# Run random forest for the training set

rf_primary_train <- randomForest(class~., data=train, proximity=TRUE) 
print(rf_primary_train)


# Prediction & Confusion Matrix - train data

p1_train <- predict(rf_primary_train, train)
confusionMatrix(p1_train, train$ class)

# Prediction & Confusion Matrix - test data

p2_test <- predict(rf_primary_train, test)
confusionMatrix(p2_test, test$ class)

# Error rate of Random Forest

plot(rf_primary_train)

# Tune mtry

t <- tuneRF(train[,-301], train[,301],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 150,
            trace = TRUE,
            improve = 0.05)

# No. of nodes for the trees

hist(treesize(rf_primary_train),
     main = "No. of Nodes for the Trees",
     col = "green")

# Variable Importance
varImpPlot(rf_primary_train,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
importance(rf_primary_train)


# Validation

p_validation<- predict(rf_primary_train, rf_data_birnbaum_t)
confusionMatrix(p_validation,samples_Birnbaum$Group)


# OOB ROC of training model

rf<-rf_primary_train

# The `pROC' package implements various AUC functions.
# Calculate the Area Under the Curve (AUC).
pROC::roc(rf$y, as.numeric(rf$predicted))

# Calculate the AUC Confidence Interval.
pROC::ci.auc(rf$y, as.numeric(rf$predicted))

# List the importance of the variables.
rn <- round(randomForest::importance(rf), 2)
head(rn[order(rn[,1], decreasing=TRUE),])

library(verification)
aucc <- verification::roc.area(as.integer(as.factor(data[ind==1, ncol(data)]))-1,
                               rf$votes[,2])$A
verification::roc.plot(as.integer(as.factor(data[ind==1, ncol(data)]))-1,
                       rf$votes[,2], main="")
legend("bottomright", bty="n",
       sprintf("Area Under the Curve (AUC) = %1.3f", aucc))
title(main="OOB ROC Curve Random Forest A0",
      sub=paste("Rattle", format(Sys.time(), "%Y-%b-%d %H:%M:%S"), Sys.info()["user"]))


# ROC of validation model

newdata<-rf_data_birnbaum_t
newdata$class<-samples_Birnbaum$Group

library(ROCR)

# ROC Curve: requires the ggplot2 package.
library(ggplot2, quietly=TRUE)

# Generate an ROC Curve for the rf model on A0 [validate].

pr <- predict(rf, newdata=newdata,
              type    = "prob")[,2]

# Remove observations with missing target.

no.miss   <- na.omit(newdata$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr[-miss.list], no.miss)
} else
{
  pred <- prediction(pr, no.miss)
}

pe <- performance(pred, "tpr", "fpr")
au <- performance(pred, "auc")@y.values[[1]]
pd <- data.frame(fpr=unlist(pe@x.values), tpr=unlist(pe@y.values))
p <- ggplot(pd, aes(x=fpr, y=tpr))
p <- p + geom_line(colour="red")
p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")
p <- p + ggtitle("ROC Curve Random Forest Validation Tissue")
p <- p + theme(plot.title=element_text(size=10))
p <- p + geom_line(data=data.frame(), aes(x=c(0,1), y=c(0,1)), colour="grey")
p <- p + annotate("text", x=0.50, y=0.00, hjust=0, vjust=0, size=5,
                  label=paste("AUC =", round(au, 2)))
print(p)

# Calculate the area under the curve for the plot.


# Remove observations with missing target.

no.miss   <- na.omit(newdata$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr[-miss.list], no.miss)
} else
{
  pred <- prediction(pr, no.miss)
}
performance(pred, "auc")



# ROC of test model

test_data<-test

library(ROCR)

# ROC Curve: requires the ggplot2 package.
library(ggplot2, quietly=TRUE)

# Generate an ROC Curve for the rf model on A0 [validate].

pr_test <- predict(rf, newdata=test_data,
              type    = "prob")[,2]

# Remove observations with missing target.

no.miss   <- na.omit(test_data$class)
miss.list <- attr(no.miss, "na.action")
attributes(no.miss) <- NULL

if (length(miss.list))
{
  pred <- prediction(pr_test[-miss.list], no.miss)
} else
{
  pred <- prediction(pr_test, no.miss)
}

pe <- performance(pred, "tpr", "fpr")
au <- performance(pred, "auc")@y.values[[1]]
pd <- data.frame(fpr=unlist(pe@x.values), tpr=unlist(pe@y.values))
p <- ggplot(pd, aes(x=fpr, y=tpr))
p <- p + geom_line(colour="red")
p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")
p <- p + ggtitle("ROC Curve Random Forest Test Set")
p <- p + theme(plot.title=element_text(size=10))
p <- p + geom_line(data=data.frame(), aes(x=c(0,1), y=c(0,1)), colour="grey")
p <- p + annotate("text", x=0.50, y=0.00, hjust=0, vjust=0, size=5,
                  label=paste("AUC =", round(au, 2)))
print(p)




## PLOT HEATMAPS FOR MCTs and ECMs IN CENTRAL CLUSTER FOR THE PRIMARY and BIRNBAUM DATASETS

# Heatmap for the primary dataset
source("./heatmap3.R")

# Data input - source is normalised and scaled (z-scored) transcriptomic data 
data_primary<-read.csv("./GSE93326_matched plus additional stroma_scaledlogTMM_data.csv", sep = ",", header=T)
data_birnbaum<-read.table("./Birnbaum_matched epithelium stroma_scaledlogTMM_data.csv", sep=",", header=T)

# Read the list of centrally clustered MCTs and ECMs
MCT_ECM_cluster<-read.csv("./MCT_ECM_names_combined.csv", sep=",", header=T)

# Subset MCTs from the MCT_ECM_cluster
MCTs<-dplyr::filter(MCT_ECM_cluster, grepl("SLC16", Gene.name))

#Subset the data_primary to capture only MCTs in central cluster
MCT_TMM_scaled<-data_primary[match(MCTs$Gene.stable.ID,rownames(data_primary)),]

#Convert gene IDs to MCT names
rownames(MCT_TMM_scaled)<-MCTs$Gene.name[match(MCTs$Gene.stable.ID,rownames(MCT_TMM_scaled))]

#Generate heat map
hm_primary<-as.matrix(MCT_TMM_scaled)

## Assign colours for tissue type and study ID. Ref: https://www.biostars.org/p/18211/
library(plyr)
tissue_colour=revalue(Group$Group, c("Stroma"="red3", "Epithelium"="seagreen")) ## ref for colors: http://applied-r.com/r-color-tables/
clab=tissue_colour
clab<-as.data.frame(clab)
colnames(clab)<-"Stroma vs Epithelium"
clab<-as.matrix(clab)

# Plot heatmap for the MCT cluster in primary dataset using the log-scaled TMM values
library(gplots)
heatmap.3(hm_primary, breaks=seq(-3,3,0.5), col=bluered, ColSideColors = clab, ColSideColorsSize=2, Rowv=TRUE, key =TRUE, scale="none",trace="none", dendrogram="both", cexRow=1, cexCol=1, density.info="none", margin=c(10,15)) 
legend("topright",legend=c("Stroma","Epithelium"),
       fill=c("red3","seagreen"), border=FALSE, bty="n", y.intersp = 1, cex=1)

# Use this when grouping stroma separately than epithelium. 
#ref for grouping of columns/samples: https://biocorecrg.github.io/CRG_RIntroduction/heatmap-2-function-from-gplots-package.html
tiff("Expression heatmap of MCTs in central cluster in GSE93326 dataset_21Feb22.tiff", units = "in", width = 10, height = 7, res=300)

heatmap.3(hm_primary, breaks=seq(-3,3,0.5), col=bluered, ColSideColors = clab, ColSideColorsSize=2, Rowv=TRUE, key =TRUE, keysize = 2, scale="none",trace="none", dendrogram="row", cexRow=2.5, cexCol=2.5, density.info="none", margin=c(16.5,15), Colv = FALSE, labRow=as.expression(lapply(rownames(hm_primary), function(a) bquote(italic(.(a))))),) #https://stackoverflow.com/questions/57207651/make-row-labels-italic-in-heatmap-2) 
legend("topright",legend=c("Stroma","Epithelium"),
       fill=c("red3","seagreen"), border=FALSE, bty="n", y.intersp = 1, cex=1)
dev.off()

# Subset ECMs from the MCT_ECM_cluster
ECMs<-MCT_ECM_cluster[-c(85,95,119,210,298,324,376),]

#Subset the data_primary to capture only ECMs in central cluster
ECM_TMM_scaled<-data_primary[match(ECMs$Gene.stable.ID,rownames(data_primary)),]

#Convert gene IDs to ECM names
rownames(ECM_TMM_scaled)<-ECMs$Gene.name[match(ECMs$Gene.stable.ID,rownames(ECM_TMM_scaled))]

#Generate heat map
hm_ECM_primary<-as.matrix(ECM_TMM_scaled)

## Assign colours for tissue type and study ID. Ref: https://www.biostars.org/p/18211/
library(plyr)
tissue_colour=revalue(Group$Group, c("Stroma"="red3", "Epithelium"="seagreen")) ## ref for colors: http://applied-r.com/r-color-tables/
clab=tissue_colour
clab<-as.data.frame(clab)
colnames(clab)<-"Stroma vs Epithelium"
clab<-as.matrix(clab)

# Use this when grouping stroma separately than epithelium. 
#ref for grouping of columns/samples: https://biocorecrg.github.io/CRG_RIntroduction/heatmap-2-function-from-gplots-package.html
heatmap.3(hm_ECM_primary, breaks=seq(-2,2,0.5), col=bluered, ColSideColors = clab, ColSideColorsSize=2, Rowv=TRUE, key =TRUE, scale="none",trace="none", dendrogram="row", cexRow=1, cexCol=1, density.info="none", margin=c(10,15), Colv = FALSE) 
legend("topright",legend=c("Stroma","Epithelium"),
       fill=c("red3","seagreen"), border=FALSE, bty="n", y.intersp = 1, cex=1)


#Subset the data_birnbaum to capture only MCTs in central cluster
MCT_TMM_scaled_birnbaum<-data_birnbaum[na.omit(match(MCTs$Gene.name,rownames(data_birnbaum))),]

#Generate heat map
hm_MCT_birnbaum<-as.matrix(MCT_TMM_scaled_birnbaum)

## Assign colours for tissue type and study ID. Ref: https://www.biostars.org/p/18211/
library(plyr)
tissue_colour=revalue(samples_Birnbaum$Group, c("Stroma"="red3", "Epithelium"="seagreen")) ## ref for colors: http://applied-r.com/r-color-tables/
clab=tissue_colour
clab<-as.data.frame(clab)
colnames(clab)<-"Stroma vs Epithelium"
clab<-as.matrix(clab)

# Use this when grouping stroma separately than epithelium. 
tiff("Expression heatmap of MCTs in Birnbaum dataset which were centrally clustered in GSE93326 dataset_21Feb22.tiff", units = "in", width = 10, height = 7, res=300)

heatmap.3(hm_MCT_birnbaum, breaks=seq(-2,2,0.5), col=bluered, ColSideColors = clab, ColSideColorsSize=2, Rowv=TRUE, key =TRUE, keysize =2, scale="none",trace="none", dendrogram="row", cexRow=2.5, cexCol=2.5, density.info="none", margin=c(10,15), Colv = FALSE, labRow=as.expression(lapply(rownames(hm_MCT_birnbaum), function(a) bquote(italic(.(a))))),) #https://stackoverflow.com/questions/57207651/make-row-labels-italic-in-heatmap-2
legend("topright",legend=c("Stroma","Epithelium"),
       fill=c("red3","seagreen"), border=FALSE, bty="n", y.intersp = 1, cex=1)
dev.off()

#Subset the data_birnbaum to capture only ECMs in central cluster
ECM_TMM_scaled_birnbaum<-data_birnbaum[na.omit(match(ECMs$Gene.name,rownames(data_birnbaum))),]

#Generate heat map
hm_ECM_birnbaum<-as.matrix(ECM_TMM_scaled_birnbaum)

## Assign colours for tissue type and study ID. Ref: https://www.biostars.org/p/18211/
library(plyr)
tissue_colour=revalue(samples$Group, c("Stroma"="red3", "Epithelium"="seagreen")) ## ref for colors: http://applied-r.com/r-color-tables/
clab=tissue_colour
clab<-as.data.frame(clab)
colnames(clab)<-"Stroma vs Epithelium"
clab<-as.matrix(clab)

# Use this when grouping stroma separately than epithelium. 
heatmap.3(hm_ECM_birnbaum, breaks=seq(-2,2,0.5), col=bluered, ColSideColors = clab, ColSideColorsSize=2, Rowv=TRUE, key =TRUE, scale="none",trace="none", dendrogram="row", cexRow=1, cexCol=1, density.info="none", margin=c(10,15), Colv = FALSE) 
legend("topright",legend=c("Stroma","Epithelium"),
       fill=c("red3","seagreen"), border=FALSE, bty="n", y.intersp = 1, cex=1)



# Direct path analysis used for Birnbaum (validation) dataset

setwd("F:/Latif Lab back up/AYSHE/PDAC/Birnbaum et al._Validation study")

# Data input - source is normalised and scaled (z-scored) transcriptomic data
data<-read.csv("./Birnbaum_matched epithelium stroma_scaledlogTMM_data.csv", sep = ",", header=T)

# Create a vector of Ensembl IDs from the row names
rownames<-rownames(data)
head(rownames)

# Create a data matrix
matrix_data <-data.matrix(data)
head(matrix_data)

## Plot matrix_data
plot(matrix_data)

# Read the differentially expressed genes from Maurer (primary) dataset 
DEGs<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/Genes with 5% FDR-adjusted p-values_matched plus additional stroma.csv", sep = ",", header=T)

# Assign the first column of DEGs as row names
rownames(DEGs)<-DEGs$X

# Create a vector of row names from DEGs and convert it to a data frame
DEGs_rowname<-rownames(DEGs)
DEGs_rowname<-as.data.frame(DEGs_rowname)
DEGs2<-cbind(DEGs_rowname,DEGs)
rownames(DEGs2)<-1:nrow(DEGs2)


# Filter out the ECMs of interest from DEGs2

# First identify integrins and keratins from the list as they are not part of the matrisome geneset
Integrins<-dplyr::filter(DEGs2, grepl("ITG", Gene_name))
Keratins<-dplyr::filter(DEGs2, grepl("KRT", Gene_name))

# List of matrisome and matrisome associated genes from Molecular Signature Database (MSigDB)
# Ref paper:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5013529/
# MSigDB ref for matrisome genes:https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=NABA_MATRISOME

# Read matrisome geneset: Ensemble of genes encoding extracellular matrix and extracellular matrix-associated proteins
# The list includes members mapped to 1026 genes
Matrisome<-read.table("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS/Matrisome_geneset.txt", sep = ",", header=T)

# Change the column name in Matrisome dataset to be able to combine with Integrins and Keratins
names(Matrisome)[1]<- "Gene_name"

# Create a new column called Matrisome.match in DEGs2 to match the differentially-expressed genes with matrisome geneset
DEGs2['Matrisome.match'] <- NA
DEGs2$Matrisome.match <- Matrisome$Gene_name[match(DEGs2$Gene_name, Matrisome$Gene_name)]

# Drop NAs from the list
library(tidyr)
DEGs_matrisome<-DEGs2 %>% drop_na()
DEGs_matrisome<-DEGs_matrisome[,1:7] #to match the column numbers of integrins and keratins dataset to combine

# Combine DE-matrisome gene set with Integrins and Keratins
All_ECM_genes<- rbind(Integrins, Keratins, DEGs_matrisome)

# Filter out the MCTs from DEGs2 
MCTs<-dplyr::filter(DEGs2, grepl("SLC16", Gene_name))
MCTs<-MCTs[,-8]

# Combine ECMs and MCTs
MCT_ECM<-rbind(MCTs,All_ECM_genes)

# Assign the names of DE-matrisome genes to "gene_list" to be used when creating a correlation matrix downstream
#gene_list<-MCT_ECM$DEGs_rowname
gene_list<-MCT_ECM$Gene_name

# Generate a correlation matrix 
# Correlation matrix based on complete dataset (epithelium and stroma). Correlating MCT_ECM against each other
cor_matrix<-cor(t(matrix_data[na.omit(match(gene_list,rownames(matrix_data))),]), # omit NA matches for instances where genes not found in matrix
                t(matrix_data[na.omit(match(gene_list,rownames(matrix_data))),])) 

# Examine the distribution of correlation matrix values for normality with a histogram
hist(cor_matrix, main = paste(""), xlab = "Correlation r-values", col = "darkgrey", cex.lab=1.5, cex.axis=1.5) #ref: https://stackoverflow.com/questions/4241798/how-to-increase-font-size-in-a-plot-in-r and https://www.datamentor.io/r-programming/histogram/
hist(cor_matrix[1:7,8:410])


# Calculate standard deviation
sd(cor_matrix)

# Generate absolute values for all correlation r-values. This matrix will be binarized
binary<-abs(cor_matrix) 

# Set r-value thresold to SD for binarisation
thresh<-0.33 # set r-value threshold as SD for binarisation
binary[which(binary>thresh)]<-1 # set any values greater than the threshold to 1
binary[which(binary!=1)]<-0 # set any values which aren't 1s as 0s 

# Hypernetwork generation and visualisation
hyp<-binary%*%t(binary) # Generate the hypernetwork by multiplying the binary matrix by the transpose of itself. ('%*%' is the operator for matrix multiplication ). This represents the adjacency matrix

# Generate heatmap (source for margin correction: https://stackoverflow.com/questions/21427863/how-to-add-more-margin-to-a-heatmap-2-plot-with-the-png-device)
# Further info on heatmap.2 function: https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2

#hm<-heatmap.2(hyp,trace="none") # generate a heatmap from the hypernetwork matrix, save the heatmap object as it contains dendrograms.Exclude the trace which aims to separates rows and columns
hm<-heatmap.2(hyp, cexRow=1, cexCol=1, margins=c(10,10), srtCol=45, trace="none")

# Customise heatmap axes labels
hm<-heatmap.2(hyp, cexRow=1, cexCol=1, margins=c(2,2), trace="none", xlab = "Differentially expressed MCT and ECM genes", ylab = "Differentially expressed MCT and ECM genes", labRow = FALSE, labCol = FALSE)

# Conduct direct path analysis

# S= (G - I + D((G-I)G))G^-1

# S = direct network (?)
# G = Hypernetwork ===>This equals to the hypernetwork/adjacency matrix
# I = Identity Matrix (diagonal ones) #poss also diagonal of hypernetwork?
# D = diagonal only

G<-cor_matrix # G is the correlation matrix
I<-matrix(0,nrow=nrow(G),ncol=ncol(G));diag(I)<-diag(G) #I is the identity matrix but the diagonals are set to the diagonals of G rather than 1s


GIG<-(G-I)%*%G
DGIG<-matrix(0,nrow=nrow(G),ncol=ncol(G));diag(DGIG)<-diag(GIG)

inv_G<-MASS::ginv(G)

S = (G - I + DGIG)%*%inv_G

colnames(S)<-rownames(S)

# Generate txt files to be imported in Cytoscape
write.table(S,"MCT_ECM_Silenced_Hypernetwork_based on cor matrix_BIRNBAUM.txt",sep="\t",col.names=NA)
write.table(G,"MCT_ECM_Comparison_Hypernetwork based on cor matrix_BIRNBAUM.txt",sep="\t",col.names=NA)

hist(S)
hist(S, main = paste(""), xlab = "Correlation r-values in S", col = "darkgrey", cex.lab=1.5, cex.axis=1.5) #ref: https://stackoverflow.com/questions/4241798/how-to-increase-font-size-in-a-plot-in-r and https://www.datamentor.io/r-programming/histogram/

hist(S[1:7,8:410])

# Run statistics
zscore <- (0.242101588104453 - mean(S))/sd(S) # The input value here is the relationship directness score for a pair of genes
pval<-2*pnorm(zscore, mean = 0, sd = 1, lower.tail = FALSE)

zscore_bottom <- (-0.004194762 - mean(S))/sd(S) # The input value here is the relationship directness score for a pair of genes
pval_bottom<-2*pnorm(zscore_bottom, mean = 0, sd = 1, lower.tail = T)


# Customise heatmap axes labels
hm_sil<-heatmap.2(S, cexRow=1, cexCol=1, margins=c(2,2), trace="none", xlab = "Differentially expressed MCT and ECM genes", ylab = "Differentially expressed MCT and ECM genes", labRow = FALSE, labCol = FALSE)

# Generate directly comparable heatmaps
hm_hyp<-heatmap.2(G,trace="none")

G_re<-G[rev(hm_hyp$rowInd),hm_hyp$rowInd]
heatmap.2(G_re,trace="none",Rowv = F,Colv = F)

S_re<-S[rev(hm_hyp$rowInd),hm_hyp$rowInd]
heatmap.2(S_re,trace="none",Rowv = F,Colv = F)
