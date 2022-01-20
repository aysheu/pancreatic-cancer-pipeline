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
thresh<-0.28 # set r-value threshold as SD for binarisation
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
I<-matrix(0,nrow=nrow(G),ncol=ncol(G));diag(I)<-1 #I is the identity matrix but the diagonals are set to the diagonals of G rather than 1s


GIG<-(G-I)*G
DGIG<-matrix(0,nrow=nrow(G),ncol=ncol(G));diag(DGIG)<-diag(GIG)

inv_G<-MASS::ginv(G)

S = (G - I + DGIG)*inv_G

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

