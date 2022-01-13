# The object 'data' is a matrix genes x samples
# The object 'gene_list' is a list of genes (eg DEGs) which form the nodes of the hypernetwork
# The object 'rf_variable' is a variable separating the samples used for random forest

setwd("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS")

library(gplots) # heatmap.2 function required
library(plyr)

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

# Assign the names of DE-matrisome genes to "gene_list" to be used when creating a correlation matrix downstream
gene_list<-All_ECM_genes$DEGs_rowname
write.csv(All_ECM_genes, "All ECM-related DEGs in stroma.csv")

# Generate a correlation matrix 
# Correlation matrix based on complete dataset (epithelium and stroma)
cor_matrix<-cor(t(matrix_data[na.omit(match(gene_list,rownames(matrix_data))),]), #correlate genes which match the list against all other genes (omit NA matches for instances where genes not found in matrix)
                t(matrix_data[-na.omit(match(gene_list,rownames(matrix_data))),])) #NB the '-' in the square brackets on this line but absent from the previous line. This means exclude these entries (entries which match the gene list)

# Examine the distribution of correlation matrix values for normality with a histogram
hist(cor_matrix, main = paste(""), xlab = "Correlation r-values", col = "darkgrey", cex.lab=1.5, cex.axis=1.5) #ref: https://stackoverflow.com/questions/4241798/how-to-increase-font-size-in-a-plot-in-r and https://www.datamentor.io/r-programming/histogram/

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
hm<-heatmap.2(hyp, cexRow=1, cexCol=1, margins=c(2,2), trace="none", xlab = "Differentially expressed ECM genes", ylab = "Differentially expressed ECM genes", labRow = FALSE, labCol = FALSE)

# Create a class hclust
dendrogram<-as.hclust(hm$rowDendrogram) 

# Cut the dendrogram to generate 2 groups (first dendrogram split)
ct<-cutree(dendrogram,k = 2) 

# Extract the names of the 1st cluster or central cluster
central_cluster_genes<-names(ct[which(ct!=1)]) #The "not equal to" operator is defined by "!=")
print(central_cluster_genes)
write.csv (central_cluster_genes, "CCG_ECM_stroma.csv")

# Extract the differentially expressed genes not in central cluster (2nd cluster)
Not_central_cluster_genes<-names(ct[which(ct==1)])
print(Not_central_cluster_genes)
write.csv(Not_central_cluster_genes, "NCG_ECM_stroma.csv")

# Extract the 'Galois' from the hypernetwork
# Use the central cluster gene list to subset the binarized correlation matrix 
galois <- binary[match(central_cluster_genes, rownames(binary)),]
galoisNC<-binary[match(Not_central_cluster_genes, rownames(binary)),]

# Identify which of the columns of this matrix now have a 1 in every row (i.e. which of the non-DE genes are correlated with every one of the central cluster genes?)
galois2 <- galois[, which(colSums(galois) == nrow(galois))] # we use colSums and nrow here as a column filled with 1s should sum to the number of rows in that column; which tells us which of the columns this is true for.
galoisNC1<- galoisNC[, which(colSums(galoisNC) == nrow(galoisNC))]
galoisNC2<- galoisNC[, which(colSums(galoisNC) > nrow(galoisNC)*0.9)]  

nonDE_genes1<-colnames(galoisNC1) #names of nonDE genes that correlate with ALL ECMs in NC cluster
nonDE_genes2<-colnames(galoisNC2) #names of nonDE genes that correlate with 90% of ECMs in NC cluster

write.csv(nonDE_genes1, "non-DE genes correlated with ALL Non-centrally clustered ECMs_stroma.csv")
write.csv(nonDE_genes2, "non-DE genes correlated with 90% of Non-centrally clustered ECMs_stroma.csv")

# Match the list of ECMs in "non-central_cluster" (that regulate MCTs) with ECMs found to be regulated by MCTs
# ECMs that regulate MCTs
ECM_regulating_MCT<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS/Analysis of ECMs using r-value threshold of 0.28 for binarisation/NCG_ECM_stroma_Ensembl IDs.csv",sep=",", header=T)

# ECMs that are regulated by MCTs (ECMs corr. with All MCTs which are equal to ECMs corr. with 90% of MCTs)
ECM_regulatedBy_MCT<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS/Analysis of MCTs using r-value threshold of 0.2 for binarisation/Non-DE Matrisome genes correlated with ALL MCTs.csv", sep=",", header=T)

# Remove below cols to match with the ECM_regulating_MCT dataframe
ECM_regulatedBy_MCT<-ECM_regulatedBy_MCT[,-1]
ECM_regulatedBy_MCT<-ECM_regulatedBy_MCT[,-3]

# Match the names of ECMs in both datasets: order matters, you need to select names from the smaller dataset to be listed in the overlap dataset, hence "ECM_regulatedBy_MCT", otherwise it won't produce exact matches
ECM_overlap<-ECM_regulatedBy_MCT$Gene.name[match(ECM_regulating_MCT$Gene.name,ECM_regulatedBy_MCT$Gene.name)]
View(ECM_overlap)
ECM_overlap<-as.data.frame(ECM_overlap)
ECM_overlap<-na.omit(ECM_overlap) # Remove NAs

# Add the Ensembl IDs to the ECM_overlap dataframe
ECM_overlap$Gene.stable.ID<-ECM_regulating_MCT$Gene.stable.ID[match(ECM_overlap$ECM_overlap,ECM_regulating_MCT$Gene.name)]
write.csv(ECM_overlap, "ECMs_overlap_between_376regulatingMCTs_255regulatedByMCTs_at90percent.csv")
