# The object 'data' is a matrix genes x samples
# The object 'gene_list' is a list of genes (eg DEGs) which form the nodes of the hypernetwork

setwd("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS")

library(gplots) # heatmap.2 function required
library(dplyr)

# Data input - source is normalised and scaled (z-scored) transcriptomic data 
data<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/GSE93326_matched plus additional stroma_scaledlogTMM_data.csv", sep = ",", header=T)

# Create a vector of Ensembl IDs from the row names
rownames<-rownames(data)
head(rownames)

# Create a data matrix
matrix_data <-data.matrix(data)
head(matrix_data)

# Plot matrix_data
plot(matrix_data)

# Read the differentially expressed genes 
DEGs<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/Genes with 5% FDR-adjusted p-values_matched plus additional stroma.csv", sep = ",", header=T)

# Assign the first column of DEGs as row names
rownames(DEGs)<-DEGs$X

# Create a vector of row names from DEGs and convert it to a data frame
DEGs_rowname<-rownames(DEGs)
DEGs_rowname<-as.data.frame(DEGs_rowname)
DEGs2<-cbind(DEGs_rowname,DEGs) # Column wise combination of DEGs_rowname and DEGs data frames
rownames(DEGs2)<-1:nrow(DEGs2) # Give numeric values to row names

# Filter out the MCTs from DEGs2 
MCTs<-dplyr::filter(DEGs2, grepl("SLC16", Gene_name))

# This assigns the names of DEGs to "gene_list" to be used when creating a correlation matrix downstream 
gene_list<-MCTs$DEGs_rowname

# Generate a correlation matrix 
# Correlation matrix based on complete dataset (epithelium and stroma)
cor_matrix<-cor(t(matrix_data[na.omit(match(gene_list,rownames(matrix_data))),]), #correlate genes which match the list against all other genes (omit NA matches for instances where genes not found in matrix)
                t(matrix_data[-na.omit(match(gene_list,rownames(matrix_data))),])) #NB the '-' in the square brackets on this line but absent from the previous line. This means exclude these entries (entries which match the gene list)

# Examine the distribution of correlation matrix values for normality with a histogram
hist(cor_matrix, main = paste(""), xlab = "Correlation r-values", col = "darkgrey", cex.lab=1.5, cex.axis=1.5) #ref: https://stackoverflow.com/questions/4241798/how-to-increase-font-size-in-a-plot-in-r

# Calculate standard deviation
sd(cor_matrix) 

# Generate absolute values for all correlation r-values. This matrix will be binarized
binary<-abs(cor_matrix) 

# Set r-value thresold to SD for binarisation
thresh<-0.2
binary[which(binary>thresh)]<-1 # Set any values greater than the threshold to 1
binary[which(binary!=1)]<-0 # Set any values which aren't 1s as 0s 

# Hypernetwork generation and visualisation
hyp<-binary%*%t(binary) # Generate the hypernetwork by multiplying the binary matrix by the transpose of itself. ('%*%' is the operator for matrix multiplication ). This represents the adjacency matrix

# Generate heatmap (source for margin correction: https://stackoverflow.com/questions/21427863/how-to-add-more-margin-to-a-heatmap-2-plot-with-the-png-device)
# Further info on heatmap.2 function: https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2
hm<-heatmap.2(hyp, cexRow=1, cexCol=1, margins=c(10,10), srtCol=45, trace="none") # Generate a heatmap from the hypernetwork matrix, save the heatmap object as it contains dendrograms.Exclude the trace which aims to separates rows and columns

# Customise heatmap axes labels
tiff("Hypernetwork heatmap of MCTs_r-value threshold of 0.2 as SD_10Feb22.tiff", units = "in", width = 10, height = 7, res=300)

hm<-heatmap.2(hyp, cexRow=1, cexCol=1, margins=c(2,2), trace="none", xlab = "Differentially expressed MCT genes", ylab = "Differentially expressed MCT genes", labRow = FALSE, labCol = FALSE)
dev.off()

# Create a class hclust
dendrogram<-as.hclust(hm$rowDendrogram) 

# Cut the dendrogram to generate 2 groups (first dendrogram split)
ct<-cutree(dendrogram,k = 2) 

# Extract the names of the 1st cluster or central cluster
central_cluster_genes<-names(ct[which(ct==1)]) 
write.csv(central_cluster_genes, "CCG_MCTs_stroma_using_r-value_threshold_of_0.2.csv")
print(central_cluster_genes)

# Extract the differentially expressed genes not in central cluster (2nd cluster)
Not_central_cluster_genes<-names(ct[which(ct!=1)])
print(Not_central_cluster_genes)
write.csv(Not_central_cluster_genes, "NCG_MCTs_stroma_using_r-value_threshold_of_0.2.csv")

# Extract the 'Galois' from the hypernetwork
# Use the central and not central cluster gene list to subset the binarized correlation matrix 
galois <- binary[match(central_cluster_genes, rownames(binary)),]
galoisNC<-binary[match(Not_central_cluster_genes, rownames(binary)),]

# Identify which of the columns of this matrix have 1 in every row (i.e. which of the non-DE genes are correlated with every one of the central cluster genes?)
galois2 <- galois[, which(colSums(galois) == nrow(galois))] # We use colSums and nrow here as a column filled with 1s should sum to the number of rows in that column; which tells us which of the columns this is true for.

# Identify which of the columns of this matrix have 1 for 90% of the central cluster genes
galois3<-galois[, which(colSums(galois) > nrow(galois)*0.9)]

# Identify which of the columns of this matrix have 1 for either all or 90% of the not central cluster genes
galoisNC<- galoisNC[, which(colSums(galoisNC) == nrow(galoisNC))]
galoisNC<- galoisNC[, which(colSums(galoisNC) > nrow(galoisNC)*0.9)]  

nonDE_genes1<-colnames(galois2) #names of nonDE genes that correlate with all MCTs in central cluster
nonDE_genes2<-colnames(galois3) #names of nonDE genes that correlate with 90% of MCTs in central cluster. This gives essentially the same list as above
nonDE_genes3<-colnames(galoisNC) #names of nonDE genes that correlate with genes in non-central cluster

#Look at folder "Analysis of MCTs using r-value threshold of 0.2 for binarisation" under Hypernetwork analysis folder
write.csv(nonDE_genes1, "non-DE genes correlated with all MCTs in central cluster_stroma.csv")
write.csv(nonDE_genes3, "non-DE genes correlated with two MCTs in non-central cluster_stroma.csv")

# Read the nonDE genes correlated with all MCTs in central cluster
# The Ensembl IDs in nonDE_genes1 were converted to gene names using Ensembl's online BioMart tool
nonDE_names1<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS/Analysis of MCTs using r-value threshold of 0.2 for binarisation/nonDE genes correlated with all MCTs in CC_stroma_Ensembl IDs.csv", sep = ",", header=T)

# List of integrins correlated with all MCTs
Integrins<-dplyr::filter(nonDE_names1, grepl("ITG", Gene.name)) #From the list that correlate with all MCTs

# List of keratins correlated with all MCTs
Keratins<-dplyr::filter(nonDE_names1, grepl("KRT", Gene.name)) #From the list that correlate with all MCTs

# List of matrisome and matrisome associated genes from Molecular Signature Database (MSigDB)
# ref paper:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5013529/
# MSigDB ref for matrisome genes:https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=NABA_MATRISOME

# Read matrisome geneset: Ensemble of genes encoding extracellular matrix and extracellular matrix-associated proteins
# The list includes members mapped to 1026 genes
Matrisome<-read.table("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS/Matrisome_geneset.txt", sep = ",", header=T)

# Add integrins and keratins to the matrisome geneset
Integrins<-Integrins[2]
Keratins<-Keratins[2] 

# Change the column name in Matrsime dataset to be able to combine with Integrins and Keratins
names(Matrisome)[1]<- "Gene.name"

# Matrisome_with_ITG.KRTs<-rbind(Matrisome,Integrins,Keratins) 
Matrisome_with_ITG.KRTs<-rbind(Matrisome,Integrins,Keratins) 
write.csv(Matrisome_with_ITG.KRTs, "All ECM genes from the primary dataset.csv")

# Match the gene names of the "nonDE_name" dataset with the genes in the matrisome dataset
nonDE_names1['Matrisome.match'] <- NA
nonDE_names1$Matrisome.match <- Matrisome_with_ITG.KRTs$Gene.name[match(nonDE_names1$Gene.name, Matrisome_with_ITG.KRTs$Gene.name)]
library(tidyr)
nonDE_names_final<-nonDE_names1 %>% drop_na()
write.csv(nonDE_names_final, "Non-DE Matrisome genes correlated with ALL MCTs.csv")
