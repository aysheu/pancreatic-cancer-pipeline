# The object 'data' is a matrix genes x samples
# The object 'gene_list' is a list of genes (eg DEGs) which form the nodes of the hypernetwork
# The object 'rf_variable' is a variable separating the samples used for random forest

setwd("F:/Latif Lab back up/AYSHE/PDAC/Birnbaum et al._Validation study")

library(gplots) # heatmap.2 function required
library(plyr)

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
gene_list<-All_ECM_genes$Gene_name
#write.csv(All_ECM_genes, "All ECM-related DEGs in stroma.csv")

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
thresh<-0.24 # set r-value threshold as SD for binarisation
binary[which(binary>thresh)]<-1 # set any values greater than the threshold to 1
binary[which(binary!=1)]<-0 # set any values which aren't 1s as 0s 

# Hypernetwork generation and visualisation
hyp<-binary%*%t(binary) # Generate the hypernetwork by multiplying the binary matrix by the transpose of itself. ('%*%' is the operator for matrix multiplication ). This represents the adjacency matrix

# Generate heatmap (source for margin correction: https://stackoverflow.com/questions/21427863/how-to-add-more-margin-to-a-heatmap-2-plot-with-the-png-device)
# Further info on heatmap.2 function: https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2

#hm<-heatmap.2(hyp,trace="none") # generate a heatmap from the hypernetwork matrix, save the heatmap object as it contains dendrograms.Exclude the trace which aims to separates rows and columns
hm<-heatmap.2(hyp, cexRow=1, cexCol=1, margins=c(10,10), srtCol=45, trace="none")

# Customise heatmap axes labels
tiff("Hypernetwork heatmap of ECMs_r-value threshold of 0.24_Birnbaum dataset.tiff", units = "in", width = 10, height = 7, res=300)

hm<-heatmap.2(hyp, cexRow=1, cexCol=1, margins=c(2,2), trace="none", xlab = "Differentially expressed ECM genes in validation dataset", ylab = "Differentially expressed ECM genes in validation dataset", labRow = FALSE, labCol = FALSE)
dev.off()

# Create a class hclust
dendrogram<-as.hclust(hm$rowDendrogram) 

# Cut the dendrogram to generate 2 groups (first dendrogram split)
ct<-cutree(dendrogram,k = 2) 

# Extract the names of the 1st cluster or central cluster
central_cluster_genes<-names(ct[which(ct==1)]) #The "not equal to" operator is defined by "!=")
print(central_cluster_genes)
write.csv (central_cluster_genes, "CCG_ECM_Birnbaum.csv")

# Extract the differentially expressed genes not in central cluster (2nd cluster)
Not_central_cluster_genes<-names(ct[which(ct!=1)])
print(Not_central_cluster_genes)
write.csv(Not_central_cluster_genes, "NCG_ECM_Birnbaum.csv")

# Extract the 'Galois' from the hypernetwork
# Use the central cluster gene list to subset the binarized correlation matrix 
galois <- binary[match(central_cluster_genes, rownames(binary)),]
#galoisNC<-binary[match(Not_central_cluster_genes, rownames(binary)),]

# Read rownames from the primary matrix data
rownames_primary<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Birnbaum et al._Validation study/rownames of matrix data from primary dataset_Gene names from Ensembl.csv", sep = ",", header = T)
MCT_list<-dplyr::filter(rownames_primary, grepl("SLC16", Gene.name))


# Compute sum of each column from the galois (total # correlations between the wider trans. and the central cluster)
percent_corr<-(colSums(galois)/nrow(galois))*100 # % correlations
MCT_percentage<-percent_corr[match(MCT_list$Gene.name,names(percent_corr))]
MCT_percentage<-na.omit(MCT_percentage)
write.csv(MCT_percentage, "F:/Latif Lab back up/AYSHE/PDAC/Birnbaum et al._Validation study/MCT percentage.csv")

hist(percent_corr)
hist(MCT_percentage)

nonMCT_percentage<-percent_corr[-na.omit(match(MCT_list$Gene.name,names(percent_corr)))]
wilcox.test(MCT_percentage,nonMCT_percentage)

hist(nonMCT_percentage)
mean(nonMCT_percentage)

mean(MCT_percentage)

length(MCT_percentage[which(MCT_percentage>75)])

length(MCT_percentage[which(MCT_percentage==100)])