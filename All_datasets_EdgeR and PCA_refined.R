library(edgeR) #https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
library(biomaRt)
library(plyr)
library(gplots)
library(mixOmics)

setwd("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/Combined datasets/Lin, Maurer, Kirby and Rashid_tumor vs healthy/ALL dataset including tumor epithelium and stroma")

#Differential gene expression (DGE) analysis

#Create a dataframe including information from all samples
samples<-read.table("./ALL dataset_design.csv",sep = ",", header=T, na.strings=c(""," ","NA")) #This is OK to use for Combined dataset (tumour vs non-tumour) and Maurer et al. (stroma vs epithelium)

#Create a design matrix for DGE analysis. 
design<-model.matrix(~samples$Group+samples$Study.ID)
rownames(design)<-samples$SampleID

# Read the data table containing all samples and genes
dat<-read.table("./ALL dataset.txt",sep=c("\t"),header=T)
counts<-dat[-1] # Remove the first column of dat which includes the Ensembl gene IDs
counts<-dat[-c(60609:60613), -1] # Omit the last 4 rows of the count table as counts are NA
rownames(counts)<-dat[-c(60609:60613),1] #Assign the first column of dat as rownames to counts table

# Read the data table containing common genes found between the 4 studies: see the R.file "Merging all datasets with common genes"
# This table is only to match the gene names with those in "dat" table
dataset.trimmed<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/PCA analysis of ALL datasets_with filter applied/Revised analysis/Common genes identified among all datasets.csv", sep = ",", header=T)
trimmed.counts<-counts[match(dataset.trimmed$x, rownames(counts)),] #Match the gene IDs of the common genes with that of those from the count table
#write.csv(trimmed.counts, "Raw count data of trimmed dataset.csv", sep=",")
counts<-trimmed.counts

ENSG<-dat[,1] # select list of gene IDs
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #Create a Mart object. This connects to biomart's ensembl database and the specified dataset
res<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'ensembl_gene_id', values = ENSG, mart = ensembl) #Retrieve gene names from the BioMart database
counts<-counts[na.omit(match(res$ensembl_gene_id,rownames(counts))),] #Remove genes with IDs from the counts table which cannot be annotated

##Create DGE object and specify counts object and grouping variable
DGE<-DGEList(counts=counts,group=samples$Group) 

##Filter low expression genes
keep <- filterByExpr(DGE) #ID low expression genes
DGE <- DGE[keep, , keep.lib.sizes=FALSE] # Remove low expression genes. The library size for each sample is recalculated to be the sum of the counts left in the rows after removing low expressed genes

##Calculate effective library sizes from raw library sizes using scaling factors
DGE <- calcNormFactors(DGE) #This gives TMM normalised counts and normalisation factors

#Extract normalised log2TMM data
cpm_data<-cpm(DGE,log=T) 
write.table(cpm_data, "ALL_dataset_tmm_data.csv",sep=",")

## Scale each row (gene) to a mean of zero and a standard deviation one
logCPM <- t(scale(t(cpm_data)))

## To avoid shift to left in CSV file, include "col.names = NA, row.names = TRUE". Adds a column header as "X" which needs to be trimmed later
write.table(logCPM, "ALL_dataset_scaled_logTMM_data_test.csv", sep=",",col.names = NA)

# Estimate dispersions 
DGE <- estimateDisp(DGE, design) 

# Identify differentialy expressed genes
et <- exactTest(DGE, pair=c("Non-tumour", "Tumour")) #Calculate differential expression. This finds genes dif. expressed in tumour vs non-tumour in the combined dataset
et$table$Gene_name<-res$external_gene_name[na.omit(match(rownames(et$table),res$ensembl_gene_id))] # assign gene names

#Extract the list that are significant p< 0.05
p<-topTags(et,n=Inf,adjust.method="BH") #Extract the DEGs with p-value adjusted for multiple testing
keep<-p$table$FDR <=0.05 #Keep only the genes with FDR<0.05
p.table<-p[keep,]

#Write a table for each analysis
write.table(p.table[keep,],"Genes with 5% FDR-adjusted p-values_T vs NT.csv",sep=",")

## PRODUCE HEATMAP FOR TUMOUR AND NON-TUMOUR SAMPLES
source("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/heatmap3.R")

# select 5% of DE-genes and produce heatmap using log-normalised and z-scored (scaled) data
selcpm <- logCPM[rownames(p$table)[p$table$FDR<0.05],1:133]

#Generate heat map
hm<-as.matrix(selcpm)

## Assign colours for tissue type and study ID. Ref: https://www.biostars.org/p/18211/
sample_ids=samples$SampleID
tissue_colour=revalue(samples$Group, c("Non-tumour"="seagreen", "Tumour"="red3")) ## ref for colors: http://applied-r.com/r-color-tables/
tissue_colour<-tissue_colour[1:133]
study_colour=revalue(samples$Study.ID, c("GSE119794"="magenta", "GSE79668"="dodgerblue", "GSE93326"="darkorchid", "GSE131050"="yellow3"))
study_colour<-study_colour[1:133]
clab=cbind(study_colour, tissue_colour)
colnames(clab)=c("Study_ID", "Tissue_type")

# As we are using log-normalised data which has been z-scored, there is no need to normalise the data below

#This one customises the x and y labels as "Samples" and "Genes" and increases the font size of the legend
#The heatmap was further edited manually for axes titles and the legend
heatmap.3(hm, breaks=seq(-1,1,0.2), col=bluered, ColSideColors = clab, ColSideColorsSize=2, Rowv=TRUE, key =TRUE, scale="none",trace="none", dendrogram="both", cexRow=1, cexCol=1, density.info="none", margin=c(15,15), xlab = "Samples", ylab = "Genes", labRow = FALSE, labCol = FALSE, distfun = dist) 
legend("bottomright",legend=c("Non-tumour","Tumour","","GSE119794","GSE79668","GSE93326","GSE131050"),
       fill=c("seagreen","red3","white", "magenta","dodgerblue","darkorchid","yellow3"), border=FALSE, bty="n", y.intersp = 0.8, cex=0.8, xy.coords(legend))


# Principal Component Analysis

#Filter only selected studies from the sample list 
library(dplyr)
samples2<-filter(samples, Study.ID == "GSE119794")
samples3<-filter(samples, Study.ID == "GSE79668")
samples4<-rbind(samples2, samples3)
samples5<- filter(samples, Study.ID =="GSE131050")
samples6<-rbind(samples4, samples5) # used to conduct PCA based on all datasets except for GSE93326 (stromal/epithelial data)

#To conduct PCA for all datasets 
counts<-trimmed.counts
counts<-counts[na.omit(match(res$ensembl_gene_id,rownames(counts))),] #Remove genes with IDs from the counts table which cannot be annotated

# To conduct PCA for all datasets excluding GSE93326
counts2<-trimmed.counts[,-c(72:86)] 
counts2<-counts2[,-c(119:307)] 
counts2<-counts2[na.omit(match(res$ensembl_gene_id,rownames(counts2))),] #Remove genes with IDs from the counts table which cannot be annotated

##Create DGE object and specify counts object and grouping variable. Use one of the options below.
DGE<-DGEList(counts=counts,group=samples$Group) #Use for all datasets
DGE<-DGEList(counts=counts2, group=samples6$Group) #use for all datasets excluding GSE93326

##Filter low expression genes
keep <- filterByExpr(DGE) #ID low expression genes
DGE <- DGE[keep, , keep.lib.sizes=FALSE] # Remove low expression genes. The library size for each sample is recalculated to be the sum of the counts left in the rows after removing low expressed genes

##Calculate effective library sizes from raw library sizes using scaling factors
DGE <- calcNormFactors(DGE) #This gives TMM normalised counts and normalisation factors

#Extract normalised log2TMM data
cpm_data<-cpm(DGE,log=T) 

#Plot PCA normalized data 
PCA_Transposed<-t(cpm_data)

#Assign tissue type or study ID as row names to the PCA_Transposed matrix. Choose either depending on the clustering criteria.
rownames(PCA_Transposed)<- samples$Group #Cluster by tissue type
rownames(PCA_Transposed)<-samples$Study.ID #Cluster by study ID
rownames(PCA_Transposed)<- samples6$Group #Cluster by tissue type
rownames(PCA_Transposed)<- samples6$Study.ID #Cluster by study ID

#Convert the PCA_Transposed matrix into a data frame
PCA_Transposed<-as.data.frame(PCA_Transposed)

#Add tissue type or study ID as a column to the PCA_Transposed data frame to be able to assign them as a factor
Group<- samples[[3]] #Use when clustering per group
Study.ID<- samples[[4]] #Use when clustering per study ID
Group<-samples6[[3]] #Use when clustering per group
Study.ID<- samples6[[4]] #Use when clustering per study ID

PCA_Transposed<-as.data.frame(cbind(Group,PCA_Transposed)) #Use when clustering per group
PCA_Transposed<-as.data.frame(cbind(Study.ID,PCA_Transposed)) #Use when clustering per study ID

PCA_Transposed$`Group`<-as.factor(PCA_Transposed$`Group`) #Use when clustering per group
PCA_Transposed$`Study.ID`<-as.factor(PCA_Transposed$`Study.ID`) #Use when clustering per study ID

clusters<-as.factor(PCA_Transposed[,1])

#Remove the first column to have only numeric values for plotting
data<-PCA_Transposed[,-1]
data<-as.data.frame(apply(data,2,as.numeric))

#PCA 
pd2<-pca(data,scale = T, center = T, ncomp = 2) 

#Plot either per Group (tissue type) or study ID
#Generate PCA plot for all datasets
plotIndiv(pd2,ind.names = F,group = samples$Group, ellipse = T, star=F, size.xlabel = rel(1.5), size.ylabel = rel(1.5), size.axis = rel(1.5), legend = T, style="ggplot2", title="Clusters based on tissue type", size.title = rel(2), legend.position = "right", legend.title = "", size.legend = rel(1.5))
plotIndiv(pd2,ind.names = F,group = samples$Study.ID, ellipse = T, star=F, size.xlabel = rel(1.5), size.ylabel = rel(1.5), size.axis = rel(1.5), legend = T, style="ggplot2", title="Clusters based on study ID", size.title = rel(2), legend.position = "right", legend.title = "", size.legend = rel(1.5))

#Generate PCA plot for all datasets excluding GSE93326
plotIndiv(pd2,ind.names = F,group = samples6$Group, ellipse = T, star=F, size.xlabel = rel(1.8), size.ylabel = rel(1.8), size.axis = rel(1.8), legend = T, style="ggplot2", title="Clusters based on tissue type", size.title = rel(2.1), legend.position = "right", legend.title = "", size.legend = rel(2.2))
plotIndiv(pd2,ind.names = F,group = samples6$Study.ID, ellipse = T, star=F, size.xlabel = rel(1.8), size.ylabel = rel(1.8), size.axis = rel(1.8), legend = T, style="ggplot2", title="Clusters based on study ID", size.title = rel(2.1), legend.position = "right", legend.title = "", size.legend = rel(2))
