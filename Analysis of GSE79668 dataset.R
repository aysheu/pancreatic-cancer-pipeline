library(edgeR) #https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
library(biomaRt)
#library(mixOmics)

setwd("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE79668")

#Differential gene expression (DGE) analysis

# Read design table
samples<-read.table("./GSE79668_design.csv",sep = ",", header=T) 

# Create a design matrix for DGE analysis
design<-model.matrix(~samples$Group)

# Read and tidy up count data 
dat<-read.table("./GSE79668_STS vs LTS_Count data_EdgeR_NEW_duplicates removed.txt",sep=c("\t"),header=T)
counts<-dat[-1] #Remove the first column of dat which includes the Ensembl gene IDs
rownames(counts)<-dat[,1] #Assign the first column of dat as rownames to counts table

ENSG<-dat[,1] # Select list of gene IDs
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #Create a Mart object. This connects to biomart's ensembl database and the specified dataset
res<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'ensembl_gene_id', values = ENSG, mart = ensembl) #Retrieve gene names from the BioMart database
counts<-counts[na.omit(match(res$ensembl_gene_id,rownames(counts))),] #Remove genes with IDs which cannot be annotated

# Create DGE object and specify counts object and grouping variable
DGE<-DGEList(counts=counts,group=samples$Group)

# Filter low expression genes
keep <- filterByExpr(DGE) #ID low expression genes
DGE <- DGE[keep, , keep.lib.sizes=FALSE] #Remove low expression genes. The library size for each sample is recalculated to be the sum of the counts left in the rows after removing low expressed genes

# Calculate effective library sizes from raw library sizes using scaling factors
DGE <- calcNormFactors(DGE)

# Identify DGEs
DGE <- estimateDisp(DGE, design) #Estimate dispersions

et <- exactTest(DGE) ##Calculate differential expression. This finds genes dif. expressed in STS vs LTS subjects in the dataset
et$table$Gene_name<-res$external_gene_name[na.omit(match(rownames(et$table),res$ensembl_gene_id))] #Assign gene names

# Extract the list that are significant p< 0.05
p<-topTags(et,n=Inf,adjust.method="BH") #Extract the DEGs with p-value adjusted for multiple testing
keep<-p$table$FDR <=0.05 #Keep only the genes with FDR<0.05
p.table<-p[keep,]
write.table(p.table[keep,],"Genes with 5% FDR-adjusted p-values.csv",sep=",") #Save FDR-adjusted data
write.table(p, "GSE79668_Unadjusted and adjusted p-values for ALL genes.csv",sep=",") #Save all data including unadjusted p-values

# Look at the number of genes that are significantly up- or down-regulated and not differentially expressed
summary(decideTests(et))

plotMDS(DGE,labels = samples$Study.ID,col=as.numeric(samples$Group))

# Extract normalised log2TMM data 
cpm_data<-cpm(DGE,log=T) 
write.table(cpm_data, "GSE79668_tmm_data_updated.csv",sep=",")

# Scale each row (gene) to a mean of zero and a standard deviation one
logCPM <- t(scale(t(cpm_data)))
write.table(logCPM, "GSE79668_scaled_logTMM_data_updated.csv", sep=",")


# GENERATE A VIOLIN PLOT FOR MCTS IN STS vs LTS

# Subset the list p to extract the table of statistical data for all all genes
table2<-p[[1]] 

# Identify all MCTs that passed low-expression threshold which you want to create the violin plot for
MCTs<-dplyr::filter(table2, grepl("SLC16", Gene_name))

# Select MCTs that were differentially expressed in stroma samples from Maurer et al. 
# Select only top 4 MCTs in terms of their FDR levels: SLC16A10, SLC16A7, SLC16A5, and SLC16A1-AS1
library(dplyr)
MCT_selected<- filter(MCTs, Gene_name %in% c('SLC16A10', 'SLC16A7', 'SLC16A5', 'SLC16A1-AS1')) 

# Subset the cpm_data table to capture only MCTs of interest
MCT_TMM<-cpm_data[match(rownames(MCT_selected),rownames(cpm_data)),]

# Convert data presentation from the current wild format to long format
# For conversion, first assign the rownames of MCT_TMM table (Ensembl IDs) into a column
# Convert MCT_TMM from matrix to data frame so that Ensembl IDs are added to the last column
MCT_TMM<-as.data.frame(MCT_TMM)
MCT_TMM['EnsemblID'] <-NA
MCT_TMM$EnsemblID=rownames(MCT_TMM) #Make sure that when this column is added, all other columns are still numeric rather than being converted to character!

#Next assign the Ensembl ID column as a factor
MCT_TMM$EnsemblID<-factor(MCT_TMM$EnsemblID)
MCT_TMM2$EnsemblID<-factor(MCT_TMM2$EnsemblID)

# Conversion to long format
library(tidyr)
data_long<-gather(MCT_TMM, sampleID, Count, SRR3308916:SRR3308930, factor_key=T) # ref: http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/

# Convert Ensembl IDs to MCT names
EnsemblID<-data_long[,1] # Select list of gene IDs
res2<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'ensembl_gene_id', values = EnsemblID, mart = ensembl) #Retrieve gene names from the BioMart database

# Add the gene IDs as a new column to data_long
data_long['GeneID'] <-NA

# Rename colname of Res2 to EnsembleID  for tidyness 
library(plyr)
res2<-plyr::rename(res2, c("ensembl_gene_id"="EnsemblID"))

# Match the EnsembleID on both dataframes and put result into GeneID of data_long
data_long$GeneID <- res2$external_gene_name[match(data_long$EnsemblID, res2$EnsemblID)]

# Add the survival status as a new column to data_long
data_long['SurvivalStatus'] <- NA

# Match the EnsembleID on both dataframes and put survival status into data_long
data_long$SurvivalStatus <- samples$Group[match(data_long$sampleID, samples$SampleID)]

# violin plot
# reference http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
library(ggplot2) 
p3 <- ggplot(data_long, aes(x=GeneID, y=Count, fill=SurvivalStatus) , trim = FALSE , position = position_dodge(0.9))

p3 + geom_violin()  + stat_summary(fun=median, geom="point", shape=23, size=2 , color="red") + 
  geom_boxplot(width=0.1, position = position_dodge(0.9)) + 
  scale_x_discrete(limits=c("SLC16A5", "SLC16A7", "SLC16A10", "SLC16A1-AS1")) +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(face = "italic"), legend.position="top", legend.key.size = unit(1, 'cm'), legend.title = element_text(size = 12), legend.text = element_text(size=12), axis.text = element_text(colour = "black", size = rel(1.0)), axis.title = element_text(size = rel(1), face = "bold"), plot.background = element_rect(colour = "white")) +
  annotate(geom="text", x=1, y=-1, label="0.784", color="black", size=4) +
  annotate(geom="text", x=2, y=-1, label="0.003", color="black", size=4) +
  annotate(geom="text", x=3, y=-1, label="0.402", color="black", size=4) +
  annotate(geom="text", x=4, y=-1, label="0.006", color="black", size=4) +
  annotate(geom="text", x=2.5, y=-1.8, label="p-Values", color="black", size=4) +
  scale_y_continuous(limits = c(-2,8)) +
  labs(title="")


# GENERATE A HEATMAP TO SHOW EXPRESSION OF DE-ECMs IN STROMA (N=502) IN STS and LTS SUBJECTS
source("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/heatmap3.R")

# Read the list of DE-ECM-related genes (n=502) in stroma-epithelium
All_ECM_genes<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS/All ECM-related DEGs in stroma.csv", sep=",", header=T)
All_ECM_genes<-All_ECM_genes2[,-1]
rownames(All_ECM_genes)<- All_ECM_genes[,1]

# Match gene names from All_ECM_genes dataset with gene names in table2 (genes that passed low-expression threshold in STS-LTS dataset)
table3 <- table2[match(rownames(All_ECM_genes), rownames(table2)),] # Not all of these ECMs are differentially expressed in STS and LTS subjects so we need to subset from this list those ECMs that are significant in STS and LTS subjects
table4 <- table3[rownames(table3)[table3$FDR<0.05],]
table4 <- na.omit(table4) # See the list of ECMs that are differentially expressed both in stroma-epithelium and STS-LTS subjects
write.csv(table4, "Differentially expressed matrisome genes in STS relative to LTS.csv")

# Produce a heatmap using log-normalised data.
#selcpm <- cpm_data[match(rownames(table3), rownames(cpm_data)),]
#selcpm <- na.omit(selcpm) # remove NAs: these are the ECMs that are differentially expressed in stroma but have been filtered out from Kirby dataset due to low expression levels
#rownames(selcpm)<-table3$Gene_name[match(rownames(selcpm),rownames(table3))]

#look at ECMs with p<0.05
#selcpm_p<-cpm_data[match(rownames(table3)[table3$PValue<0.05], rownames(cpm_data)),]
#selcpm_p<-na.omit(selcpm_p)
#rownames(selcpm_p)<-table3$Gene_name[match(rownames(selcpm_p),rownames(table3))]

# Generate a heatmap using ECMs with FDR<0.05
selcpm_FDR <- cpm_data[match(rownames(table3)[table3$FDR<0.05], rownames(cpm_data)),]
selcpm_FDR <- na.omit(selcpm_FDR) # remove NAs: these are the ECMs that are differentially expressed in stroma but have been filtered out from Kirby dataset due to low expression levels
rownames(selcpm_FDR)<-table3$Gene_name[match(rownames(selcpm_FDR),rownames(table3))]

# Generate a heatmap object
hm<-as.matrix(selcpm_FDR)

## Assign colours for STS and LTS. Ref: https://www.biostars.org/p/18211/
library(plyr)
tissue_colour=revalue(samples$Group, c("STS"="magenta", "LTS"="seagreen")) ## ref for colors: http://applied-r.com/r-color-tables/
clab=tissue_colour
clab<-as.data.frame(clab)
colnames(clab)<-"STS vs LTS"
clab<-as.matrix(clab)
# sample_type<- as.factor(samples$Group)
# dendog = cluster_within_group(hm, sample_type)
order.dendrogram(dendog)

# As we are using log-normalised data which has not been z-scored, we are normalising the data by row below
library(gplots)

#Use this when you are not separating STS and LTS into 2 groups. This way, dendrogram can be used for columns.
heatmap.3(hm, col=bluered, ColSideColors = clab, ColSideColorsSize=2, Rowv=TRUE, key =TRUE, scale="row",trace="none", dendrogram="both", cexRow=1, cexCol=1, density.info="none", margin=c(10,15)) 
legend("topright",legend=c("STS","LTS"),
       fill=c("magenta","seagreen"), border=FALSE, bty="n", y.intersp = 1, cex=1)

# Use this when grouping STS separately than LTS.
#ref for grouping of columns/samples: https://biocorecrg.github.io/CRG_RIntroduction/heatmap-2-function-from-gplots-package.html
heatmap.3(hm, col=bluered, ColSideColors = clab, ColSideColorsSize=2, Rowv=TRUE, key =TRUE, scale="row",trace="none", dendrogram="row", cexRow=1, cexCol=1, density.info="none", margin=c(10,15), Colv = FALSE) 
legend("topright",legend=c("STS","LTS"),
       fill=c("magenta","seagreen"), border=FALSE, bty="n", y.intersp = 1, cex=1)
