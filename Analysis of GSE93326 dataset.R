library(edgeR) #https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
library(biomaRt)
library(mixOmics)
library(plyr)
library(gplots)

setwd("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326")

#Differential gene expression (DGE) analysis

# Read design table for matched epithelial and stromal data with additional stromal data
samples<-read.table("./GSE93326_combined 57 stroma_design.csv",sep = ",", header=T) #to be prepared

#Rename AJCC stages - convert normal numbers to Roman numbers. This will be used in PCA
samples$AJCC.stage<-gsub("2B","IIB",samples$AJCC.stage) 
samples$AJCC.stage<-gsub("2A","IIA",samples$AJCC.stage) 
samples$AJCC.stage<-gsub("3","III",samples$AJCC.stage) 
samples$AJCC.stage<-gsub("4","IV",samples$AJCC.stage) 
samples$AJCC.stage<-gsub("1A","IA",samples$AJCC.stage) 

# Create a design matrix for DGE analysis
design<-model.matrix(~samples$Group)

# Read count data for marched data with additional stroma
dat<-read.table("./GSE93326,Count_S vs E_combined with additional 57 stroma_EdgeR.txt",sep=c("\t"),header=T)

# Read the count table containing all samples and genes and clean data to remove NAs
counts<-dat[-1] # Remove the first column of dat which includes the Ensembl gene IDs
counts<-dat[-c(60609:60613), -1] # Omit the last 4 rows of the count table as counts are NA
rownames(counts)<-dat[-c(60609:60613),1] #Assign the first column of dat as rownames to counts table

ENSG<-dat[,1] # select list of gene IDs
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #Create a Mart object. This connects to biomart's ensembl database and the specified dataset
res<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'ensembl_gene_id', values = ENSG, mart = ensembl) #Retrieve gene names from the BioMart database
counts<-counts[na.omit(match(res$ensembl_gene_id,rownames(counts))),] #Remove genes with IDs which cannot be annotated

#Create DGE object and specify counts object and grouping variable
DGE<-DGEList(counts=counts,group=samples$Group) #Specify counts object and grouping variable

#Filter low expression genes
keep <- filterByExpr(DGE) #ID low expression genes
DGE <- DGE[keep, , keep.lib.sizes=FALSE] # Remove low expression genes. The library size for each sample is recalculated to be the sum of the counts left in the rows after removing low expressed genes

#Calculate effective library sizes from raw library sizes using scaling factors
DGE <- calcNormFactors(DGE) #This gives TMM normalised counts and normalisation factors

#Identify DGEs
DGE <- estimateDisp(DGE, design) # Estimate dispersions

et <- exactTest(DGE) #Calculate differential expression. This finds genes dif. expressed in stroma vs epithelium in the dataset
et$table$Gene_name<-res$external_gene_name[na.omit(match(rownames(et$table),res$ensembl_gene_id))] # assign gene names

#Extract the list that are significant p< 0.05
p<-topTags(et,n=Inf,adjust.method="BH") #Extract the DEGs with p-value adjusted for multiple testing
keep<-p$table$FDR <=0.05 #Keep only the genes with FDR<0.05
p.table<-p[keep,]

#Save FDR-adjusted data
write.table(p.table[keep,],"Genes with 5% FDR-adjusted p-values_matched plus additional stroma.csv",sep=",")

# Save all data
write.table(p, "Unadjusted and adjusted p-values for ALL genes_matched plus additional stroma.csv",sep=",")

#Look at the number of genes that are significantly up- or down-regulated and not differentially expressed
summary(decideTests(et))

#Extract normalised log2TMM data 
cpm_data<-cpm(DGE,log=T) 
write.table(cpm_data, "GSE93326_matched plus additional stroma_logTMM_data.csv",sep=",")

#Scale each row (gene) to a mean of zero and a standard deviation one
logCPM <- t(scale(t(cpm_data)))
write.table(logCPM, "GSE93326_matched plus additional stroma_scaledlogTMM_data.csv", sep=",")


# Principal Component Analysis

# Read count data for marched data with additional stroma
dat<-read.table("./GSE93326,Count_S vs E_combined with additional 57 stroma_EdgeR.txt",sep=c("\t"),header=T)

#Use only when clustering based on tumour overall stage, grade or N.score
library(dplyr)
library(tidyr)
samples2<-samples[-c(63,129),] #First remove samples without metadata from the original data frame. Use when clustering based on tumour stage and grade 
samples3<- samples2 %>% drop_na(N.score) # Then drop samples from the list without N.score info. Only to be used when clustering based on N-score

# Read the count table containing all samples and genes and clean data to remove NAs
counts<-dat[-1] # Remove the first column of dat which includes the Ensembl gene IDs
counts2<-dat[-c(60609:60613), -c(1,63,129)] # Omit the last 4 rows of the count table as counts are NA and the samples without any metadata related to tumour stage and grade
counts3<-dat[-c(60609:60613), -c(1,63,129,12,48,52,78,114,118,137,143,149,151)] # Omit the last 4 rows of the count table as counts are NA and the samples without tumour N.score info
rownames(counts2)<-dat[-c(60609:60613),1] #Assign the first column of dat as rownames to counts table
rownames(counts3)<-dat[-c(60609:60613),1] #Assign the first column of dat as rownames to counts table

ENSG<-dat[,1] # select list of gene IDs
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #Create a Mart object. This connects to biomart's ensembl database and the specified dataset
res<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'ensembl_gene_id', values = ENSG, mart = ensembl) #Retrieve gene names from the BioMart database
counts<-counts[na.omit(match(res$ensembl_gene_id,rownames(counts))),] #Remove genes with IDs which cannot be annotated. Use when conducting PCA for all samples regardless of tumour stage, grade, or lymph node metastasis
counts2<-counts2[na.omit(match(res$ensembl_gene_id,rownames(counts2))),] #Remove genes with IDs which cannot be annotated. Use when conducting PCA for samples with tumour stage and grade info
counts3<-counts3[na.omit(match(res$ensembl_gene_id,rownames(counts3))),] #Remove genes with IDs which cannot be annotated. Use when conducting PCA for samples with tumour lymph node metastasis (N.score) info

#Create DGE object and specify counts object and grouping variable. Choose from below options depending on the type of PCA plot that will be generated
DGE<-DGEList(counts=counts,group=samples$Group) #Choose when conducting DGE analysis and PCA for all samples
DGE<-DGEList(counts=counts2, group=samples2$Group) #Choose when conducting PCA based on available tumour stage and grade info
DGE<-DGEList(counts=counts3, group=samples3$Group) #Choose when conducting PCA based on available tumour N.score info

#Filter low expression genes
keep <- filterByExpr(DGE) #ID low expression genes
DGE <- DGE[keep, , keep.lib.sizes=FALSE] # Remove low expression genes. The library size for each sample is recalculated to be the sum of the counts left in the rows after removing low expressed genes

#Calculate effective library sizes from raw library sizes using scaling factors
DGE <- calcNormFactors(DGE) #This gives TMM normalised counts and normalisation factors

#Extract normalised log2TMM data 
cpm_data<-cpm(DGE,log=T) 

#Take transpose of the cpm data
PCA_Transposed<-t(cpm_data)

#Use either of below to assign either group, N.score, AJCC.stage or grade as rownames. Use only one at a time.
rownames(PCA_Transposed)<- samples$Group #Use when clustering per group
rownames(PCA_Transposed)<- samples2$AJCC.stage #Use when clustering per AJCC stage
rownames(PCA_Transposed)<- samples2$Grade #Use when clustering per grade
rownames(PCA_Transposed)<- samples3$N.score #Use when clustering N.score

#Convert the matrix into dataframe
PCA_Transposed<-as.data.frame(PCA_Transposed)

#Identify different factors to be combined with the PCA_Transposed dataframe. Use only one at a time.
Group<- samples[[2]] ## Use this line to cluster per group
AJCC.stage<-samples2[[3]] ## Cluster per AJCC.stage
Grade<-samples2[[7]] ## Cluster per grade
N.score<-samples3[[5]] ## Cluster per N.score

#Combine the factors identified above with the PCA_Transposed dataframe.Use only one at a time
PCA_Transposed<-as.data.frame(cbind(Group,PCA_Transposed))
PCA_Transposed<-as.data.frame(cbind(AJCC.stage,PCA_Transposed))
PCA_Transposed<-as.data.frame(cbind(Grade,PCA_Transposed))
PCA_Transposed<-as.data.frame(cbind(N.score,PCA_Transposed))

#Assign factor columns as a factor. Use only one at a time
PCA_Transposed$`Group`<-as.factor(PCA_Transposed$`Group`)
PCA_Transposed$`AJCC.stage`<-as.factor(PCA_Transposed$`AJCC.stage`)
PCA_Transposed$`Grade`<-as.factor(PCA_Transposed$`Grade`)
PCA_Transposed$`N.score`<-as.factor(PCA_Transposed$`N.score`)

#Identify clusters as a factor
clusters<-as.factor(PCA_Transposed[,1])

#Remove the first column from the PCA_Transposed dataframe to have only numeric data
data<-PCA_Transposed[,-1]
data<-as.data.frame(apply(data,2,as.numeric))#make values numeric

#PCA 
pd2<-pca(data,scale = T, center = T, ncomp = 2)
plotIndiv(pd2,ind.names = F,group = samples$Group, ellipse = T, star=F, size.xlabel = rel(1.5), size.ylabel = rel(1.5), size.axis = rel(1.5), legend = T, style="ggplot2", title="Tumour epithelium & stroma", legend.position = "right", legend.title = "", size.legend = rel(1.5))
plotIndiv(pd2,ind.names =F, group = samples2$AJCC.stage, ellipse = T, star=F, size.xlabel = rel(1.5), size.ylabel = rel(1.5), size.axis = rel(1.5), legend=T, style="ggplot2", title="Epithelium & stroma distribution based on tumour stage", legend.position = "right", legend.title = "Tumour stage", size.legend.title = rel(1.5), size.legend = rel(1.5))
plotIndiv(pd2,ind.names =F, group = samples2$Grade, ellipse = T, star=F, size.xlabel = rel(1.5), size.ylabel = rel(1.5), size.axis = rel(1.5), legend=T, style="ggplot2", title="Epithelium & stroma distribution based on tumour grade",legend.position = "right", legend.title = "Tumour grade", size.legend.title = rel(1.5), size.legend = rel(1.5))
plotIndiv(pd2,ind.names =F, group = samples3$N.score, ellipse = T, star=F, size.xlabel = rel(1.5), size.ylabel = rel(1.5), size.axis = rel(1.5),legend=T, style="ggplot2", title="Epithelium & stroma distribution based on N-score", legend.position = "right", legend.title = "N-score", size.legend.title = rel(1.5), size.legend = rel(1.5))


# GENERATE A VOLCANO PLOT FOR ALL GENES AND ASSOCIATED P-VALUES

# Read design table for matched epithelial and stromal data with additional stromal data
samples<-read.table("./GSE93326_combined 57 stroma_design.csv",sep = ",", header=T) #to be prepared

# Read the count table containing all samples and genes and clean data to remove NAs
counts<-dat[-1] # Remove the first column of dat which includes the Ensembl gene IDs
counts<-dat[-c(60609:60613), -1] # Omit the last 4 rows of the count table as counts are NA
rownames(counts)<-dat[-c(60609:60613),1] #Assign the first column of dat as rownames to counts table

counts<-counts[na.omit(match(res$ensembl_gene_id,rownames(counts))),] #Remove genes with IDs which cannot be annotated

#Create DGE object and specify counts object and grouping variable
DGE<-DGEList(counts=counts,group=samples$Group) #Specify counts object and grouping variable

#Filter low expression genes
keep <- filterByExpr(DGE) #ID low expression genes
DGE <- DGE[keep, , keep.lib.sizes=FALSE] # Remove low expression genes. The library size for each sample is recalculated to be the sum of the counts left in the rows after removing low expressed genes

#Calculate effective library sizes from raw library sizes using scaling factors
DGE <- calcNormFactors(DGE) #This gives TMM normalised counts and normalisation factors

#Identify DGEs
DGE <- estimateDisp(DGE, design) # Estimate dispersions

et <- exactTest(DGE) #Calculate differential expression. This finds genes dif. expressed in stroma vs epithelium in the dataset
et$table$Gene_name<-res$external_gene_name[na.omit(match(rownames(et$table),res$ensembl_gene_id))] # assign gene names

#Extract the list that are significant p< 0.05
p<-topTags(et,n=Inf,adjust.method="BH") #Extract the DEGs with p-value adjusted for multiple testing
keep<-p$table$FDR <=0.05 #Keep only the genes with FDR<0.05
p.table<-p[keep,]

# Subset the list p to extract the table of statistical data for all all genes
table2<-p[[1]] 
logFC<-table2[[1]]
FDR<-table2[[5]]

#Create a data frame from the vector elements for the volcano plot
DF<-data.frame(x=logFC, y=-log10(FDR))

#Rename the column names
names(DF)[1] <- "logFC"
names(DF)[2] <- "-log10(FDR)"
ggplot(data=DF, aes(x=logFC, y=-log10(FDR))) + geom_point(aes(colour=FDR <0.05))


## PRODUCE A HEATMAP FOR STROMAL AND EPITHELIAL SAMPLES
source("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/heatmap3.R")

# Read design table for matched epithelial and stromal data with additional stromal data
samples<-read.table("./GSE93326_combined 57 stroma_design.csv",sep = ",", header=T) #to be prepared

# Read the count table containing all samples and genes and clean data to remove NAs
counts<-dat[-1] # Remove the first column of dat which includes the Ensembl gene IDs
counts<-dat[-c(60609:60613), -1] # Omit the last 4 rows of the count table as counts are NA
rownames(counts)<-dat[-c(60609:60613),1] #Assign the first column of dat as rownames to counts table

counts<-counts[na.omit(match(res$ensembl_gene_id,rownames(counts))),] #Remove genes with IDs which cannot be annotated

#Create DGE object and specify counts object and grouping variable
DGE<-DGEList(counts=counts,group=samples$Group) #Specify counts object and grouping variable

#Filter low expression genes
keep <- filterByExpr(DGE) #ID low expression genes
DGE <- DGE[keep, , keep.lib.sizes=FALSE] # Remove low expression genes. The library size for each sample is recalculated to be the sum of the counts left in the rows after removing low expressed genes

#Calculate effective library sizes from raw library sizes using scaling factors
DGE <- calcNormFactors(DGE) #This gives TMM normalised counts and normalisation factors

#Identify DGEs
DGE <- estimateDisp(DGE, design) # Estimate dispersions

et <- exactTest(DGE) #Calculate differential expression. This finds genes dif. expressed in stroma vs epithelium in the dataset
et$table$Gene_name<-res$external_gene_name[na.omit(match(rownames(et$table),res$ensembl_gene_id))] # assign gene names

#Extract the list that are significant p< 0.05
p<-topTags(et,n=Inf,adjust.method="BH") #Extract the DEGs with p-value adjusted for multiple testing
keep<-p$table$FDR <=0.05 #Keep only the genes with FDR<0.05
p.table<-p[keep,]

#Extract normalised log2TMM data 
cpm_data<-cpm(DGE,log=T) 

#Scale each row (gene) to a mean of zero and a standard deviation one
logCPM <- t(scale(t(cpm_data)))

# Select 5% of DE-genes and produce heatmap using log-normalised and z-scored (scaled) data.
selcpm <- logCPM[rownames(p$table)[p$table$FDR<0.05],]

#Generate heat map
hm<-as.matrix(selcpm)

## Assign colours for tissue type and study ID. Ref: https://www.biostars.org/p/18211/
tissue_colour=revalue(samples$Group, c("Stroma"="red3", "Epithelium"="seagreen")) ## ref for colors: http://applied-r.com/r-color-tables/
clab=tissue_colour
clab<-as.data.frame(clab)
colnames(clab)<-"Stroma vs Epithelium"
clab<-as.matrix(clab)

# As we are using log-normalised data which has been z-scored, we are not normalising the data by row below
# Ref for removing row and col names: https://stackoverflow.com/questions/37241205/remove-row-and-column-names-heatmap-2
# The legend and axes labels have been edited manually which required manual adjustment of the margins
heatmap.3(hm, breaks=seq(-2,2,0.5), col=bluered, ColSideColors = clab, ColSideColorsSize=2, Rowv=TRUE, key =TRUE, scale="none",trace="none", dendrogram="both", cexRow=1, cexCol=1, density.info="none", margin=c(15,12), xlab = "Samples", ylab = "Genes", labRow = FALSE, labCol = FALSE) 
legend("topright",legend=c("Stroma","Epithelium"),
       fill=c("red3","seagreen"), border=FALSE, bty="n", y.intersp = 1, cex=1) 


# GENERATE VIOLIN PLOT FOR MCTS IN STROMA VS EPITHELIUM

# Read design table for matched epithelial and stromal data with additional stromal data
samples<-read.table("./GSE93326_combined 57 stroma_design.csv",sep = ",", header=T) #to be prepared

# Read the count table containing all samples and genes and clean data to remove NAs
dat<-read.table("./GSE93326,Count_S vs E_combined with additional 57 stroma_EdgeR.txt",sep=c("\t"),header=T)
counts<-dat[-1] # Remove the first column of dat which includes the Ensembl gene IDs
counts<-dat[-c(60609:60613), -1] # Omit the last 4 rows of the count table as counts are NA
rownames(counts)<-dat[-c(60609:60613),1] #Assign the first column of dat as rownames to counts table

ENSG<-dat[,1] # select list of gene IDs
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #Create a Mart object. This connects to biomart's ensembl database and the specified dataset
res<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'ensembl_gene_id', values = ENSG, mart = ensembl) #Retrieve gene names from the BioMart database
counts<-counts[na.omit(match(res$ensembl_gene_id,rownames(counts))),] #Remove genes with IDs which cannot be annotated

# Create a design matrix for DGE analysis
design<-model.matrix(~samples$Group)

#Create DGE object and specify counts object and grouping variable
DGE<-DGEList(counts=counts,group=samples$Group) #Specify counts object and grouping variable

#Filter low expression genes
keep <- filterByExpr(DGE) #ID low expression genes
DGE <- DGE[keep, , keep.lib.sizes=FALSE] # Remove low expression genes. The library size for each sample is recalculated to be the sum of the counts left in the rows after removing low expressed genes

#Calculate effective library sizes from raw library sizes using scaling factors
DGE <- calcNormFactors(DGE) #This gives TMM normalised counts and normalisation factors

#Identify DGEs
DGE <- estimateDisp(DGE, design) # Estimate dispersions

et <- exactTest(DGE) #Calculate differential expression. This finds genes dif. expressed in stroma vs epithelium in the dataset
et$table$Gene_name<-res$external_gene_name[na.omit(match(rownames(et$table),res$ensembl_gene_id))] # assign gene names

#Extract normalised log2TMM data 
cpm_data<-cpm(DGE,log=T)

#Extract the list that are significant p< 0.05
p<-topTags(et,n=Inf,adjust.method="BH") #Extract the DEGs with p-value adjusted for multiple testing
keep<-p$table$FDR <=0.05 #Keep only the genes with FDR<0.05
p.table<-p[keep,]

# Subset the list p to extract the table of statistical data for all all genes
table2<-p[[1]] 

# Identify all MCTs that passed low-expression threshold which you want to create the violin plot for
MCTs<-dplyr::filter(table2, grepl("SLC16", Gene_name))

# Subset the cpm_data table to capture only MCTs
MCT_TMM<-cpm_data[match(rownames(MCTs),rownames(cpm_data)),]

# Convert data presentation from the current wild format to long format
# For conversion, first assign the rownames of MCT_TMM table(Ensembl IDs) into a column
# Convert MCT_TMM from matrix to data frame so that Ensembl IDs are added to the last column
MCT_TMM<-as.data.frame(MCT_TMM)
MCT_TMM['EnsemblID'] <-NA # Add a new column for Ensembl IDs
MCT_TMM$EnsemblID=rownames(MCT_TMM) #Make sure that when this column is added, all other columns are still numeric rather than being converted to character!

#Next, assign the Ensembl ID column as a factor
MCT_TMM$EnsemblID<-factor(MCT_TMM$EnsemblID)

library(tidyr)

# Conversion to long format
data_long<-gather(MCT_TMM, sampleID, Count, SRR5163409.10:SRR7963431, factor_key=T) # ref: http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/
write.csv(data_long, "Data_long.csv",sep=",")

#Convert Ensembl IDs to MCT names
EnsemblID<-data_long[,1] # select list of gene IDs
res2<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'ensembl_gene_id', values = EnsemblID, mart = ensembl) #Retrieve gene names from the BioMart database
write.table(res2, "Res2_ensemblid_genename.csv",sep=",")

# Add a new column to the data_long data frame
data_long['GeneID'] <-NA

# Rename Res2 colname to EnsembleID  for tidyness not required
library(plyr)
res2<-plyr::rename(res2, c("ensembl_gene_id"="EnsemblID")) # ref: https://statisticsglobe.com/r-error-cant-rename-columns-that-dont-exist

# Match the EnsembleID on both dataframes and put result into Gene ID of data_long
data_long$GeneID <- res2$external_gene_name[match(data_long$EnsemblID, res2$EnsemblID)]

# Add new column to the data_long data frame
data_long['TissueType'] <-NA

# In samples data frame, replace "/" with "." in the SampleID column
samples$SampleID <- gsub('/', '.', samples$SampleID)

# Match the EnsembleID on both data frames and put result into Gene ID of data_long
data_long$TissueType <- samples$Group[match(data_long$sampleID, samples$SampleID)]

# Add a new column to data_long
data_long['Age'] <- NA

# Match sample IDs on data_long and samples dataframes and add age info to data_long
data_long$Age <- samples$Age[match(data_long$sampleID, samples$SampleID)]

# Generate the violin plot
# reference http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
library(ggplot2) 

tiff("Violin plot of MCTs in stroma-epithelium_edited for paper_10Feb22.tiff", units="in", width=12, height=6, res=300)

p3 <- ggplot(data_long, aes(x=GeneID, y=Count, fill=TissueType) , trim = FALSE , position = position_dodge(0.9))

p3 + geom_violin()  + stat_summary(fun=median, geom="point", shape=23, size=2 , color="red") + 
  geom_boxplot(width=0.1, position = position_dodge(0.9)) + 
  scale_x_discrete(limits=c("SLC16A1", "SLC16A2","SLC16A3", "SLC16A4","SLC16A5", "SLC16A6", "SLC16A7", "SLC16A9", "SLC16A10", "SLC16A13", "SLC16A14", "SLC16A1-AS1")) +
  theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(face = "italic"), legend.position="top", legend.key.size = unit(1, 'cm'), legend.title = element_text(size = 12), legend.text = element_text(size=12), axis.text = element_text(colour = "black", size = rel(1.0)), axis.title = element_text(size = rel(1), face = "bold"), plot.background = element_rect(colour = "white")) +
  annotate(geom="text", x=1, y=-3.4, label="1.61e-1", color="black") +
  annotate(geom="text", x=2, y=-3.4, label="2.25e-11", color="black") +
  annotate(geom="text", x=3, y=-3.4, label="2.05e-1", color="black") +
  annotate(geom="text", x=4, y=-3.4, label="7.97e-9", color="black") +
  annotate(geom="text", x=5, y=-3.4, label="2.50e-12", color="black") +
  annotate(geom="text", x=6, y=-3.4, label="6.54e-5", color="black") +
  annotate(geom="text", x=7, y=-3.4, label="2.02e-16", color="black") +
  annotate(geom="text", x=8, y=-3.4, label="1.03e-3", color="black") +
  annotate(geom="text", x=9, y=-3.4, label="5.95e-24", color="black") +
  annotate(geom="text", x=10, y=-3.4, label="2.43e-3", color="black") +
  annotate(geom="text", x=11, y=-3.4, label="3.89e-1", color="black") +
  annotate(geom="text", x=12, y=-3.4, label="1.77e-11", color="black") +
  annotate(geom="text", x=7, y=-4.0, label="FDR Values", color="black") +
  scale_y_continuous(limits = c(-4,10)) +
  labs(title="")
dev.off()

# PEARSON CORRELATION ANALYSIS

# For the correlation analysis, need to correlate age at diagnosis with MCTs and ECMs.

# Read design table for matched epithelial and stromal data with additional stromal data
samples<-read.table("./GSE93326_combined 57 stroma_design.csv",sep = ",", header=T) #to be prepared

# Read the count table containing all samples and genes and clean data to remove NAs
counts<-dat[-1] # Remove the first column of dat which includes the Ensembl gene IDs
counts<-dat[-c(60609:60613), -1] # Omit the last 4 rows of the count table as counts are NA
rownames(counts)<-dat[-c(60609:60613),1] #Assign the first column of dat as rownames to counts table

counts<-counts[na.omit(match(res$ensembl_gene_id,rownames(counts))),] #Remove genes with IDs which cannot be annotated

#Create DGE object and specify counts object and grouping variable
DGE<-DGEList(counts=counts,group=samples$Group) #Specify counts object and grouping variable

#Filter low expression genes
keep <- filterByExpr(DGE) #ID low expression genes
DGE <- DGE[keep, , keep.lib.sizes=FALSE] # Remove low expression genes. The library size for each sample is recalculated to be the sum of the counts left in the rows after removing low expressed genes

#Calculate effective library sizes from raw library sizes using scaling factors
DGE <- calcNormFactors(DGE) #This gives TMM normalised counts and normalisation factors

#Identify DGEs
DGE <- estimateDisp(DGE, design) # Estimate dispersions

et <- exactTest(DGE) #Calculate differential expression. This finds genes dif. expressed in stroma vs epithelium in the dataset
et$table$Gene_name<-res$external_gene_name[na.omit(match(rownames(et$table),res$ensembl_gene_id))] # assign gene names

#Extract the list that are significant p< 0.05
p<-topTags(et,n=Inf,adjust.method="BH") #Extract the DEGs with p-value adjusted for multiple testing
keep<-p$table$FDR <=0.05 #Keep only the genes with FDR<0.05
p.table<-p[keep,]

# Subset the list p to extract the table of statistical data for all all genes
table2<-p[[1]] 

# Identify all MCTs that passed low-expression threshold which you want to create the violin plot for
MCTs<-dplyr::filter(table2, grepl("SLC16", Gene_name))

# First subset the cpm_data table to capture only MCTs
MCT_TMM<-cpm_data[match(rownames(MCTs),rownames(cpm_data)),]

# Then transpose MCT_TMM data frame
MCT_TMM_transposed<- t(MCT_TMM)
MCT_TMM_transposed <- as.data.frame(MCT_TMM_transposed)

# Need to add age of each patient to the MCT_TMM_transposed data frame
# Add two new columns and label as age and tissue type
MCT_TMM_transposed['Age'] <-1

# In samples data frame, replace "/" with "." in the SampleID column 
samples$SampleID <- gsub('/', '.', samples$SampleID)

# Add age and tissue type data by matching sample IDs from samples dataset
MCT_TMM_transposed$Age <- samples$Age[match(rownames(MCT_TMM_transposed),samples$SampleID)]

# Remove rows with NAs (no age info) from the dataset
MCT_TMM_transposed <- na.omit(MCT_TMM_transposed) 
write.csv(MCT_TMM_transposed, "Maurer_Age vs MCT expression_logTMM.csv")

# Identify the ECMs that showed the strongest association with MCTs in hypernetwork analysis
ECMs<- read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS/Analysis of MCTs using r-value threshold of 0.2 for binarisation/Non-DE Matrisome genes correlated with ALL MCTs.csv", sep=",", header=TRUE)
ECMs<-ECMs[,-1]

# Subset the cpm_data table to capture only the ECMs that showed the strongest correlation with MCTs
ECM_TMM<-cpm_data[match(ECMs$Gene.stable.ID,rownames(cpm_data)),]

# Rename rownames of ECM_TMM matrix
# First read the ECMs/matrisome related genes that strongly correlate with MCTs
ECM_cor_MCTs<-read.csv("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE93326/HYPERNETWORK ANALYSIS/Analysis of MCTs using r-value threshold of 0.2 for binarisation/Non-DE Matrisome genes correlated with ALL MCTs.csv", sep=",", header=T)
ECM_cor_MCTs<-ECM_cor_MCTs[,-1]
rownames(ECM_TMM)<- ECM_cor_MCTs$Gene.name[match(rownames(ECM_TMM), ECM_cor_MCTs$Gene.stable.ID)]

#Transpose the ECM_TMM dataset
ECM_TMM_transposed<- t(ECM_TMM)
ECM_TMM_transposed<- as.data.frame(ECM_TMM_transposed)

# Add two new columns and label as age and tissue type
ECM_TMM_transposed['Age'] <-NA

# In samples data frame, replace "/" with "." in the SampleID column 
samples$SampleID <- gsub('/', '.', samples$SampleID)

# Add age and tissue type data by matching sample IDs from samples dataset
ECM_TMM_transposed$Age <- samples$Age[match(rownames(ECM_TMM_transposed),samples$SampleID)]

# Remove rows with NAs (no age info) from the dataset
ECM_TMM_transposed <- na.omit(ECM_TMM_transposed)
#ECM_TMM_transposed <- as.matrix(ECM_TMM_transposed) # TO BE USED FOR CORRELATION ANALYSIS
write.csv(ECM_TMM_transposed, "Maurer_Age vs ECM expression_logTMM.csv") 
write.csv(ECM_TMM_transposed, "Maurer_Age vs ECM expression_logTMM_version2.csv")


# Conduct correlation analysis between MCTs/ECMs and Age at diagnosis
# Extract correlation p-values
library(ggpubr)
library(Hmisc)


# Produce correlation r-values for the MCT_TMM_transposed dataframe
corAll<- cor(MCT_TMM_transposed, method = "pearson")

# Correlation r- and p-values can also be produced using "rcorr" function 
# This requires conversion of dataframe into matrix
MCT_TMM_transposed2<-as.matrix(MCT_TMM_transposed)

# Convert Ensembl IDs to gene names
#MCTlist<-read.csv("./MCTlist.csv", sep=",", header=T)
colnames(MCT_TMM_transposed2) <- c("SLC16A10", "SLC16A7", "SLC16A5", "SLC16A1-AS1", "SLC16A2", "SLC16A4", "SLC16A6", "SLC16A9", "SLC16A13", "SLC16A1", "SLC16A3", "SLC16A14", "Age")

pall<- rcorr(MCT_TMM_transposed2, type = c("pearson"))

corrvalues<-pall[["r"]]
corpvalues<-pall[["P"]]

write.csv(corrvalues, "Maurer_Correlation r-values between MCTs and Age at diagnosis.csv")
write.csv(corpvalues, "Maurer_Correlation p-values between MCTs and Age at diagnosis.csv")


# Produce correlation r- and p-values for the ECM_TMM_transposed dataframe
# This requires conversion of dataframe into matrix
ECM_TMM_transposed<-as.matrix(ECM_TMM_transposed)

pall_ECM<- rcorr(ECM_TMM_transposed, type = c("pearson"))

corrvalues_ECM<-pall_ECM[["r"]]
corpvalues_ECM<-pall_ECM[["P"]]

write.csv(corrvalues_ECM, "Maurer_Correlation r-values between ECMs and Age at diagnosis.csv")
write.csv(corpvalues_ECM, "Maurer_Correlation p-values between ECMs and Age at diagnosis.csv")

