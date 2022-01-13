setwd("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/PCA analysis of ALL datasets_with filter applied/Revised analysis")

#Read all datasets with low expression genes removed
GSE119794<-read.table("F:/Latif Lab back up/AYSHE/PDAC/Differential expression analysis/GSE119794/GSE119794_cpm_data_filtered_geneIDversions removed.csv",sep = ",", header=T)
GSE79668<-read.table("F:/Latif Lab back up/AYSHE/PDAC/Other datasets for PCA/GSE79668_all 51 tumour cpm data_filtered.csv",sep = ",", header=T)
GSE93326<-read.table("F:/Latif Lab back up/AYSHE/PDAC/Other datasets for PCA/GSE93326_ALL samples_bulk tumour, epithelium and stroma_cpm_data_filtered.csv",sep = ",", header=T)
GSE131050_47bulk<-read.table("F:/Latif Lab back up/AYSHE/PDAC/Other datasets for PCA/GSE131050_47 FF bulk tumor_cpm_data_filtered_revised.csv",sep = ",", header=T)

#To identify the common genes (as long as row names are genes) between 4 datasets
common_genes <- intersect(intersect(intersect(rownames(GSE119794),rownames(GSE79668)),rownames(GSE93326)),rownames(GSE131050_47bulk))
write.table(common_genes, "Common genes identified among all datasets.csv", sep=",")