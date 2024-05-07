####IMPORTANT! Please download the github version to use the most recent version of TCGAbiolinks
#devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
###Use devtools::install_github(repo = "ELELAB/TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(recount)
###conversion of uuids to TCGA barcodes
#library(TCGAutils)

# have to exclude
# adrenal gland
# adipose tissue
# blood vessel
# bone marrow
# cervix uteri
# fallopian tube
# salivary gland
# small intestine
gtex_datasets = c("bladder","blood", "brain", "breast", "colon", "esophagus", "heart", "kidney", "liver", "lung", "muscle", "nerve", "ovary","pancreas", "pituitary", "prostate", "skin", "spleen", "stomach", "testis", "thyroid", "uterus", "vagina")
saved_gtex_datasets = vector()

min_size = 100

for (tissue in gtex_datasets) {
	raw_data <-TCGAquery_recount2(project="GTEX", tissue=tissue)[[1]]
	n_samples <- dim(raw_data)[2]
	print(paste(tissue, n_samples, sep='    '))
	if (n_samples >= min_size) {
		saved_gtex_datasets <- append(saved_gtex_datasets, tissue)
		save(raw_data, file=paste0('GTEX_', tissue, '.Rdata'))	
	}
}

write.table(saved_gtex_datasets, file='saved_gtex_datasets.txt', sep='\t')

