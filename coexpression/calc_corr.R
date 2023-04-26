library(TCGAbiolinks)
library(biomaRt)
library(recount)

sessionInfo()

corr_co         = 0.6          # correlation cut_off
dataNorm_method = "gcContent"  # data normalisation method for dataNorm function
dataFilt_qnco   = 0.25         # quantile cut-off for dataFilt function
out_dir         = "./"
platform        = "gtex"
avail_datasets  = "saved_gtex_datasets.txt"

# load list of downloaded tissues
tissues = read.table(avail_datasets)[[1]]

# Download correspondence between Ensembl ID and HGCN symbols
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listOfGenes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "gene_biotype"),
                     filters = c("biotype"),values = list(biotype="protein_coding"), 
                     mart = mart)

for (tissue in tissues) {

	print(paste("Now doing", tissue))

	# Download GTEX data
	raw_data <- get(load(paste0('GTEX_', tissue, '.Rdata')))

	# Scale counts
	data <- scale_counts(raw_data, round = TRUE)

	# remove isoforms from Ensembl IDs
	rownames(data)<-gsub("\\..*", "",rownames(data))

	# keep only the genes that are in the protein coding gene list
	data <- data[is.element(rownames(data), listOfGenes$ensembl_gene_id),]

	# normalisation to remove bias for GC-enriched genes (GC content) and others,
	dataNorm <- TCGAanalyze_Normalization(tabDF = data,
                                      geneInfo = geneInfoHT,
                                      method = dataNorm_method)                

	# quantile filtering - remove genes that are very lowly expressed in many samples
	dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  dataFilt_qnco)

	# switch rownames to gene names according to the previously downloaded matching
	rownames(dataFilt) <- listOfGenes$hgnc_symbol[match(rownames(dataFilt), listOfGenes$ensembl_gene_id)]

	# transpose gene matrix and calculate correlation between expression profiles of genes
	corrs = cor(t(dataFilt), method='spearman')

	# coexpression analysis starts
	out_fname = paste(out_dir, "corrs_", tissue, ".txt", sep="")
	write.table(corrs, file=out_fname, quote=FALSE, sep="\t", row.names=FALSE, col.names = TRUE)
}
