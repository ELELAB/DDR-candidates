library(TCGAbiolinks)
library(SummarizedExperiment)
library(recount)
library(biomaRt)
library(CEMiTool)
library(stringr)
library(limma)
###conversion of uuids to TCGA barcodes
#library(TCGAutils)

corr_co         = 0.6          # correlation cut_off
dataNorm_method = "gcContent"  # data normalisation method for dataNorm function
dataFilt_qnco   = 0.25         # quantile cut-off for dataFilt function
out_dir         = "./"
platform        = "gtex"
t               = "uterus"

# function to generate fictional barcode for GTEX studies
generate_barcode <- function(i) {
    #TCGA-93-A4JP-01A-11H-A24S-13
    return(paste('TCGA', '00', '0000', '00A', '00D', str_pad(i, 4, pad='0'), '00', sep='-' ))
}

# Download correspondence between Ensembl ID and HGCN symbols
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listOfGenes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "gene_biotype"),
                     filters = c("biotype"),values = list(biotype="protein_coding"), 
                     mart = mart)



# load saved datasets
datasets = read.table('saved_gtex_datasets.txt')

for (dataset in datasets$x) {

	print(paste("currently on", dataset))

	# Load GTEX data
	raw_data <- get(load(paste0('GTEX_', dataset, '.Rdata')))

	# Scale counts
	data <- scale_counts(raw_data, round = TRUE)

	# remove isoforms from Ensembl IDs
	rownames(data)<-gsub("\\..*", "",rownames(data))

	# keep only the genes that are in the protein coding gene list
	data <- data[is.element(rownames(data), listOfGenes$ensembl_gene_id),]

	# add fictional barcode
	data$barcode <- sapply(1:dim(colData(data))[1], generate_barcode)

	# add tumor definition
	data$definition <- rep('Normal', dim(colData(data))[1])

	# CEMiTool pipeline starts here
	dataPrep <- TCGAanalyze_Preprocessing(object = data, cor.cut = corr_co, filename = NULL)

	# normalisation to remove bias for GC-enriched genes (GC content) and others,
	dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = dataNorm_method)                

	# quantile filtering - remove genes that are very lowly expressed in many samples
	dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  dataFilt_qnco)

	# switch rownames to gene names according to the previously downloaded matching
	rownames(dataFilt) <- listOfGenes$hgnc_symbol[match(rownames(dataFilt), listOfGenes$ensembl_gene_id)]

	# data normalised with the voom methodology
	dataFilt_voom = voom(dataFilt)$E
	dataFilt_ok <- as.data.frame(dataFilt_voom)

	#here the coexpression analysis starts
	cem <- cemitool(dataFilt_ok)
	cem #to have basic information on the modules and selected genes
	# inspect modules
	print(nmodules(cem))
	modules <- module_genes(cem)
	print(head(modules))
	out_fname = paste(out_dir, "modules_", dataset, ".txt", sep="")
	write.table(modules, file=out_fname, sep="\t", quote = TRUE)
}
