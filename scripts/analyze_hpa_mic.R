hpa_file <- "/home/quim/Databases/human_protein_atlas/normal_tissue.csv"
hk_file <- "/home/quim/project/tissue_specificity/housekeeping/tissue_specificity_rna_any_expressed.tab"
nrf2_file <- "/home/quim/project/tissue_specificity/nrf2ome/07172014-nrf2-sBC6cR.csv"


hpadata <- read.csv(hpa_file, header=TRUE) # Read data file
hkdata <- read.csv(hk_file, header=TRUE, sep="\t") # Read house-keeping data
nrf2data <- read.csv(nrf2_file, header=TRUE, sep=";") # Read nrf2 data

ensembls <- hkdata$Ensembl

hpa_hk <- hpadata[ hpadata$Gene %in% ensembls, ]

hpa_hk_liver <- hpa_hk[ hpa_hk$Tissue == 'liver',]

levels <- table(hpa_hk_liver$Level)
lvls <- c( levels[names(levels)=="Not detected"], levels[names(levels)=="Low"], levels[names(levels)=="Medium"], levels[names(levels)=="High"] ) 

png(filename="expression_liver.png")
barplot(lvls,
	main="Expression levels in HK genes in liver")
dev.off()

levels <- table(hpa_hk$Level)
lvls <- c( levels[names(levels)=="Not detected"], levels[names(levels)=="Low"], levels[names(levels)=="Medium"], levels[names(levels)=="High"] ) 
png(filename="expression_all_tissues.png")
barplot(lvls,
	main="Expression levels in HK genes in all tissues")
dev.off()

rel <- hpa_hk_liver$Reliability
lvl <- hpa_hk_liver$Level
xtabs(~ rel + lvl, hpa_hk_liver)

ensembls_nrf2_1 <- as.vector(nrf2data$source_speciesID) # Interactor 1
ensembls_nrf2_2 <- as.vector(nrf2data$target_speciesID) # Interactor 2
ensembls_nrf2 <- c(ensembls_nrf2_1, ensembls_nrf2_2) # Add the two columns in a single vector
ensembls_nrf2 <- unique(ensembls_nrf2) # Unify the repeated elements

vector <- c()
for (item in ensembls_nrf2){
	str <- strsplit(item, ",") # Split multiple ensembls in one row
	for (element in str){
		vector <- c(vector, element) # Add them in a separate vector
	}
}
ensembls_nrf2 <- c(ensembls_nrf2, vector) # Add the vector to the main vector
ensembls_nrf2 <- unique(ensembls_nrf2) # Unify the repeated elements

hpa_nrf2 <- hpadata[ hpadata$Gene %in% ensembls_nrf2, ] # Join HPA data with nrf2 data
hpa_nrf2_liver <- hpa_nrf2[ hpa_nrf2$Tissue == 'liver',] # Filter by liver
levels <- table(hpa_nrf2_liver$Level) # Classify by levels
lvls <- c( levels[names(levels)=="Not detected"], levels[names(levels)=="Low"], levels[names(levels)=="Medium"], levels[names(levels)=="High"] ) 
png(filename="nrf2_expression_liver.png")
barplot(lvls,
	main="Expression levels in nrf2 genes in liver")
dev.off()

