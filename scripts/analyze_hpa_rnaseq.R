hpa_file <- "/home/quim/Databases/human_protein_atlas/rna_tissue.csv"
hk_file <- "/home/quim/project/tissue_specificity/housekeeping/tissue_specificity_rna_any_expressed.tab"
nrf2_file <- "/home/quim/project/tissue_specificity/nrf2ome/07172014-nrf2-sBC6cR.csv"


hpadata <- read.csv(hpa_file, header=TRUE) # Read data file
hkdata <- read.csv(hk_file, header=TRUE, sep="\t") # Read house-keeping data
nrf2data <- read.csv(nrf2_file, header=TRUE, sep=";") # Read nrf2 data

ensembls <- hkdata$Ensembl

hpa_hk <- hpadata[ hpadata$Gene %in% ensembls, ]

hpa_hk_liver <- hpa_hk[ hpa_hk$Sample == 'liver',]

value_liver <- hpa_hk_liver[hpa_hk_liver$Unit == "TPM", "Value"]

png(filename="rnaseq_expression_liver.png")
hist( value_liver,
      breaks = 50,
      main="RNAseq expression levels in HK genes in liver",
      xlab="RNAseq value (TPM)" )
dev.off()

value_tissues <- hpa_hk[hpa_hk$Unit == "TPM", "Value"]

png(filename="rnaseq_expression_all_tissues.png")
hist( value_tissues,
      breaks = 50,
      main="RNAseq expression levels in HK genes in all tissues",
      xlab="RNAseq value (TPM)" )
dev.off()


value_liver_less_1000 <- hpa_hk_liver[hpa_hk_liver$Value <= 1000, "Value"]
print(value_liver_less_1000)

png(filename="rnaseq_expression_liver_less_1000.png")
hist( value_liver_less_1000,
      breaks = 10,
      main="RNAseq expression levels <= 1000 TPM in HK genes in liver",
      xlab="RNAseq value (TPM)" )
dev.off()

value_liver_less_100 <- hpa_hk_liver[hpa_hk_liver$Value <= 100, "Value"]
print(value_liver_less_100)

png(filename="rnaseq_expression_liver_less_100.png")
hist( value_liver_less_100,
      breaks = 10,
      main="RNAseq expression levels <= 100 TPM in HK genes in liver",
      xlab="RNAseq value (TPM)" )
dev.off()

value_liver_less_10 <- hpa_hk_liver[hpa_hk_liver$Value <= 10, "Value"]
print(value_liver_less_10)

png(filename="rnaseq_expression_liver_less_10.png")
hist( value_liver_less_10,
      breaks = 10,
      main="RNAseq expression levels <= 10 TPM in HK genes in liver",
      xlab="RNAseq value (TPM)" )
dev.off()

