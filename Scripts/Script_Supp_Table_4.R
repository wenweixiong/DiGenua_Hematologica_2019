# Load packages
library(org.Mm.eg.db)
library(GO.db)
library(GOstats)

# Set working directory
setwd("/Users/WIMM/Documents/Cristina_2018/Github/")

# Retrieve Entrez Gene ID for all genes
    # Read comparison file
    diff <- read.table("Data/comparison_analysis_AEKvsAE_all_genes_new.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

    # Retrieve
    ID_all <- select(org.Mm.eg.db, keys=diff$Geneid, columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL")

    # Remove non-matches
    ID_all <- ID_all[which(!is.na(ID_all$ENTREZID)), "ENTREZID"]

# Subset significant genes
up <- diff[which(diff$FDR < 0.05 & diff$logFC > 0), "Geneid"]

# Retrieve Entrez Gene ID for significant genes
    # Retrieve
    ID <- select(org.Mm.eg.db, keys=up, columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL")

    # Remove non-matches
    ID <- ID[which(!is.na(ID$ENTREZID)), "ENTREZID"]

# Test of over-representation
    # Set up parameters for analysis
    params <- new('GOHyperGParams', geneIds=ID, universeGeneIds=ID_all, annotation='org.Mm.eg.db', ontology="BP", pvalueCutoff=0.05, conditional=TRUE, testDirection="over")

    # Analyse
    go <- hyperGTest(params)

    # Generate result table
    go.table <- summary(go)

    # Filter significant terms after adjustment
    go.table$bonferroni <- p.adjust(go.table$Pvalue, method="bonferroni")
    
    # Indicate no. of genes up-regulated
    go.table$Total_genes_upregulated <- length(up)
    
    # Calculate precentage of hits
    go.table$Count_pct <- round((go.table$Count/go.table$Total_genes_upregulated)*100, 2)
    
    # Include comparison label
    go.table$Comparison <- "AKMvsAM"
    
    # Reorder columns
    go.table <- go.table[, c(1, 3:5, 9:10, 6:7, 2, 8, 11)]

# Keep adjusted p-values < 0.05
go.table <- go.table[which(go.table$bonferroni < 0.05), ]

# Write file
write.table(go.table, "Tables/Supp_Table_4.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
