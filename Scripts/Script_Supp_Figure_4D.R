# Load packages
library(genefilter)
library(ggplot2)

# Set working directory
setwd("/Users/WIMM/Documents/Cristina_2018/Github/")

# Read transformed RPKM table
df <- read.table("Data/rpkm_table_log2.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Set factor levels
df$Group <- factor(df$Group, levels=c("CON", "AM", "KM", "AKM"))

# Identify significant gene sets
  # Perform F-tests
  ftest <- colFtests(as.matrix(df[,-1]), fac=df$Group, var.equal=FALSE)

  # Correct for multiple testing
  ftest$fdr <- p.adjust(ftest$p.value, method="fdr")

  # Subset for significant genes
  ftest <- ftest[which(ftest$fdr < 0.05),]
  ftest$Gene <- row.names(ftest)

# Retrieve significant genes
df <- df[, names(df) %in% c("Group", ftest$Gene)]

# PCA
    # Reduce dimensions
    pca <- prcomp(df[,-1], center=TRUE, scale.=TRUE)

    # Retrieve proportion of variance explained
    pca_var <- summary(pca)
    
    # Scatterplot
        # Definitions
        data <- data.frame(pca$x[,1], pca$x[,1])
        x <- pca$x[,1]
        y <- pca$x[,2]
        group <- df$Group
        xtitle <- paste("PC1 ", "(", signif(pca_var$importance[2,1]*100, 2), "%)", sep="")
        ytitle <- paste("PC2 ", "(", signif(pca_var$importance[2,2]*100, 2), "%)", sep="")
        legend.title <- "Group"
        color_group <- c("CON", "AM", "KM", "AKM")
        color <- c("blue", "red", "green", "purple")

        # Plot
        plot <- ggplot() +
           geom_point(data, mapping=aes(x=x, y=y, color=group), size=2) +
           labs(x=xtitle, y=ytitle) +
           scale_color_manual(breaks=color_group, values=color, name=legend.title) +
           theme_classic() +
           theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background=element_blank(),
                axis.line=element_line(colour="black"),
                axis.text=element_text(size=13),
                axis.title=element_text(size=18),
                axis.text.x=element_text(colour="black"),
                plot.title=element_text(size=18),
                legend.title=element_text(size=18),
                legend.text=element_text(size=18))


        # Save plot
        ggsave("Figures/Supp_Figure_4D.pdf", plot, width=7, height=7)
