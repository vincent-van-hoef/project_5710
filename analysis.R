# This is the analysis script for project #5710. This should be sourced before compiling the .Rmd report.

# Load packages
suppressMessages(library("ggplot2"))
suppressMessages(library("tidyr"))
suppressMessages(library("DESeq2"))
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("org.Dr.eg.db"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("enrichplot"))

# Set working directory to script location
proj_dir <- dirname(sys.frame(1)$ofile)
setwd(proj_dir)

# Load several custom functions
source("helpers.R")

#############
# Load data #
#############

# Set up a report folder, remove existing ones first
res_dir <- paste0(proj_dir, "/Report/")
unlink(res_dir, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE)

# Load datasets, make sure sample names are the same in data and meta
cts   <- read.csv(paste0(proj_dir, "/data/merged_gene_counts_rest.txt"), sep = "\t", row.names = 1, check.names = FALSE)
meta   <- read.csv(paste0(proj_dir, "/data/20210302_metadata.csv"), sep = ";", row.names = 1)
meta$Group <- paste0("Group_", meta$Group )
meta$Genotype <- gsub("-", "", meta$Genotype)
meta$samples     <- rownames(meta)

# make sure samples are in the same order in metadata and count data
cts    <- cts[,rownames(meta)]
if(!identical(colnames(cts), rownames(meta))){stop("Samplenames do not match between counts and metadata")}

# Make and Save DESeq object
dds <- DESeqDataSetFromMatrix(countData = cts, colData = meta, design = ~ 1)

##############
# Pre Filter #
##############

# Remove all zero genes
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

###################
# Quality Control #
###################

# Create a QC folder in the Report dir
qc_dir <- paste0(res_dir, "QC/")
dir.create(qc_dir, showWarnings = FALSE)

# Create barplot of the total number of reads per sample
pdf(paste0(qc_dir, "barplot_raw_counts.pdf"))
par(mar=c(7,6,2,2), mgp=c(5,1,0))
barplot(colSums(counts(dds)), 
        col = ifelse(colData(dds)$Genotype == "wildtype", "green", "red"),
        ylab = "Number of mapped reads", las =2)
abline(h = mean(colSums(counts(dds))))
dev.off()

# Create boxplot of logged raw counts: are there big differences in the library size?
pdf(paste0(qc_dir, "boxplot_raw_counts.pdf"))
par(mar=c(6,8,2,2))
boxplot(log2(counts(dds)+1), 
        col = ifelse(colData(dds)$Genotype == "wildtype", "green", "red"),
        pch = ".",
        horizontal = TRUE,
        las = 1,
        xlab = "log2(Counts +1)")
dev.off()

# PCA plot of 100 most variable genes: how is the data structured?
dds   <- estimateSizeFactors(dds)
rld   <- rlog(dds)
p1    <- plotPCA(rld, intgroup = "Group", ntop = 100)
ggsave(paste0(qc_dir, "PCA_rlog_top100.pdf"), plot = p1)
p2    <- plotPCA(rld, intgroup = "Group", ntop = 500)
ggsave(paste0(qc_dir, "PCA_rlog_top500.pdf"), plot = p2)
p3    <- plotPCA(rld, intgroup = "Group", ntop = 1000)
ggsave(paste0(qc_dir, "PCA_rlog_top1000.pdf"), plot = p3)

# Create boxplot of logged raw counts: are there big differences in the library size?
pdf(paste0(qc_dir, "boxplot_normalized_counts.pdf"))
par(mar=c(6,8,2,2))
boxplot(assay(rld), 
        col = ifelse(colData(dds)$Genotype == "wildtype", "green", "red"),
        pch = ".",
        horizontal = TRUE,
        las = 1,
        xlab = "Counts (normalized by rlog)")
dev.off()

# Sample Distances: how are the samples related to each other?
sampleDists                 <- dist(t(assay(rld)))
sampleDistMatrix            <- as.matrix(sampleDists)
rownames(sampleDistMatrix)  <- paste(meta$samples, meta$number.of.eggs, sep = "_" )
colnames(sampleDistMatrix)  <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf(paste0(qc_dir, "HeatmapDistances.pdf"))
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

###########################
# Differential Expression #
###########################

# Create a differential expression result folder
de_dir <- paste0(res_dir, "Differential_Expression/")
dir.create(de_dir, showWarnings = FALSE)

# Exclude certain samples from downstream analysis
colData(dds)$Group1_Excl <- ifelse(colData(dds)$samples == "TG-2625-C3", "Exclude", as.vector(colData(dds)$Group))
colData(dds)$Group2_Excl <- ifelse(colData(dds)$samples == "TG-2625-M6", "Exclude", as.vector(colData(dds)$Group))
colData(dds)$Group12_Excl <- ifelse(colData(dds)$samples %in% c("TG-2625-C3", "TG-2625-M6"), "Exclude", as.vector(colData(dds)$Group))

designList   <- list(Genotype = list(Model = "~ 0 + Genotype",
                                       Contrast = list("Mutant_vs_WT" = c("Genotype", "mutant", "wildtype"))),
                     EggNumber = list(Model = "~ 0 + number.of.eggs",
                                      Contrast = list("High_vs_Low" = c("number.of.eggs", "high", "low"))),
                     Groups_AllSamp = list(Model = "~ 0 + Group",
                                   Contrast = list("Group1_vs_Group2" = c("Group", "Group_1", "Group_2"),
                                                       "Group1_vs_Group3" = c("Group", "Group_1", "Group_3"),
                                                       "Group2_vs_Group3" = c("Group", "Group_2", "Group_3"))),
                     Group1_Excl = list(Model = "~ 0 + Group1_Excl",
                                        Contrast = list("Group1_Excl_vs_Group2" = c("Group1_Excl", "Group_1", "Group_2"))),
                     Group2_Excl = list(Model = "~ 0 + Group2_Excl",
                                        Contrast = list("Group1_vs_Group2_Excl" = c("Group2_Excl", "Group_1", "Group_2"))),
                     Group12_Excl = list(Model = "~ 0 + Group12_Excl",
                                        Contrast = list("Group1_Excl_vs_Group2_Excl" = c("Group12_Excl", "Group_1", "Group_2"))))

# Run DE and plotting function over each comparison in the designlist
for (design in names(designList)){
  
  design_dir <- paste0(de_dir, design, "/")
  dir.create(design_dir, showWarnings = FALSE)
  
  # Update design and run DESeq
  dds <- DESeqDataSet(dds, as.formula(designList[[design]][["Model"]]))
  dds <- DESeq(dds)
  
  # Create results for each contrast in design parameter
  for(contrast in names(designList[[design]][["Contrast"]])) {
    
    # Make separate folder for each contrast
    contrast_dir <- paste0(design_dir, contrast, "/")
    dir.create(contrast_dir, showWarnings = FALSE)
    
    # create DE folder in each contrast folder
    diff_exp_dir <- paste0(contrast_dir, "DE/")
    dir.create(diff_exp_dir, showWarnings = FALSE)
    
    # Extract the comparison for plotting purposes
    comp <- designList[[design]][["Contrast"]][[contrast]]

    # Extract results
    res <- results(dds, contrast = comp)

    # Convert ensembl id to symbols
    res$Symbol <- mapIds(org.Dr.eg.db, keys = rownames(res), keytype = "ENSEMBL", column = "SYMBOL")
    res <- res[order(res$padj),]
    write.csv(res, paste0(diff_exp_dir, paste(comp, collapse="_"), ".csv"))
    
    # Create volcano plot
    png(paste0(diff_exp_dir, contrast, "_volcano.png"))
    p1 <- volcano_fun_p(res, fc = 1, sig = 0.05)
    print(p1)
    dev.off()

    ##################
    # GSEA Enrichment#
    ##################

    # Log2Fc statistic
    ##################
    
    genelist <- sort(setNames(res$log2FoldChange, rownames(res)), decreasing = TRUE)
    
    # genelist names should be entrezid
    names(genelist) <- mapIds(org.Dr.eg.db, 
                              keys = names(genelist), 
                              column = "ENTREZID", 
                              keytype = "ENSEMBL")
    genelist <- genelist[!is.na(names(genelist))]
    
    # GO BP Enrichment
    gobp_dir <- paste0(contrast_dir, "GSEA/Log2FC/GO/BP/")
    dir.create(gobp_dir, showWarnings = FALSE, recursive = TRUE)
  
    gsea_go_viz(geneList = genelist, 
                go_class = "BP", 
                n_terms = 20, 
                outdir = gobp_dir, 
                comp=comp,
                orgdb = org.Dr.eg.db)
    
    # GO MF Enrichment
    gomf_dir <- paste0(contrast_dir, "GSEA/Log2FC/GO/MF/")
    dir.create(gomf_dir, showWarnings = FALSE, recursive = TRUE)
    
    gsea_go_viz(geneList = genelist, 
                go_class = "MF", 
                n_terms = 20, 
                outdir = gomf_dir, 
                comp=comp,
                orgdb = org.Dr.eg.db)
    
    # GO CC Enrichment
    gocc_dir <- paste0(contrast_dir, "GSEA/Log2FC/GO/CC/")
    dir.create(gocc_dir, showWarnings = FALSE, recursive = TRUE)
    
    gsea_go_viz(geneList = genelist, 
                go_class = "CC", 
                n_terms = 20, 
                outdir = gocc_dir, 
                comp=comp,
                orgdb = org.Dr.eg.db)
    
    # KEGG Enrichment
    kegg_dir <- paste0(contrast_dir, "GSEA/Log2FC/KEGG/")
    dir.create(kegg_dir, showWarnings = FALSE, recursive = TRUE)
    
    gsea_kegg_viz(geneList = genelist, 
                  n_terms = 20, 
                  outdir = kegg_dir, 
                  comp=comp,
                  org = "dre",
                  orgdb = org.Dr.eg.db)
    
    # sign(Log2Fc)*-log10(Pvalue)
    #############################

    genelist <- sort(setNames(sign(res$log2FoldChange)*-log10(res$pvalue), rownames(res)), decreasing = TRUE)
    # genelist names should be entrezid
    names(genelist) <- mapIds(org.Dr.eg.db, 
                                keys = names(genelist), 
                                column = "ENTREZID", 
                                keytype = "ENSEMBL")
    genelist <- genelist[!is.na(names(genelist))]
    
    # GO BP Enrichment
    gobp_dir <- paste0(contrast_dir, "GSEA/signed_Pval/GO/BP/")
    dir.create(gobp_dir, showWarnings = FALSE, recursive = TRUE)
    
    gsea_go_viz(geneList = genelist, 
                go_class = "BP", 
                n_terms = 20, 
                outdir = gobp_dir, 
                comp=comp,
                orgdb = org.Dr.eg.db)
    
    # GO MF Enrichment
    gomf_dir <- paste0(contrast_dir, "GSEA/signed_Pval/GO/MF/")
    dir.create(gomf_dir, showWarnings = FALSE, recursive = TRUE)
    
    gsea_go_viz(geneList = genelist, 
                go_class = "MF", 
                n_terms = 20, 
                outdir = gomf_dir, 
                comp=comp,
                orgdb = org.Dr.eg.db)
    
    # GO CC Enrichment
    gocc_dir <- paste0(contrast_dir, "GSEA/signed_Pval/GO/CC/")
    dir.create(gocc_dir, showWarnings = FALSE, recursive = TRUE)
    
    gsea_go_viz(geneList = genelist, 
                go_class = "CC", 
                n_terms = 20, 
                outdir = gocc_dir, 
                comp=comp,
                orgdb = org.Dr.eg.db)
    
    # KEGG Enrichment
    kegg_dir <- paste0(contrast_dir, "GSEA/signed_Pval/KEGG/")
    dir.create(kegg_dir, showWarnings = FALSE, recursive = TRUE)
    
    gsea_kegg_viz(geneList = genelist, 
                  n_terms = 20, 
                  outdir = kegg_dir, 
                  comp=comp,
                  org = "dre",
                  orgdb = org.Dr.eg.db)
  }
}