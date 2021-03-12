# input res
volcano_fun_p <- function(res, 
                          sig, 
                          fc
                          ) {
  
  require("ggplot2") #Best plots
  
  df <- as.data.frame(res)
  df$id <- rownames(df)
  
  # reduce too high significances
  df$pvalue <- ifelse(-log10(df$pvalue) > 10, 0.0000000001, df$pvalue)
  
  n_up <- nrow(subset(df, df$log2FoldChange > fc & df$pvalue < sig))
  n_down <- nrow(subset(df, df$log2FoldChange < -fc & df$pvalue < sig))
  
  df$col <- ifelse(abs(df$log2FoldChange) > fc & df$padj < sig, "Significant", "Not Significant")
  df$col[is.na(df$col)] <- "Not Significant"
  
  p1 <- ggplot(df, aes(log2FoldChange, -log10(pvalue))) + 
    geom_point(size = 0.5, aes(col=col)) +
    scale_color_manual(values=c("black", "red")) + 
    geom_hline(yintercept = -log10(sig), col="grey") +
    geom_vline(xintercept = c(-fc,fc), col="grey") +
    annotate("text", x = max(df$log2FoldChange), y = -log10(sig/2), label = n_up) +
    annotate("text", x = min(df$log2FoldChange), y = -log10(sig/2), label = n_down) +
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle(gsub(".*:\\s[A-Za-z0-9]*\\s", "", res@elementMetadata@listData$description[2]))
  
  p1
}


# GO Enrichment
gsea_go_viz <- function(geneList = genelist, 
                     go_class = "BP", 
                     n_terms = 20, 
                     outdir = gobp_dir, 
                     comp=comp,
                     orgdb = org.Dr.eg.db
                     ) {
  
  require("enrichplot") #Best plots
  
  res_gsea <- gseGO(geneList     = geneList,
                     OrgDb        = org.Dr.eg.db,
                     keyType      = "ENTREZID", 
                     ont          = go_class,
                     minGSSize    = 10,
                     maxGSSize    = 500,
                     pvalueCutoff = 1,
                     verbose      = FALSE)

  res_gsea <- setReadable(res_gsea, orgdb, 'ENTREZID')
  write.csv(res_gsea, paste0(outdir, paste(comp, collapse="_"), "_GO_enrichment.csv"))
  
  
  # Dotplot
  pdf(paste0(outdir, paste(comp, collapse="_"), "_GO_", go_class, "_dotplot.pdf"))
  print(dotplot(res_gsea, showCategory=n_terms, font.size=10, color = "pvalue", x = "NES"))
  dev.off()
  
  # Enrichment map
  emap <- enrichplot::pairwise_termsim(res_gsea, method = "JC", semData = NULL, showCategory = n_terms)
  pdf(paste0(outdir, paste(comp, collapse="_"), "_GO_", go_class, "_emap.pdf"))
  print(emapplot(emap, showCategory=n_terms, layout = "kk", color = "pvalue"))
  dev.off()
  
  # GSEA plots
  dir_gsea <- paste0(outdir, "GSEA_PLOTS/")
  dir.create(dir_gsea, showWarnings = FALSE, recursive = TRUE)
  for(i in 1:n_terms){
    pdf(paste0(dir_gsea, paste(comp, collapse="_"), "_GO_", go_class, "_gseaplot_", gsub(" |/", "_", substr(res_gsea$Description[i], 0, 15)),  ".pdf"))
    print(gseaplot2(res_gsea, geneSetID = i, title = res_gsea$Description[i]))
    dev.off()
  }
}


# KEGG Enrichment
gsea_kegg_viz <- function(geneList = genelist, 
                        n_terms = "20", 
                        outdir = kegg_dir, 
                        comp=comp,
                        org = "dre",
                        orgdb = org.Dr.eg.db
) {
  
  res_gsea <- gseKEGG(geneList   = geneList,
                    organism     = org,
                    minGSSize    = 10,
                    maxGSSize    = 500,
                    pvalueCutoff = 1,
                    verbose      = FALSE)
  
  res_gsea <- setReadable(res_gsea, orgdb, 'ENTREZID')
  write.csv(res_gsea, paste0(outdir, paste(comp, collapse="_"), "_KEGG_enrichment.csv"))
  
  
  # Dotplot
  pdf(paste0(outdir, paste(comp, collapse="_"), "_KEGG_", go_class, "_dotplot.pdf"))
  print(dotplot(res_gsea, showCategory=n_terms, font.size=10, color = "pvalue", x = "NES"))
  dev.off()
  
  # Enrichment map
  emap <- enrichplot::pairwise_termsim(res_gsea, method = "JC", semData = NULL, showCategory = n_terms)
  pdf(paste0(outdir, paste(comp, collapse="_"), "_KEGG_", go_class, "_emap.pdf"))
  print(emapplot(emap, showCategory=n_terms, layout = "kk", color = "pvalue"))
  dev.off()
  
  # GSEA plots
  dir_gsea <- paste0(outdir, "GSEA_PLOTS/")
  dir.create(dir_gsea, showWarnings = FALSE, recursive = TRUE)
  for(i in 1:n_terms){
    pdf(paste0(dir_gsea, paste(comp, collapse="_"), "_KEGG_", go_class, "_gseaplot_", gsub(" |/", "_", substr(res_gsea$Description[i], 0, 15)),  ".pdf"))
    print(gseaplot2(res_gsea, geneSetID = i, title = res_gsea$Description[i]))
    dev.off()
  }
}
