###### Do the top 50 genes#############

load("/oak/stanford/scg/lab_twc/JinglinXiong/BAP/Output/workspace.RData")
#####Top 50 genes#########
gene_top_50 <- rownames(metaOutputFDR_OrderbyPval)[1:50]
Rows_Interest<-MetaAnalysis_FoldChanges_ForMeta$x%in%gene_top_50
Log2FC_Subsetted <- MetaAnalysis_FoldChanges_ForMeta[Rows_Interest, ]
Log2FC_Subsetted_Matrix <- as.matrix(Log2FC_Subsetted[,-1])
row.names(Log2FC_Subsetted_Matrix)<-Log2FC_Subsetted$x
str(Log2FC_Subsetted_Matrix)

library(pheatmap)

pdf("/labs/twc/JinglinXiong/BAP/Output/Heatmap_TopMetaGenes_50_color_scale.pdf",
    height = 11, width = 8.5)

pheatmap(Log2FC_Subsetted_Matrix,
         color           = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         scale         = "none",
         breaks        = seq(-1.5, 1.5, length.out = 101),
         cluster_rows  = TRUE,
         cluster_cols  = TRUE,
         fontsize_row  = 8,
         fontsize_col  = 8,  
         width         = 8.5,
         height        = 11,
         border_color  = NA)

dev.off()

########With smaller width##############
pdf("/labs/twc/JinglinXiong/BAP/Output/Heatmap_TopMetaGenes_50_color_scale_narrow.pdf",
    height = 11, width = 5)

pheatmap(Log2FC_Subsetted_Matrix,
         color           = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         scale         = "none",
         breaks        = seq(-1.5, 1.5, length.out = 101),
         cluster_rows  = TRUE,
         cluster_cols  = TRUE,
         fontsize_row  = 8,
         fontsize_col  = 8,  
         width         = 8.5,
         height        = 11,
         border_color  = NA)

dev.off()

########All Genes#############

gene_all <- rownames(metaOutputFDR_OrderbyPval)[1:171]
Rows_Interest<-MetaAnalysis_FoldChanges_ForMeta$x%in%gene_all
Log2FC_Subsetted <- MetaAnalysis_FoldChanges_ForMeta[Rows_Interest, ]
Log2FC_Subsetted_Matrix <- as.matrix(Log2FC_Subsetted[,-1])
row.names(Log2FC_Subsetted_Matrix)<-Log2FC_Subsetted$x
str(Log2FC_Subsetted_Matrix)

library(pheatmap)

pdf("/labs/twc/JinglinXiong/BAP/Output/Heatmap_TopMetaGenes_all_color_scale.pdf",
    height = 11, width = 8.5)

pheatmap(Log2FC_Subsetted_Matrix,
         color           = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         scale         = "none",
         breaks        = seq(-1.5, 1.5, length.out = 101),
         cluster_rows  = TRUE,
         cluster_cols  = TRUE,
         fontsize_row  = 4,
         fontsize_col  = 8,  
         width         = 8.5,
         height        = 11,
         border_color  = NA)

dev.off()

pdf("/labs/twc/JinglinXiong/BAP/Output/Heatmap_TopMetaGenes_all_color_scale_narrow.pdf",
    height = 11, width = 5)

pheatmap(Log2FC_Subsetted_Matrix,
         color           = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         scale         = "none",
         breaks        = seq(-1.5, 1.5, length.out = 101),
         cluster_rows  = TRUE,
         cluster_cols  = TRUE,
         fontsize_row  = 4,
         fontsize_col  = 8,  
         width         = 8.5,
         height        = 11,
         border_color  = NA)

dev.off()

###########all genes, but in all comparisons#########
gene_set         <- rownames(metaOutputSigNoNa)          # 171 meta‑sig genes
present_mat      <- sapply(ListOfDEResults, \(x) gene_set %in% rownames(x$Log2FC))
genes_all6       <- gene_set[rowSums(present_mat) == ncol(present_mat)]  # keep #136

Rows_Interest<-MetaAnalysis_FoldChanges_ForMeta$x%in%genes_all6
Log2FC_Subsetted <- MetaAnalysis_FoldChanges_ForMeta[Rows_Interest, ]
Log2FC_Subsetted_Matrix <- as.matrix(Log2FC_Subsetted[,-1])
row.names(Log2FC_Subsetted_Matrix)<-Log2FC_Subsetted$x
str(Log2FC_Subsetted_Matrix)

library(pheatmap)

pdf("/labs/twc/JinglinXiong/BAP/Output/Heatmap_TopMetaGenes_all_color_scale_allcomparison.pdf",
    height = 11, width = 8.5)

pheatmap(Log2FC_Subsetted_Matrix,
         color           = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         scale         = "none",
         breaks        = seq(-1.5, 1.5, length.out = 101),
         cluster_rows  = TRUE,
         cluster_cols  = TRUE,
         fontsize_row  = 4,
         fontsize_col  = 8,  
         width         = 8.5,
         height        = 11,
         border_color  = NA)

dev.off()


pdf("/labs/twc/JinglinXiong/BAP/Output/Heatmap_TopMetaGenes_all_color_scale_allcomparison_narrow.pdf",
    height = 11, width = 5)

pheatmap(Log2FC_Subsetted_Matrix,
         color           = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         scale         = "none",
         breaks        = seq(-1.5, 1.5, length.out = 101),
         cluster_rows  = TRUE,
         cluster_cols  = TRUE,
         fontsize_row  = 4,
         fontsize_col  = 8,  
         width         = 8.5,
         height        = 11,
         border_color  = NA)

dev.off()


# 1. Select the top 50 genes by FDR (adjust the column name to match your data: e.g. "FDR", "padj", etc.)
#    Here I assume your metaOutputSigNoNa has a column called 'FDR'
fdr_vec <- metaOutputSigNoNa[ , "FDR" ]          # change name if needed

top50_genes <- rownames(metaOutputSigNoNa)[
  order(fdr_vec, decreasing = FALSE)
][ seq_len( min(50, nrow(metaOutputSigNoNa)) ) ]

# 2. Of those top 50, keep only the ones present in every dataset in ListOfDEResults
present_mat <- sapply(ListOfDEResults, function(res) {
  top50_genes %in% rownames(res$Log2FC)
})
# genes present in all comparisons:
genes_to_plot <- top50_genes[ rowSums(present_mat) == length(ListOfDEResults) ]

# (If you prefer to force exactly 50 even if some drop out, you could skip the presence filter 
#  or widen it to e.g. present in ≥X datasets.)

# 3. Subset your fold‐change table
Rows_Interest <- MetaAnalysis_FoldChanges_ForMeta$x %in% genes_to_plot
Log2FC_Subsetted <- MetaAnalysis_FoldChanges_ForMeta[Rows_Interest, ]
Log2FC_Matrix <- as.matrix(Log2FC_Subsetted[ , -1])
rownames(Log2FC_Matrix) <- Log2FC_Subsetted$x

# 4. Draw the heatmap
library(pheatmap)
pdf("/labs/twc/JinglinXiong/BAP/Output/Heatmap_Top50_MetaGenes_AllComparisons.pdf",
    height = 11, width = 8.5)
pheatmap(
  Log2FC_Matrix,
  color          = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  scale          = "none",
  breaks         = seq(-1.5, 1.5, length.out = 101),
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  fontsize_row   = 8,
  fontsize_col   = 8,
  border_color   = NA
)
dev.off()

pdf("/labs/twc/JinglinXiong/BAP/Output/Heatmap_Top50_MetaGenes_AllComparisons_narrow.pdf",
    height = 11, width = 5)
pheatmap(
  Log2FC_Matrix,
  color          = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  scale          = "none",
  breaks         = seq(-1.5, 1.5, length.out = 101),
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  fontsize_row   = 8,
  fontsize_col   = 8,
  border_color   = NA
)
dev.off()
