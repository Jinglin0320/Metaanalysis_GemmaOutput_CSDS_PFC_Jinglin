####Clean Up the GSEA table, pick out only significant genes from the metaanalysis###########
######Jinglin Xiong, Jun 23, 2025############################################################

library(dplyr)
library(stringr)  
# Read and filter the GSEA table
GSEA_Results <- read.csv("..../BAP_Updated_2025/Output/GSEA_Results.csv")
GSEA_Results_sig <- subset(GSEA_Results, padj < 0.05)

# Get the meta-analysis significant genes
sig_genes <- rownames(metaOutputSigNoNa) 

# Turn the leadingEdge column into lists of genes
dge_list <- strsplit(GSEA_Results_sig$leadingEdge, ",")  
dge_list <- lapply(dge_list, trimws)     
str(dge_list)
# List of 53
# $ : chr [1:46] "Egr2" "Fa2h" "Mobp" "Plp1" ...
# $ : chr [1:123] "Fos" "Ccl21b" "Dusp1" "Egr2" ...
# $ : chr [1:216] "Fos" "Ccn1" "Egr2" "Asb4" ...
# $ : chr [1:30] "Fa2h" "Nkx2-1" "Mobp" "Plp1" ...
# $ : chr [1:4] "Lcn2" "S100a8" "S100a9" "Slc30a7"
#.....

# Keep only genes that are also significant in my meta-analysis
dge_list_filt <- lapply(dge_list,
                         function(vec) vec[vec %in% sig_genes])
str(dge_list_filt)
# List of 53
# $ : chr [1:2] "Fa2h" "Lpar1"
# $ : chr [1:9] "Fos" "Dusp1" "Btg2" "Egr1" ...
# $ : chr [1:9] "Fos" "Ccn1" "Dnai3" "Sox9" ...
# $ : chr [1:3] "Fa2h" "Lpar1" "Sox9"
# $ : chr(0) 
# $ : chr [1:9] "Fos" "Otx2" "Egr1" "Zfp979" ...

# Add the filtered gene lists back to the data frame
GSEA_Results_sig$leadingEdge_sig <- sapply(dge_list_filt,
                                       paste, collapse = ",")
write.csv(GSEA_Results_sig, "..../Output/GSEA_Results_sig.csv")
