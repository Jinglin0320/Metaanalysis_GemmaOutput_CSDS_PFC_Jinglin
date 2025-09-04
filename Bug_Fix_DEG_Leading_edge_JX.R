#Pulling out validated genes
#Adapted from Dr. Megan Hagenauer
#Sept 3, 2025

setwd("/labs/twc/JinglinXiong/BAP_Updated_2025/Output")
fGSEA_Results<-read.csv("/labs/twc/JinglinXiong/BAP_Updated_2025/Output/GSEA_Results.csv", header=TRUE, stringsAsFactors = FALSE)

str(fGSEA_Results)
#'data.frame':	10847 obs. of  9 variables:

#Thresholding the fGSEA results by FDR:
fGSEA_Sig<-fGSEA_Results[fGSEA_Results$padj<0.05,]
nrow(fGSEA_Sig)

#Reading in some sort of data frame of DEGs from your meta-analysis:
MetaFDR05_wValidInfo_Exploratory_fGSEA <- read.csv("metaOutputFDR_orderedByPval.csv", header = TRUE, stringsAsFactors = FALSE)


str(MetaFDR05_wValidInfo_Exploratory_fGSEA)
#'data.frame':	21379 obs. of  8 variables:

MetaFDR05_wValidInfo_Exploratory_fGSEA <-na.omit(MetaFDR05_wValidInfo_Exploratory_fGSEA[MetaFDR05_wValidInfo_Exploratory_fGSEA$FDR < 0.05, ])

str(MetaFDR05_wValidInfo_Exploratory_fGSEA)
# 'data.frame':	133 obs. of  8 variables:

colnames(MetaFDR05_wValidInfo_Exploratory_fGSEA)



for (j in 1:length(fGSEA_Sig$leadingEdge)) {
  # Initialize vector to collect matching gene symbols
  MatchingGenes <- c()
  LeadingEdgeVector<-strsplit(fGSEA_Sig$leadingEdge[j], ",")
  print(LeadingEdgeVector[[1]])
  for (i in 1:nrow(MetaFDR05_wValidInfo_Exploratory_fGSEA)) {
    gene <- MetaFDR05_wValidInfo_Exploratory_fGSEA$X[i]
    GeneWithAnchors <- paste("^", gene, "$", sep = "")
    print(GeneWithAnchors)
    InGeneSet<- grepl(GeneWithAnchors, LeadingEdgeVector[[1]])
    print(InGeneSet)
    
    
    if(sum(InGeneSet>0)){
      MatchingGenes <- c(MatchingGenes, gene)
      print(MatchingGenes)
    }else{}
  }
  
  fGSEA_Sig$MatchedDEGs[j] <- paste(MatchingGenes, collapse = ";")
}


# View the new column
head(fGSEA_Sig$MatchedDEGs)


# Save the updated results
write.csv(fGSEA_Sig, "/labs/twc/JinglinXiong/BAP_Updated_2025/Output/Bug_Test_fGSEA_Sig_MatchedDEGs.csv")




