# Messing around with meta-analysis of four studies on CSDS
# Jinglin Xiong
# Aug 1, 2022

#6) Combine together the relevant results from different studies into a single 
#data-frame for the effect sizes (Log2FC) and a single data.frame for the sampling variances. 

#7) Make a correlation matrix to compare the overall results from the different studies.
#Further visualize the comparison using a hierarchically-clustered heatmap of the correlation matrix. 
#Which studies seem to have similar results? 

#8) Run a meta-analysis using all of the effect sizes for each gene that has data from at least 2 studies. 

#9) Correct the meta-analysis output to take into account the fact that we are running
#the statistical calculations many times and therefore have a heightened risk of false discovery 
#(false discovery rate correction) 

#10) Determine which are the top differentially expressed genes and create 
#forest plots to visualize the effect sizes for those top differentially expressed 
#genes across the different studies. 


#Code packages used (may require installation & installation of dependencies):
load("/oak/stanford/scg/lab_twc/JinglinXiong/BAP_Updated_2025/Output/DEextraction_Workspace.RData")
library(plyr)
library(metafor)
library(reshape)
library(multtest)

#Set working directory:
setwd("/oak/stanford/scg/lab_twc/JinglinXiong/BAP_Updated_2025/Output")

#6) Combine together the relevant results from different studies into a single data-frame 
#for the effect sizes (Log2FC) and a single data.frame for the sampling variances. 

#The Log2FC values are the first element in the differential expression result object for each study
#The gene symbols are the row.names - 75% of gene symbols for rats and mice are the same, 
#so for simplicity sake, instead of using a gene orthology database to inform equivalency, 
#we're just going to align the datasets by gene symbol.


AligningDatasets<-function(ListOfDEResults){
  
  MetaAnalysis_FoldChange_Dfs<-list()
  
  for(i in c(1:length(ListOfDEResults))){
    MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(x=row.names(ListOfDEResults[[i]][[1]]),ListOfDEResults[[i]][[1]], stringsAsFactors=FALSE)
  }
  
  print("MetaAnalysis_FoldChange_Dfs:")
  print(str(MetaAnalysis_FoldChange_Dfs))
  
  MetaAnalysis_FoldChanges<<-join_all(MetaAnalysis_FoldChange_Dfs, by="x", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("MetaAnalysis_FoldChanges:")
  print(str(MetaAnalysis_FoldChanges))
  
  MetaAnalysis_SV_Dfs<-list()
  
  for(i in c(1:length(ListOfDEResults))){
    MetaAnalysis_SV_Dfs[[i]]<-data.frame(x=row.names(ListOfDEResults[[i]][[4]]),ListOfDEResults[[i]][[4]], stringsAsFactors=FALSE)
  }
  
  print("MetaAnalysis_SV_Dfs:")
  print(str(MetaAnalysis_SV_Dfs))
  
  MetaAnalysis_SV<<-join_all(MetaAnalysis_SV_Dfs, by="x", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("MetaAnalysis_SV:")
  print(str(MetaAnalysis_SV))
  
  rm(MetaAnalysis_SV_Dfs, MetaAnalysis_FoldChange_Dfs)
}


#Example Usage;

ListOfDEResults<-list(DEResults_GSE109315, DEResults_GSE81672, DEResults_GSE85472,DEResults_GSE81587,DEResults_GSE102556.1,DEResults_GSE151807)



AligningDatasets(ListOfDEResults)
# [1] "MetaAnalysis_FoldChange_Dfs:"
# List of 6
# $ :'data.frame':	32953 obs. of  3 variables:
#   ..$ x                                     : chr [1:32953] "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# ..$ GSE109315_StressResilient_Vs_Control  : num [1:32953] -0.084 -0.0858 0.0425 -0.408 0.2268 ...
# ..$ GSE109315_StressSusceptible_Vs_Control: num [1:32953] 0.1929 -0.0733 -0.0991 -0.0613 0.3685 ...
# $ :'data.frame':	21628 obs. of  3 variables:
#   ..$ x                                    : chr [1:21628] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ GSE81672_StressResilient_Vs_Control  : num [1:21628] -0.7165 -0.0636 0.6613 0.1336 0.0292 ...
# ..$ GSE81672_StressSusceptible_Vs_Control: num [1:21628] -1.667 -0.0998 -1.264 -0.2959 -0.0832 ...
# $ :'data.frame':	20237 obs. of  2 variables:
#   ..$ x                       : chr [1:20237] "0610009B22Rik" "0610009L18Rik" "0610010K14Rik" "0610012G03Rik" ...
# ..$ GSE85472_CSDS_vs_Control: num [1:20237] -0.03498 0.03654 -0.116 0.01691 0.00807 ...
# $ :'data.frame':	23318 obs. of  2 variables:
#   ..$ x                       : chr [1:23318] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ GSE81587_CUMS_vs_Control: num [1:23318] -1.536 0.0853 0.1828 0.1766 -0.0161 ...
# $ :'data.frame':	24789 obs. of  2 variables:
#   ..$ x                          : chr [1:24789] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ GSE102556.1_CUMS_vs_Control: num [1:24789] -0.3705 0.00622 -0.2275 0.2992 0.06806 ...
# $ :'data.frame':	18537 obs. of  2 variables:
#   ..$ x                        : chr [1:18537] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010K14Rik" ...
# ..$ GSE151807_CUMS_vs_Control: num [1:18537] 0.033 0.1549 -0.1254 -0.036 0.0529 ...
# NULL
# [1] "MetaAnalysis_FoldChanges:"
# 'data.frame':	35210 obs. of  9 variables:
#   $ x                                     : chr  "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# $ GSE109315_StressResilient_Vs_Control  : num  -0.084 -0.0858 0.0425 -0.408 0.2268 ...
# $ GSE109315_StressSusceptible_Vs_Control: num  0.1929 -0.0733 -0.0991 -0.0613 0.3685 ...
# $ GSE81672_StressResilient_Vs_Control   : num  -0.7165 NA -0.0636 0.6613 0.1336 ...
# $ GSE81672_StressSusceptible_Vs_Control : num  -1.667 NA -0.0998 -1.264 -0.2959 ...
# $ GSE85472_CSDS_vs_Control              : num  NA NA -0.035 NA 0.0365 ...
# $ GSE81587_CUMS_vs_Control              : num  -1.536 NA 0.0853 0.1828 0.1766 ...
# $ GSE102556.1_CUMS_vs_Control           : num  -0.3705 NA 0.00622 -0.2275 0.2992 ...
# $ GSE151807_CUMS_vs_Control             : num  0.033 NA 0.155 NA -0.125 ...
# NULL
# [1] "MetaAnalysis_SV_Dfs:"
# List of 6
# $ :'data.frame':	32953 obs. of  3 variables:
#   ..$ x                                     : chr [1:32953] "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# ..$ GSE109315_StressResilient_Vs_Control  : num [1:32953] 0.1704 0.0202 0.0063 0.1143 0.0344 ...
# ..$ GSE109315_StressSusceptible_Vs_Control: num [1:32953] 0.09173 0.01089 0.00339 0.06156 0.01856 ...
# $ :'data.frame':	21628 obs. of  3 variables:
#   ..$ x                                    : chr [1:21628] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ GSE81672_StressResilient_Vs_Control  : num [1:21628] 0.77289 0.00638 0.23854 0.04249 0.011 ...
# ..$ GSE81672_StressSusceptible_Vs_Control: num [1:21628] 0.91575 0.00756 0.28301 0.0504 0.01303 ...
# $ :'data.frame':	20237 obs. of  2 variables:
#   ..$ x                       : chr [1:20237] "0610009B22Rik" "0610009L18Rik" "0610010K14Rik" "0610012G03Rik" ...
# ..$ GSE85472_CSDS_vs_Control: num [1:20237] 0.00185 0.00309 0.00189 0.00142 0.0056 ...
# $ :'data.frame':	23318 obs. of  2 variables:
#   ..$ x                       : chr [1:23318] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ GSE81587_CUMS_vs_Control: num [1:23318] 0.6262 0.0111 0.2658 0.0684 0.0259 ...
# $ :'data.frame':	24789 obs. of  2 variables:
#   ..$ x                          : chr [1:24789] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ GSE102556.1_CUMS_vs_Control: num [1:24789] 0.3037 0.0344 0.2408 0.0368 0.0304 ...
# $ :'data.frame':	18537 obs. of  2 variables:
#   ..$ x                        : chr [1:18537] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010K14Rik" ...
# ..$ GSE151807_CUMS_vs_Control: num [1:18537] 0.02327 0.00238 0.00343 0.00324 0.00259 ...
# NULL
# [1] "MetaAnalysis_SV:"
# 'data.frame':	35210 obs. of  9 variables:
#   $ x                                     : chr  "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# $ GSE109315_StressResilient_Vs_Control  : num  0.1704 0.0202 0.0063 0.1143 0.0344 ...
# $ GSE109315_StressSusceptible_Vs_Control: num  0.09173 0.01089 0.00339 0.06156 0.01856 ...
# $ GSE81672_StressResilient_Vs_Control   : num  0.77289 NA 0.00638 0.23854 0.04249 ...
# $ GSE81672_StressSusceptible_Vs_Control : num  0.91575 NA 0.00756 0.28301 0.0504 ...
# $ GSE85472_CSDS_vs_Control              : num  NA NA 0.00185 NA 0.00309 ...
# $ GSE81587_CUMS_vs_Control              : num  0.6262 NA 0.0111 0.2658 0.0684 ...
# $ GSE102556.1_CUMS_vs_Control           : num  0.3037 NA 0.0344 0.2408 0.0368 ...
# $ GSE151807_CUMS_vs_Control             : num  0.02327 NA 0.00238 NA 0.00343 ...
# NULL
####################################

#7) Make a correlation matrix to compare the overall results from the different studies. 
#Further visualize the comparison using a hierarchically-clustered heatmap of the correlation matrix. 
#Which studies seem to have similar results? 

#We can generally compare the differential expression associated with different datasets 
#or variable comparisons using a scatterplot and correlation analysis.

#The code for making scatterplots is very similar to the boxplot code that we used earlier. 
#It uses a y~x formula, and can include many of the same parameters (e.g., x and y labels, color)

plot(MetaAnalysis_FoldChanges$GSE109315_StressResilient_Vs_Control~
       MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Control, 
     ylab="Stress Resilient Log2FC", xlab="Stress Susceptible Log2FC" )

#Within this plot, each data point represents the differential expression results (log2FC)
#for a particular gene for two different comparisons: stress susceptible vs. no stress (x-axis) 
#and stress resilient vs. no stress (y axis)

#From looking at this plot, we can see that, in general, if a gene shows a positive log2FC 
#for the stress susceptible vs. no stress comparison (i.e., the stress susceptible group has 
#greater log2 expression for that gene than the no stress group), then that gene is 
#also likely to have a positive log2FC for the stress resilient vs. no stress comparison.

#Similarly, if a gene shows a negative log2FC for the stress susceptible vs. 
#no stress comparison (i.e., the stress susceptible group has lower log2 expression 
#for that gene than the no stress group), then that gene is also likely to have 
#a negative log2FC for the stress resilient vs. no stress comparison.

#This means that the differential expression results associated with the stress 
#susceptible and stress resilient comparisons are positively correlated - 
#they show a similar direction of effect.

#We can illustrate this by adding a trendline to the plot:

#We use a linear regression model to calculate the intercept and slope for the linear relationship between the two variables using the y~x formula above:
Trendline<-lm(MetaAnalysis_FoldChanges$GSE109315_StressResilient_Vs_Control~
                MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Control)

#And then add that line to our scatterplot:
abline(Trendline, col=2, lwd=3)

#If we want to know whether that linear relationship is stronger than we might expect due to random chance:
summary.lm(Trendline)
# Call:
#   lm(formula = MetaAnalysis_FoldChanges$GSE109315_StressResilient_Vs_Control ~ 
#        MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Control)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -3.13469 -0.09483 -0.00301  0.09229  2.03775 
# 
# Coefficients:
#   Estimate Std. Error t value             Pr(>|t|)    
# (Intercept)                                                     0.008518   0.001479   5.758        0.00000000858 ***
#   MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Control 0.260519   0.006695  38.913 < 0.0000000000000002 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.2683 on 32951 degrees of freedom
# (2257 observations deleted due to missingness)
# Multiple R-squared:  0.04393,	Adjusted R-squared:  0.0439 
# F-statistic:  1514 on 1 and 32951 DF,  p-value: < 0.00000000000000022

#The estimate for intercept is where the trend line crosses the y-axis
#The estimate for "Stress Susceptible vs. Control" is the slope - how much of an increase in Log2FC you should expect in the "Stress Resilient vs. Control" if you see a one unit increase in Log2FC for "Stress Susceptible vs. Control" 
#The Pr(>|t|) is the p-value for that relationship, in this case it is smaller than R is willing to display (<2e-16)

#If we want to gene a normalized correlation coefficient for this relationship (ranging between -1 to 1, with -1 being a perfect negative correlation and +1 being a perfect positive correlation), we can run a correlation analysis:
#While running this correlation analysis, we have to tell R to ignore any rows of differential expression output that don't have Log2FC for one of our variables (use "pairwise complete observations")
cor(MetaAnalysis_FoldChanges$GSE109315_StressResilient_Vs_Control,
    MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Control, use="pairwise.complete.obs")
#[1] 0.2096038

#So here's where things get particularly cool: 
#We can actually run a correlation analysis comparing the Log2FC for every set 
#of differential expression results in our meta-analysis using a single line of code.
#This is called a correlation matrix.
cor(as.matrix(MetaAnalysis_FoldChanges[,-1]), use="pairwise.complete.obs")
# GSE109315_StressResilient_Vs_Control GSE109315_StressSusceptible_Vs_Control
# GSE109315_StressResilient_Vs_Control                            1.000000000                            0.209603799
# GSE109315_StressSusceptible_Vs_Control                          0.209603799                            1.000000000
# GSE81672_StressResilient_Vs_Control                             0.041576271                            0.058634459
# GSE81672_StressSusceptible_Vs_Control                           0.053830157                           -0.022506798
# GSE85472_CSDS_vs_Control                                       -0.020812101                            0.088833896
# GSE81587_CUMS_vs_Control                                       -0.001838177                           -0.041118824
# GSE102556.1_CUMS_vs_Control                                     0.005626461                           -0.010055882
# GSE151807_CUMS_vs_Control                                      -0.004870402                            0.001522263
# GSE81672_StressResilient_Vs_Control GSE81672_StressSusceptible_Vs_Control
# GSE109315_StressResilient_Vs_Control                            0.04157627                           0.053830157
# GSE109315_StressSusceptible_Vs_Control                          0.05863446                          -0.022506798
# GSE81672_StressResilient_Vs_Control                             1.00000000                           0.455978932
# GSE81672_StressSusceptible_Vs_Control                           0.45597893                           1.000000000
# GSE85472_CSDS_vs_Control                                        0.05527107                          -0.015546804
# GSE81587_CUMS_vs_Control                                        0.01957326                           0.002982487
# GSE102556.1_CUMS_vs_Control                                     0.01466969                           0.074005675
# GSE151807_CUMS_vs_Control                                      -0.00992222                          -0.052802083
# GSE85472_CSDS_vs_Control GSE81587_CUMS_vs_Control GSE102556.1_CUMS_vs_Control
# GSE109315_StressResilient_Vs_Control               -0.020812101             -0.001838177                 0.005626461
# GSE109315_StressSusceptible_Vs_Control              0.088833896             -0.041118824                -0.010055882
# GSE81672_StressResilient_Vs_Control                 0.055271070              0.019573262                 0.014669693
# GSE81672_StressSusceptible_Vs_Control              -0.015546804              0.002982487                 0.074005675
# GSE85472_CSDS_vs_Control                            1.000000000              0.039247153                -0.005804135
# GSE81587_CUMS_vs_Control                            0.039247153              1.000000000                 0.041525439
# GSE102556.1_CUMS_vs_Control                        -0.005804135              0.041525439                 1.000000000
# GSE151807_CUMS_vs_Control                           0.094564384              0.110097305                 0.034211134
# GSE151807_CUMS_vs_Control
# GSE109315_StressResilient_Vs_Control                -0.004870402
# GSE109315_StressSusceptible_Vs_Control               0.001522263
# GSE81672_StressResilient_Vs_Control                 -0.009922220
# GSE81672_StressSusceptible_Vs_Control               -0.052802083
# GSE85472_CSDS_vs_Control                             0.094564384
# GSE81587_CUMS_vs_Control                             0.110097305
# GSE102556.1_CUMS_vs_Control                          0.034211134
# GSE151807_CUMS_vs_Control                            1.000000000
#I find that these correlation matrices can be a little easier to look at in Excel, 
#so I often output them into a file:
write.csv(cor(as.matrix(MetaAnalysis_FoldChanges[,-1]), use="pairwise.complete.obs"), 
          "CorMatrix_PFC_StressComparisons_Log2FC.csv")

#In the output, each cell includes the correlation coefficient reflecting the similarity of the 
#effects for the variable in the row and column.  
#So at the intersection of the row for "StressSusceptible_Vs_Control" and column 
#"StressResilient_Vs_Control" we can see the correlation coefficient that we calculated earlier (0.5155)
#The intersection of a variable with itself creates a coefficient of 1 (i.e., identical)

#Looking at the full matrix, most of the correlation coefficients are very close to 0.
#The only correlation coefficients that are larger, positive numbers are comparisons that 
#come from the same dataset originally. e.g., "StressSusceptible_Vs_Control" and  "StressResilient_Vs_Control"

#Disappointing, but not surprising.  The comparisons that come from the same dataset (e.g.) 
#often reflect comparisons with the same reference group. Therefore, any random variation in the reference group 
#is going to be artificially shared in the differential expression results for both comparisons.


#You can also visualize that correlation matrix using a hierarchically-clustered heatmap:
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[,-1]), use="pairwise.complete.obs"))
#In this heatmap, white/yellow indicates a more positive correlation
#The groups are placed in order by similarity, as determined by hierarchical clustering.
#The lines ("tree branches") on the left and top illustrate that similarity (clustering) using a "dendrogram"

#################################


#8) Run a meta-analysis using all of the effect sizes for each gene that has data from at least 2 studies. 


#We can only run a meta-analysis if there are differential expression results 
#from more than one comparison.
#Since I know that the differential expression results from the same study (dataset) 
#are artificially correlated, I would actually prefer that there are results from more than one dataset.

#How many genes satisfy this criteria?

#This code caculates the number of NAs in each row:
MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges, 1, function(y) sum(is.na(y)))

#I'm going to make a histogram of the results because I'm curious to see how they are distributed
hist(MetaAnalysis_FoldChanges_NAsPerRow)


#Or, since there are a limited number of integer answers (0-11), I could make a table of the results:
table(MetaAnalysis_FoldChanges_NAsPerRow)
# MetaAnalysis_FoldChanges_NAsPerRow
# 0     1     2     3     4     5     6     7 
# 14733  2104  4542  1051  2104  3176  5500  2000 


#5 NAs is too many

RunBasicMetaAnalysis<-function(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV){
  
  MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges, 1, function(y) sum(is.na(y)))
  
  print("Table of # of NAs per Row (Gene):")
  print(table(MetaAnalysis_FoldChanges_NAsPerRow))
  
  MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  
  print("MetaAnalysis_FoldChanges_ForMeta:")
  print(str(MetaAnalysis_FoldChanges_ForMeta))
  
  #I'm going to make an empty matrix to store the results of my meta-analysis:
  metaOutput<-matrix(NA, length(MetaAnalysis_FoldChanges_ForMeta$x), 6)
  
  #And then run a loop that run's a meta-analysis on the differential expression results (columns 2-10) for each gene (row):
  for(i in c(1:length(MetaAnalysis_FoldChanges_ForMeta$x))){
    
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-1])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-1])
    
    #I added a function tryCatch that double-checks that the meta-analysis function (rma) doesn't produce errors (which breaks the loop):
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(effect, var), error = function(e) {skip_to_next <<- TRUE})
    
    if(skip_to_next){}else{
      TempMeta<-rma(effect, var)
      metaOutput[i, 1]<-TempMeta$b #gives estimate Log2FC
      metaOutput[i, 2]<-TempMeta$se #gives standard error
      metaOutput[i, 3]<-TempMeta$pval #gives pval
      metaOutput[i, 4]<-TempMeta$ci.lb #gives confidence interval lower bound
      metaOutput[i, 5]<-TempMeta$ci.ub #gives confidence interval upper bound
      metaOutput[i, 6]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
      rm(TempMeta)
    }
    rm(effect, var)
  }
  
  colnames(metaOutput)<-c("Log2FC_estimate", "SE", "pval", "CI_lb", "CI_ub", "Number_Of_Comparisons")
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta$x
  
  metaOutput<<-metaOutput
  return(metaOutput)
  
  print("metaOutput:")
  print(str(metaOutput))
  
  print("Top of metaOutput:")
  print(head(metaOutput))
  
  print("Bottom of metaOutput")
  print(tail(metaOutput))
  
}

#Example Usage:
NumberOfComparisons=8
CutOffForNAs=3
#I want at least 5 datasets (so 3 na's is too many)

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

# [1] "Table of # of NAs per Row (Gene):"
# MetaAnalysis_FoldChanges_NAsPerRow
# 0     1     2     3     4     5     6     7 
# 14733  2104  4542  1051  2104  3176  5500  2000 
# [1] "MetaAnalysis_FoldChanges_ForMeta:"
# 'data.frame':	21379 obs. of  9 variables:
#   $ x                                     : chr  "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# $ GSE109315_StressResilient_Vs_Control  : num  -0.084 0.0425 -0.408 0.2268 -0.0261 ...
# $ GSE109315_StressSusceptible_Vs_Control: num  0.1929 -0.0991 -0.0613 0.3685 0.0872 ...
# $ GSE81672_StressResilient_Vs_Control   : num  -0.7165 -0.0636 0.6613 0.1336 0.0292 ...
# $ GSE81672_StressSusceptible_Vs_Control : num  -1.667 -0.0998 -1.264 -0.2959 -0.0832 ...
# $ GSE85472_CSDS_vs_Control              : num  NA -0.035 NA 0.0365 -0.116 ...
# $ GSE81587_CUMS_vs_Control              : num  -1.536 0.0853 0.1828 0.1766 -0.0161 ...
# $ GSE102556.1_CUMS_vs_Control           : num  -0.3705 0.00622 -0.2275 0.2992 0.06806 ...
# $ GSE151807_CUMS_vs_Control             : num  0.033 0.155 NA -0.125 -0.036 ...
# NULL
# There were 50 or more warnings (use warnings() to see the first 50)

str(metaOutput)
# num [1:21379, 1:6] -0.050652 -0.000666 -0.170171 0.08503 -0.047527 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:21379] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ : chr [1:6] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

head(metaOutput)
# Log2FC_estimate         SE      pval       CI_lb      CI_ub Number_Of_Comparisons
# 0610005C13Rik   -0.0506520554 0.12215410 0.6783930 -0.29006968 0.18876557                     7
# 0610009B22Rik   -0.0006663673 0.03893387 0.9863446 -0.07697535 0.07564262                     8
# 0610009E02Rik   -0.1701706846 0.19878017 0.3919564 -0.55977265 0.21943128                     6
# 0610009L18Rik    0.0850303704 0.07872635 0.2801087 -0.06927044 0.23933118                     8
# 0610010K14Rik   -0.0475266528 0.03319646 0.1522359 -0.11259052 0.01753721                     8
# 0610012G03Rik    0.0272287434 0.02131852 0.2015195 -0.01455478 0.06901227                     8

tail(metaOutput)
# Log2FC_estimate         SE        pval       CI_lb      CI_ub Number_Of_Comparisons
# Rbfox1     -0.01417218 0.06336024 0.823009220 -0.13835597 0.11001161                     6
# Rdm1       -0.01736512 0.03794070 0.647174460 -0.09172753 0.05699730                     6
# Sptbn2     -0.03186377 0.02651631 0.229492682 -0.08383479 0.02010725                     6
# Strbp       0.06556348 0.02492778 0.008535033  0.01670594 0.11442103                     6
# Syne1       0.01073551 0.03582372 0.764424320 -0.05947769 0.08094871                     6
# Usp9x       0.00698564 0.03183765 0.826327214 -0.05541500 0.06938628                     6
########################################

## Multiple Comparison corrections
#The following code applies two different types of multiple-comparison 
#corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli) 
#Meta-analysis output with adjusted p-values is then outputted along with effect size information.

#9) Correct the meta-analysis output to take into account the fact that 
#we are running the statistical calculations many times and therefore 
#have a heightened risk of false discovery (false discovery rate correction) 

colnames(metaOutput)
# [1] "Log2FC_estimate"       "SE"                    "pval"                  "CI_lb"                
# [5] "CI_ub"                 "Number_Of_Comparisons"
#The p-values are in column 3

#I'm going to grab just the column with the pvalues and run a multiple comparisons correction using the Benjamini-Hochberg method ("FDR" or "q-value")
tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))

#unfortunately, the output from that function re-orders the FDR-corrected p-values 
#in order of "significance"
#We would like the FDR-corrected p-values to be in the their original order 
#(i. the order of the rest of our statistical output!). 
#This order is recorded in the index (row numbers) for the p-values:
metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]

str(metaPvalAdj)
# num [1:21379, 1:2] 0.678 0.986 0.392 0.28 0.152 ...
# - attr(*, "dimnames")=List of 2
# ..$ : NULL
# ..$ : chr [1:2] "rawp" "BH"

#adjusted pvalue object is in same orientation as metaOutput so can simply be binded together
#We can double-check that it is in the same order by comparing the uncorrected p-values in the metaPvalAdj object to the pvalues in the metaOutput matrix

plot(metaPvalAdj[,1]~metaOutput[,3])
#Perfectly straight line - they are exactly the same! *phew*

metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])

str(metaOutputFDR)
# num [1:21379, 1:7] -0.050652 -0.000666 -0.170171 0.08503 -0.047527 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:21379] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

#Let's functionalize it!
FalseDiscoveryCorrection<-function(metaOutput){
  
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))
  
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  colnames(metaOutputFDR)[7]<-"FDR"
  
  metaOutputFDR<<-metaOutputFDR
  
  print("metaOutputFDR:")
  print(str(metaOutputFDR))
  
  write.csv(metaOutputFDR, "metaOutputFDR.csv")
  
  #a version of the output in order by p-value:
  metaOutputFDR_OrderbyPval<<-metaOutputFDR[order(metaOutputFDR[,3]),]
  
  #Let's write out a version of the output in order by p-value:
  write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_orderedByPval.csv")
  
  print("Do we have any genes that are statistically significant following false discovery rate correction?")
  print(sum(metaOutputFDR[,7]<0.10, na.rm=TRUE))
  
  print("What are the top results?")
  print(head(metaOutputFDR[order(metaOutputFDR[,3]),]))
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
}

#Example usage:

FalseDiscoveryCorrection(metaOutput)
# [1] "metaOutputFDR:"
# num [1:21379, 1:7] -0.050652 -0.000666 -0.170171 0.08503 -0.047527 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:21379] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following false discovery rate correction?"
# [1] 237
# [1] "What are the top results?"
# Log2FC_estimate         SE              pval       CI_lb       CI_ub Number_Of_Comparisons           FDR
# Hepacam      -0.1090172 0.01806603 0.000000001595820 -0.14442595 -0.07360840                     8 0.00002802899
# Elovl5       -0.1075961 0.01816007 0.000000003125414 -0.14318913 -0.07200298                     8 0.00002802899
# Dusp1        -0.5071662 0.08615048 0.000000003933157 -0.67601803 -0.33831435                     8 0.00002802899
# Mid1ip1      -0.1138308 0.01980716 0.000000009086755 -0.15265208 -0.07500944                     8 0.00004856644
# Egfem1        0.1361996 0.02390941 0.000000012230138  0.08933805  0.18306121                     7 0.00005057839
# Csrp1        -0.1576198 0.02791261 0.000000016336256 -0.21232756 -0.10291214                     8 0.00005057839
# ############################

#10) Determine which are the top differentially expressed genes and create forest 
#plots to visualize the effect sizes for those top differentially expressed genes 
#across the different studies. 

row.names(metaOutputFDR_OrderbyPval)[c(1:100)]
# [1] "Hepacam"       "Elovl5"        "Dusp1"         "Mid1ip1"       "Egfem1"        "Csrp1"         "Fbxo32"        "Cyp39a1"      
# [9] "Scd2"          "Tcf7l2"        "Tnni1"         "Fads6"         "P2rx7"         "LOC108167320"  "Syndig1"       "Gna12"        
# [17] "Otx2"          "Ppp1r16b"      "Rbms1"         "Xdh"           "Aldh3b1"       "Dync1i1"       "Nsl1"          "Fa2h"         
# [25] "Rcan2"         "Best2"         "Kif1c"         "Unc5b"         "Btg2"          "Olfml1"        "Glmp"          "Zfp979"       
# [33] "Adgra3"        "S100b"         "Fam163a"       "Sox9"          "Slc25a15"      "Siae"          "Atf4"          "Cpn1"         
# [41] "Map3k11"       "Col6a2"        "Tshr"          "Ggnbp2"        "9430085M18Rik" "Mxra7"         "4930455C13Rik" "Gm14434"      
# [49] "Icmt"          "Stpg4"         "Ccn1"          "Myo1h"         "Gm30539"       "Gm17762"       "Cfap276"       "Prr18"        
# [57] "Cd274"         "Gm41875"       "Ube2d1"        "Bfsp2"         "Csrnp1"        "Psrc1"         "Kdr"           "Nipal4"       
# [65] "Tek"           "Stk31"         "Desi1"         "Gm15850"       "Ddr1"          "Rnd2"          "Dnai3"         "Slc39a8"      
# [73] "Lpar1"         "Gm973"         "Slc24a4"       "Gm46405"       "Tedc2"         "Aqp4"          "Ier2"          "Gm34421"      
# [81] "Alpl"          "Zfp362"        "Gm40468"       "5430421F17Rik" "Lgi3"          "Gdf11"         "Chrd"          "Trp53inp2"    
# [89] "Ubqln1"        "Sft2d2"        "Klf7"          "Rgs8"          "Gm30716"       "Dancr"         "Jsrp1"         "BC034090"     
# [97] "Arc"           "Ppp1r21"       "Mmgt2"         "Anxa5"         

#Let's plot some of those top results!

#Quickly looking at the range of Log2FC values to figure out the limits for the x-axis for the forest plots:
hist(metaOutputFDR[,1], breaks=40)
#Range is mostly -0.5 to 0.5, but there are a few with Log2FC as big as -1 or 1.5

MakeForestPlots<-function(GeneSymbol){
  
  pdf(paste("ForestPlot_", GeneSymbol, ".pdf", sep=""), height=5, width=8)
  
  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
  var<-as.numeric(MetaAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
  
  forest.rma(rma(effect, var),slab=colnames(MetaAnalysis_FoldChanges_ForMeta)[-1],  xlim=c(-3, 3))
  
  mtext(paste(GeneSymbol), line=-1.5, cex=2)
  dev.off()
}

setwd("/labs/twc/JinglinXiong/BAP_Updated_2025/Output")

MakeForestPlots("Hepacam") 
MakeForestPlots("Elovl5") 
MakeForestPlots("Dusp1")
MakeForestPlots("Mid1ip1") 
MakeForestPlots("Egfem1") 
MakeForestPlots("Csrp1")
MakeForestPlots("Fbxo32") 
MakeForestPlots("Cyp39a1")
MakeForestPlots("Scd2")
MakeForestPlots("Tcf7l2") 
MakeForestPlots("Tnni1") 
MakeForestPlots("Fads6") 
MakeForestPlots("P2rx7") 
MakeForestPlots("LOC108167320") 
MakeForestPlots("Syndig1") 
MakeForestPlots("Gna12") 

#In general, it seems to me like the results that have both a significant FDR 
#& large effect size tend to be the most convincing 

#Here's a summary of the distribution for the Log2FCs
summary(metaOutputFDR[,1])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -1.49021 -0.03419 -0.00058 -0.00135  0.03142  1.54745       83 
# how many are statistically significant
sum(metaOutputFDR[,7]<0.05, na.rm = TRUE)
# [1] 133
# store significant data
metaOutputSig <- metaOutputFDR[metaOutputFDR[,7]<0.05,]
# remove na
metaOutputSigNoNa <- na.omit(metaOutputSig)

# count upregulated genes
metaOutputupregulated <- metaOutputSigNoNa[metaOutputSigNoNa[,1]>0,]

str(metaOutputupregulated) #36
# num [1:36, 1:7] 0.105 0.211 0.168 0.212 0.182 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:36] "4732471J01Rik" "4930455C13Rik" "5430421F17Rik" "Best2" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
metaOutputdownregulated <- metaOutputSigNoNa[metaOutputSigNoNa[,1]<0,]
str(metaOutputdownregulated)
# num [1:97, 1:7] -0.7164 -0.167 -0.1109 -0.1332 -0.0885 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:97] "1700017M07Rik" "9430085M18Rik" "Adgra3" "Aldh3b1" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

#Let's see how many results are both statistically significant 
#(using a more conservative FDR) and have a larger estimated Log2FC (>0.10 or <-0.10)
sum(metaOutputFDR[,7]<0.05 & abs(metaOutputFDR[,1])>0.1, na.rm=TRUE)
#[1] 93

#What are their gene symbols?
row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & 
                                       abs(metaOutputFDR_OrderbyPval[,1])>0.1 & 
                                       is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]
# [1] "Hepacam"       "Elovl5"        "Dusp1"         "Mid1ip1"       "Egfem1"        "Csrp1"         "Fbxo32"        "Cyp39a1"      
# [9] "Tcf7l2"        "Tnni1"         "Fads6"         "LOC108167320"  "Syndig1"       "Gna12"         "Otx2"          "Ppp1r16b"     
# [17] "Rbms1"         "Xdh"           "Aldh3b1"       "Nsl1"          "Fa2h"          "Rcan2"         "Best2"         "Unc5b"        
# [25] "Btg2"          "Olfml1"        "Zfp979"        "Adgra3"        "S100b"         "Fam163a"       "Sox9"          "Siae"         
# [33] "Atf4"          "Cpn1"          "Col6a2"        "Tshr"          "9430085M18Rik" "Mxra7"         "4930455C13Rik" "Gm14434"      
# [41] "Ccn1"          "Myo1h"         "Gm30539"       "Gm17762"       "Cfap276"       "Prr18"         "Cd274"         "Gm41875"      
# [49] "Bfsp2"         "Csrnp1"        "Psrc1"         "Kdr"           "Nipal4"        "Tek"           "Desi1"         "Gm15850"      
# [57] "Ddr1"          "Rnd2"          "Dnai3"         "Slc39a8"       "Lpar1"         "Gm973"         "Slc24a4"       "Gm46405"      
# [65] "Aqp4"          "Ier2"          "Gm34421"       "Gm40468"       "5430421F17Rik" "Lgi3"          "Gdf11"         "Chrd"         
# [73] "Gm30716"       "Dancr"         "Jsrp1"         "BC034090"      "Arc"           "Mmgt2"         "Anxa5"         "Egr1"         
# [81] "4732471J01Rik" "Qprt"          "Atad2"         "Ngf"           "Nphs1os"       "Nrros"         "Cox4i2"        "Sntb1"        
# [89] "1700017M07Rik" "S1pr5"         "LOC102640693"  "Fos"           "Glt28d2" 



#Instead of effect sizes, that if I focus on genes represented in more than 5 datasets?

colnames(metaOutputFDR)
# [1] "Log2FC_estimate"       "SE"                    "pval"                  "CI_lb"                 "CI_ub"                
# [6] "Number_Of_Comparisons" "FDR"  

sum(metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>5, na.rm=TRUE)
#[1] 133

row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & 
                                       metaOutputFDR_OrderbyPval[,6]>5 & 
                                       is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]

#Adding in a p-value threshold too:
sum(metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>5 & metaOutputFDR[,3]<0.00001, na.rm=TRUE)
#[1] 42

row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & 
                                       metaOutputFDR_OrderbyPval[,6]>5 & 
                                       metaOutputFDR_OrderbyPval[,3]<0.00001 & 
                                        is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]
# [1] "Hepacam"      "Elovl5"       "Dusp1"        "Mid1ip1"      "Egfem1"       "Csrp1"        "Fbxo32"       "Cyp39a1"     
# [9] "Scd2"         "Tcf7l2"       "Tnni1"        "Fads6"        "P2rx7"        "LOC108167320" "Syndig1"      "Gna12"       
# [17] "Otx2"         "Ppp1r16b"     "Rbms1"        "Xdh"          "Aldh3b1"      "Dync1i1"      "Nsl1"         "Fa2h"        
# [25] "Rcan2"        "Best2"        "Kif1c"        "Unc5b"        "Btg2"         "Olfml1"       "Glmp"         "Zfp979"      
# [33] "Adgra3"       "S100b"        "Fam163a"      "Sox9"         "Slc25a15"     "Siae"         "Atf4"         "Cpn1"        
# [41] "Map3k11"      "Col6a2"  



#What if we look at the correlation between datasets again, 
#but only using the top genes found in all 7 datasets?
sum(metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>6, na.rm=TRUE)
#[1] 119

TopGenesInAll11Datasets<-row.names(metaOutputFDR)[metaOutputFDR[,7]<0.05 & 
                                                    metaOutputFDR[,6]>6 & 
                                                    is.na(metaOutputFDR[,6])==FALSE]

cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges$x%in%TopGenesInAll11Datasets),-1]), use="pairwise.complete.obs")
# GSE109315_StressResilient_Vs_Control GSE109315_StressSusceptible_Vs_Control
# GSE109315_StressResilient_Vs_Control                            1.000000000                           -0.002076656
# GSE109315_StressSusceptible_Vs_Control                         -0.002076656                            1.000000000
# GSE81672_StressResilient_Vs_Control                             0.209615728                            0.261492475
# GSE81672_StressSusceptible_Vs_Control                           0.532757515                            0.080748422
# GSE85472_CSDS_vs_Control                                        0.293433600                            0.441771625
# GSE81587_CUMS_vs_Control                                        0.387986926                            0.126548566
# GSE102556.1_CUMS_vs_Control                                     0.353700938                            0.067179672
# GSE151807_CUMS_vs_Control                                       0.462142021                            0.218890674
# GSE81672_StressResilient_Vs_Control GSE81672_StressSusceptible_Vs_Control
# GSE109315_StressResilient_Vs_Control                             0.2096157                            0.53275752
# GSE109315_StressSusceptible_Vs_Control                           0.2614925                            0.08074842
# GSE81672_StressResilient_Vs_Control                              1.0000000                            0.40499018
# GSE81672_StressSusceptible_Vs_Control                            0.4049902                            1.00000000
# GSE85472_CSDS_vs_Control                                         0.5810925                            0.47963286
# GSE81587_CUMS_vs_Control                                         0.1469892                            0.49890084
# GSE102556.1_CUMS_vs_Control                                      0.5763801                            0.56627224
# GSE151807_CUMS_vs_Control                                        0.6270959                            0.64181127
# GSE85472_CSDS_vs_Control GSE81587_CUMS_vs_Control GSE102556.1_CUMS_vs_Control
# GSE109315_StressResilient_Vs_Control                  0.2934336                0.3879869                  0.35370094
# GSE109315_StressSusceptible_Vs_Control                0.4417716                0.1265486                  0.06717967
# GSE81672_StressResilient_Vs_Control                   0.5810925                0.1469892                  0.57638005
# GSE81672_StressSusceptible_Vs_Control                 0.4796329                0.4989008                  0.56627224
# GSE85472_CSDS_vs_Control                              1.0000000                0.5148619                  0.31821068
# GSE81587_CUMS_vs_Control                              0.5148619                1.0000000                  0.27730371
# GSE102556.1_CUMS_vs_Control                           0.3182107                0.2773037                  1.00000000
# GSE151807_CUMS_vs_Control                             0.7537755                0.5784865                  0.57439787
# GSE151807_CUMS_vs_Control
# GSE109315_StressResilient_Vs_Control                   0.4621420
# GSE109315_StressSusceptible_Vs_Control                 0.2188907
# GSE81672_StressResilient_Vs_Control                    0.6270959
# GSE81672_StressSusceptible_Vs_Control                  0.6418113
# GSE85472_CSDS_vs_Control                               0.7537755
# GSE81587_CUMS_vs_Control                               0.5784865
# GSE102556.1_CUMS_vs_Control                            0.5743979
# GSE151807_CUMS_vs_Control                              1.0000000

pdf("Heatmap_CorMatrix_PFC_StressDatasets_GenesFDR05inAll.pdf", height=10, width=10)
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges$x%in%TopGenesInAll11Datasets),
                                               -1]), use="pairwise.complete.obs"))
dev.off()

CorMatrixTopGenes<-cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges
                                                                $x%in%TopGenesInAll11Datasets),-1]),
                       use="pairwise.complete.obs")

write.csv(CorMatrixTopGenes, "CorMatrix_PFC_StressDatasets_GenesFDR05inAll.csv")

