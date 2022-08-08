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

library(plyr)
library(metafor)
library(reshape)
library(multtest)

#Set working directory:
setwd("~/Documents/2022 Summer/Gemma data/PFC-social defeat")

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

ListOfDEResults<-list(DEResults_GSE109315, DEResults_GSE114224_Left, DEResults_GSE114224_Right, DEResults_GSE81672, DEResults_GSE85472)

#I just discovered that if I have datasets with *comparisons with the same name*, when I join the datasets some of those comparisons disappear. :(

AligningDatasets(ListOfDEResults)
# [1] "MetaAnalysis_FoldChange_Dfs:"
# List of 5
# $ :'data.frame':	32953 obs. of  3 variables:
#   ..$ x                                     : chr [1:32953] "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# ..$ GSE109315_StressResilient_Vs_Control  : num [1:32953] -0.084 -0.0858 0.0425 -0.408 0.2268 ...
# ..$ GSE109315_StressSusceptible_Vs_Control: num [1:32953] 0.1929 -0.0733 -0.0991 -0.0613 0.3685 ...
# $ :'data.frame':	15951 obs. of  2 variables:
#   ..$ x                             : chr [1:15951] "0610009B22Rik" "0610010K14Rik" "0610012G03Rik" "0610040J01Rik" ...
# ..$ GSE114224_Left_CSDS_vs_Control: num [1:15951] -0.0093 0.02176 -0.01592 0.00432 -0.021 ...
# $ :'data.frame':	15949 obs. of  2 variables:
#   ..$ x                              : chr [1:15949] "0610009B22Rik" "0610010K14Rik" "0610012G03Rik" "0610040J01Rik" ...
# ..$ GSE114224_Right_CSDS_vs_Control: num [1:15949] 0.002769 0.02463 0.02097 0.000499 -0.003534 ...
# $ :'data.frame':	21628 obs. of  3 variables:
#   ..$ x                                    : chr [1:21628] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ GSE81672_StressResilient_Vs_Control  : num [1:21628] -0.7165 -0.0636 0.6613 0.1336 0.0292 ...
# ..$ GSE81672_StressSusceptible_Vs_Control: num [1:21628] -1.667 -0.0998 -1.264 -0.2959 -0.0832 ...
# $ :'data.frame':	20237 obs. of  2 variables:
#   ..$ x                       : chr [1:20237] "0610009B22Rik" "0610009L18Rik" "0610010K14Rik" "0610012G03Rik" ...
# ..$ GSE85472_CSDS_vs_Control: num [1:20237] -0.03498 0.03654 -0.116 0.01691 0.00807 ...
# NULL
# [1] "MetaAnalysis_FoldChanges:"
# 'data.frame':	33841 obs. of  8 variables:
#   $ x                                     : chr  "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# $ GSE109315_StressResilient_Vs_Control  : num  -0.084 -0.0858 0.0425 -0.408 0.2268 ...
# $ GSE109315_StressSusceptible_Vs_Control: num  0.1929 -0.0733 -0.0991 -0.0613 0.3685 ...
# $ GSE114224_Left_CSDS_vs_Control        : num  NA NA -0.0093 NA NA ...
# $ GSE114224_Right_CSDS_vs_Control       : num  NA NA 0.00277 NA NA ...
# $ GSE81672_StressResilient_Vs_Control   : num  -0.7165 NA -0.0636 0.6613 0.1336 ...
# $ GSE81672_StressSusceptible_Vs_Control : num  -1.667 NA -0.0998 -1.264 -0.2959 ...
# $ GSE85472_CSDS_vs_Control              : num  NA NA -0.035 NA 0.0365 ...
# NULL
# [1] "MetaAnalysis_SV_Dfs:"
# List of 5
# $ :'data.frame':	32953 obs. of  3 variables:
#   ..$ x                                     : chr [1:32953] "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# ..$ GSE109315_StressResilient_Vs_Control  : num [1:32953] 0.1704 0.0202 0.0063 0.1143 0.0344 ...
# ..$ GSE109315_StressSusceptible_Vs_Control: num [1:32953] 0.09173 0.01089 0.00339 0.06156 0.01856 ...
# $ :'data.frame':	15951 obs. of  2 variables:
#   ..$ x                             : chr [1:15951] "0610009B22Rik" "0610010K14Rik" "0610012G03Rik" "0610040J01Rik" ...
# ..$ GSE114224_Left_CSDS_vs_Control: num [1:15951] 1.91e-04 6.80e-05 1.20e-04 4.90e-05 6.32e-05 ...
# $ :'data.frame':	15949 obs. of  2 variables:
#   ..$ x                              : chr [1:15949] "0610009B22Rik" "0610010K14Rik" "0610012G03Rik" "0610040J01Rik" ...
# ..$ GSE114224_Right_CSDS_vs_Control: num [1:15949] 1.24e-04 1.30e-04 1.78e-04 1.63e-04 7.44e-05 ...
# $ :'data.frame':	21628 obs. of  3 variables:
#   ..$ x                                    : chr [1:21628] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ GSE81672_StressResilient_Vs_Control  : num [1:21628] 0.77289 0.00638 0.23854 0.04249 0.011 ...
# ..$ GSE81672_StressSusceptible_Vs_Control: num [1:21628] 0.91575 0.00756 0.28301 0.0504 0.01303 ...
# $ :'data.frame':	20237 obs. of  2 variables:
#   ..$ x                       : chr [1:20237] "0610009B22Rik" "0610009L18Rik" "0610010K14Rik" "0610012G03Rik" ...
# ..$ GSE85472_CSDS_vs_Control: num [1:20237] 0.00185 0.00309 0.00189 0.00142 0.0056 ...
# NULL
# [1] "MetaAnalysis_SV:"
# 'data.frame':	33841 obs. of  8 variables:
#   $ x                                     : chr  "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# $ GSE109315_StressResilient_Vs_Control  : num  0.1704 0.0202 0.0063 0.1143 0.0344 ...
# $ GSE109315_StressSusceptible_Vs_Control: num  0.09173 0.01089 0.00339 0.06156 0.01856 ...
# $ GSE114224_Left_CSDS_vs_Control        : num  NA NA 0.000191 NA NA ...
# $ GSE114224_Right_CSDS_vs_Control       : num  NA NA 0.000124 NA NA ...
# $ GSE81672_StressResilient_Vs_Control   : num  0.77289 NA 0.00638 0.23854 0.04249 ...
# $ GSE81672_StressSusceptible_Vs_Control : num  0.91575 NA 0.00756 0.28301 0.0504 ...
# $ GSE85472_CSDS_vs_Control              : num  NA NA 0.00185 NA 0.00309 ...
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
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                     0.008518   0.001479   5.758 8.58e-09 ***
#   MetaAnalysis_FoldChanges$GSE109315_StressSusceptible_Vs_Control 0.260519   0.006695  38.913  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2683 on 32951 degrees of freedom
# (888 observations deleted due to missingness)
# Multiple R-squared:  0.04393,	Adjusted R-squared:  0.0439 
# F-statistic:  1514 on 1 and 32951 DF,  p-value: < 2.2e-16

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
#                                             GSE109315_StressResilient_Vs_Control
# GSE109315_StressResilient_Vs_Control                            1.000000000
# GSE109315_StressSusceptible_Vs_Control                          0.209603799
# GSE114224_Left_CSDS_vs_Control                                 -0.009182396
# GSE114224_Right_CSDS_vs_Control                                 0.062960251
# GSE81672_StressResilient_Vs_Control                             0.041576271
# GSE81672_StressSusceptible_Vs_Control                           0.053830157
# GSE85472_CSDS_vs_Control                                       -0.020812101
#                                               GSE109315_StressSusceptible_Vs_Control GSE114224_Left_CSDS_vs_Control
# GSE109315_StressResilient_Vs_Control                               0.20960380                   -0.009182396
# GSE109315_StressSusceptible_Vs_Control                             1.00000000                    0.089912682
# GSE114224_Left_CSDS_vs_Control                                     0.08991268                    1.000000000
# GSE114224_Right_CSDS_vs_Control                                   -0.02448630                    0.219358022
# GSE81672_StressResilient_Vs_Control                                0.05863446                    0.099781605
# GSE81672_StressSusceptible_Vs_Control                             -0.02250680                   -0.011179858
# GSE85472_CSDS_vs_Control                                           0.08883390                    0.094358420
#                                                      GSE114224_Right_CSDS_vs_Control GSE81672_StressResilient_Vs_Control
# GSE109315_StressResilient_Vs_Control                        0.06296025                          0.04157627
# GSE109315_StressSusceptible_Vs_Control                     -0.02448630                          0.05863446
# GSE114224_Left_CSDS_vs_Control                              0.21935802                          0.09978160
# GSE114224_Right_CSDS_vs_Control                             1.00000000                          0.05037023
# GSE81672_StressResilient_Vs_Control                         0.05037023                          1.00000000
# GSE81672_StressSusceptible_Vs_Control                       0.06718874                          0.45597893
# GSE85472_CSDS_vs_Control                                   -0.05733548                          0.05527107
#                                                 GSE81672_StressSusceptible_Vs_Control GSE85472_CSDS_vs_Control
# GSE109315_StressResilient_Vs_Control                              0.05383016              -0.02081210
# GSE109315_StressSusceptible_Vs_Control                           -0.02250680               0.08883390
# GSE114224_Left_CSDS_vs_Control                                   -0.01117986               0.09435842
# GSE114224_Right_CSDS_vs_Control                                   0.06718874              -0.05733548
# GSE81672_StressResilient_Vs_Control                               0.45597893               0.05527107
# GSE81672_StressSusceptible_Vs_Control                             1.00000000              -0.01554680
# GSE85472_CSDS_vs_Control                                         -0.01554680               1.00000000

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
# 0     1     2     3     4     5     6 
# 12786   794  4842  5699  1213  7881   626 


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
NumberOfComparisons=7
CutOffForNAs=5
#I want at least 7 datasets (so 5 na's is too many)

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

# [1] "Table of # of NAs per Row (Gene):"
# MetaAnalysis_FoldChanges_NAsPerRow
# 0     1     2     3     4     5     6 
# 12786   794  4842  5699  1213  7881   626 
# [1] "MetaAnalysis_FoldChanges_ForMeta:"
# 'data.frame':	25334 obs. of  8 variables:
#   $ x                                     : chr  "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# $ GSE109315_StressResilient_Vs_Control  : num  -0.084 0.0425 -0.408 0.2268 -0.0261 ...
# $ GSE109315_StressSusceptible_Vs_Control: num  0.1929 -0.0991 -0.0613 0.3685 0.0872 ...
# $ GSE114224_Left_CSDS_vs_Control        : num  NA -0.0093 NA NA 0.0218 ...
# $ GSE114224_Right_CSDS_vs_Control       : num  NA 0.00277 NA NA 0.02463 ...
# $ GSE81672_StressResilient_Vs_Control   : num  -0.7165 -0.0636 0.6613 0.1336 0.0292 ...
# $ GSE81672_StressSusceptible_Vs_Control : num  -1.667 -0.0998 -1.264 -0.2959 -0.0832 ...
# $ GSE85472_CSDS_vs_Control              : num  NA -0.035 NA 0.0365 -0.116 ...
# NULL
# There were 50 or more warnings (use warnings() to see the first 50)

str(metaOutput)
# num [1:25334, 1:6] -0.07545 -0.00625 -0.23821 0.11246 0.0194 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:25334] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ : chr [1:6] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

head(metaOutput)
#                  Log2FC_estimate          SE        pval        CI_lb      CI_ub Number_Of_Comparisons
# 0610005C13Rik    -0.075454183 0.239525899 0.752750480 -0.544916319 0.39400795                     4
# 0610009B22Rik    -0.006248948 0.008292547 0.451112289 -0.022502041 0.01000414                     7
# 0610009E02Rik    -0.238211306 0.330250337 0.470721720 -0.885490072 0.40906746                     4
# 0610009L18Rik     0.112463747 0.094580203 0.234406858 -0.072910044 0.29783754                     5
# 0610010K14Rik     0.019400199 0.006558148 0.003094588  0.006546465 0.03225393                     7
# 0610012G03Rik     0.002279618 0.013155172 0.862425986 -0.023504045 0.02806328                     7

tail(metaOutput)
#            Log2FC_estimate         SE        pval         CI_lb      CI_ub Number_Of_Comparisons
# Cltc      -0.002639249 0.02094763 0.899737712 -0.0436958547 0.03841736                     3
# Gm7467     0.093811796 0.03631940 0.009795486  0.0226270753 0.16499652                     3
# Macf1      0.037881310 0.03188043 0.234742473 -0.0246031837 0.10036580                     3
# Malat1     0.145622767 0.15662612 0.352501756 -0.1613587864 0.45260432                     3
# Rn7sk      0.050970116 0.02809653 0.069661275 -0.0040980795 0.10603831                     3
# Sptan1     0.037420419 0.01903013 0.049255022  0.0001220497 0.07471879                     3

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
# num [1:25334, 1:2] 0.75275 0.45111 0.47072 0.23441 0.00309 ...
# - attr(*, "dimnames")=List of 2
# ..$ : NULL
# ..$ : chr [1:2] "rawp" "BH"

#adjusted pvalue object is in same orientation as metaOutput so can simply be binded together
#We can double-check that it is in the same order by comparing the uncorrected p-values in the metaPvalAdj object to the pvalues in the metaOutput matrix

plot(metaPvalAdj[,1]~metaOutput[,3])
#Perfectly straight line - they are exactly the same! *phew*

metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])

str(metaOutputFDR)
# num [1:25334, 1:7] -0.07545 -0.00263 -0.23821 0.11246 -0.02447 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:25334] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
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
# num [1:25334, 1:7] -0.07545 -0.00625 -0.23821 0.11246 0.0194 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:25334] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following false discovery rate correction?"
# [1] 695
# [1] "What are the top results?"
# Log2FC_estimate          SE         pval       CI_lb       CI_ub Number_Of_Comparisons          FDR
# Col5a1       0.04411235 0.003047610 1.758368e-47  0.03813914  0.05008556                     7 4.454649e-43
# Nars        -0.03133029 0.002948264 2.239181e-26 -0.03710878 -0.02555180                     7 2.836371e-22
# Dus4l        0.02993351 0.003049209 9.530948e-23  0.02395717  0.03590985                     7 8.048568e-19
# Shisa7       0.03189163 0.003702804 7.127634e-18  0.02463427  0.03914900                     7 4.514287e-14
# Zmat3       -0.02480211 0.002962386 5.648426e-17 -0.03060828 -0.01899594                     7 2.861945e-13
# Gpr37l1     -0.03552813 0.004378199 4.866597e-16 -0.04410924 -0.02694702                     7 2.054839e-12

############################

#10) Determine which are the top differentially expressed genes and create forest 
#plots to visualize the effect sizes for those top differentially expressed genes 
#across the different studies. 

row.names(metaOutputFDR_OrderbyPval)[c(1:100)]
# [1] "Col5a1"        "Nars"          "Dus4l"         "Shisa7"        "Zmat3"         "Gpr37l1"      
# [7] "Tmem267"       "Extl1"         "Tceal3"        "Stk32c"        "Ramp3"         "Prr7"         
# [13] "Ahi1"          "Usp32"         "Znrf1"         "Asb11"         "Rgs10"         "Npy"          
# [19] "Tdrd7"         "Pou1f1"        "Npas4"         "P2ry12"        "Kcnip4"        "Fam214a"      
# [25] "Actr1a"        "Clstn1"        "Acsbg1"        "Rad23b"        "Cmpk2"         "Pcdhb20"      
# [31] "Lor"           "Klhl7"         "Lrrtm3"        "Ezr"           "Tas2r129"      "Map10"        
# [37] "Zfp36l1"       "Ptk2"          "Stard3nl"      "Slc7a10"       "Lonrf2"        "Exosc7"       
# [43] "Med20"         "Stmn4"         "Arc"           "Rwdd4a"        "Adam9"         "Lrpap1"       
# [49] "Zfp397"        "Slitrk4"       "Med6"          "Zfx"           "Nisch"         "Hax1"         
# [55] "S100b"         "Hdgfl2"        "Ccdc85b"       "Spink10"       "Slc39a12"      "Slc30a4"      
# [61] "Pex26"         "Hsd11b1"       "Dgkg"          "Kmt5c"         "Srcin1"        "Myo6"         
# [67] "Opn3"          "Rora"          "Clpp"          "Ube2e1"        "Cdh13"         "Tceal6"       
# [73] "Nudt18"        "Atp13a2"       "Rps4l"         "Crym"          "Thy1"          "Ankrd39"      
# [79] "Pycr2"         "Casp8"         "Selenon"       "Egfem1"        "Tsc22d2"       "Hnrnpa2b1"    
# [85] "Amz2"          "2410018L13Rik" "Tubg2"         "Pkia"          "Lsm8"          "Ptdss1"       
# [91] "Nxnl2"         "Anapc4"        "Ssx2ip"        "Cyth3"         "Trmt44"        "Gm5617"       
# [97] "Dbpht2"        "H4c8"          "Smg6"          "Hsph1"       

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


MakeForestPlots("Col5a1") 
MakeForestPlots("Nars") 
MakeForestPlots("Dus4l") #This one looks strange
# Warning message:
#   Fisher scoring algorithm may have gotten stuck at a local maximum.
# Setting tau^2 = 0. Check the profile likelihood plot with profile(). 
MakeForestPlots("Shisa7")
MakeForestPlots("Zmat3") # all over the place
MakeForestPlots("Gpr37l1")
MakeForestPlots("Tmem267") #all over the place
MakeForestPlots("Extl1")
MakeForestPlots("Tceal3")
# Warning message:
#   Studies with NAs omitted from model fitting. 



#In general, it seems to me like the results that have both a significant FDR 
#& large effect size tend to be the most convincing 

#Here's a summary of the distribution for the Log2FCs
summary(metaOutputFDR[,1])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
# -1.027083 -0.010951  0.000305  0.001064  0.010438  1.685594         8 

#Let's see how many results are both statistically significant 
#(using a more conservative FDR) and have a larger estimated Log2FC (>0.10 or <-0.10)
sum(metaOutputFDR[,7]<0.05 & abs(metaOutputFDR[,1])>0.1, na.rm=TRUE)
#[1] 66

#What are their gene symbols?
row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & 
                                       abs(metaOutputFDR_OrderbyPval[,1])>0.1 & 
                                       is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]
# [1] "Egfem1"        "2410018L13Rik" "Gm5617"        "Syndig1"       "Gm30539"       "4930429B21Rik"
# [7] "Kcnj10"        "Gm14124"       "F930017D23Rik" "Gm46920"       "A830052D11Rik" "Zfp862-ps"    
# [13] "Gm46405"       "Gm973"         "Tmem243"       "B3glct"        "D330050G23Rik" "Gm14199"      
# [19] "Gm38418"       "4732471J01Rik" "Gm31465"       "Gcnt4"         "4930480E11Rik" "Gm15870"      
# [25] "Gm30392"       "Gm26782"       "Gm34140"       "Gm34922"       "Gm36163"       "Cd99"         
# [31] "2900052N01Rik" "1700017M07Rik" "AU015836"      "B930018H19Rik" "Gm30716"       "Adgra3"       
# [37] "Gm30687"       "5033406O09Rik" "Ccdc146"       "Gm35572"       "Gm26994"       "LOC108168680" 
# [43] "Gm46519"       "Gm39805"       "B230217O12Rik" "Car10"         "Sycp2l"        "Siae"         
# [49] "Gm26555"       "Gm39653"       "Wt1os"         "LOC102640693"  "Mirt2"         "Gm38983"      
# [55] "Gm34004"       "Gm29994"       "Morf4l1b"      "Spinkl"        "Muc21"         "Colec11"      
# [61] "Gm16169"       "Gm6525"        "Gtpbp10"       "Gm35081"       "Gm6600"        "Gm20385"    

MakeForestPlots("Egfem1")
# Warning message:
#   Studies with NAs omitted from model fitting. 
MakeForestPlots("2410018L13Rik")
# Warning message:
#   Studies with NAs omitted from model fitting. 
MakeForestPlots("Gm5617") 
# Warning message:
#   Studies with NAs omitted from model fitting. 
MakeForestPlots("Syndig1") 
# Warning message:
#   Studies with NAs omitted from model fitting. 

#They all look good

#Instead of effect sizes, that if I focus on genes represented in more than 5 datasets?

colnames(metaOutputFDR)

sum(metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>5, na.rm=TRUE)
#[1] 376

row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & 
                                       metaOutputFDR_OrderbyPval[,6]>5 & 
                                       is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]

#Adding in a p-value threshold too:
sum(metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>5 & metaOutputFDR[,3]<0.00001, na.rm=TRUE)
#[1] 114

row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & 
                                       metaOutputFDR_OrderbyPval[,6]>5 & 
                                       metaOutputFDR_OrderbyPval[,3]<0.00001 & 
                                       is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]
# [1] "Col5a1"    "Nars"      "Dus4l"     "Shisa7"    "Zmat3"     "Gpr37l1"   "Tmem267"   "Extl1"    
# [9] "Tceal3"    "Stk32c"    "Ramp3"     "Prr7"      "Ahi1"      "Usp32"     "Znrf1"     "Asb11"    
# [17] "Rgs10"     "Npy"       "Tdrd7"     "Npas4"     "P2ry12"    "Kcnip4"    "Fam214a"   "Actr1a"   
# [25] "Acsbg1"    "Rad23b"    "Cmpk2"     "Pcdhb20"   "Lor"       "Klhl7"     "Lrrtm3"    "Ezr"      
# [33] "Map10"     "Zfp36l1"   "Ptk2"      "Stard3nl"  "Slc7a10"   "Lonrf2"    "Exosc7"    "Med20"    
# [41] "Stmn4"     "Arc"       "Rwdd4a"    "Adam9"     "Lrpap1"    "Zfp397"    "Slitrk4"   "Med6"     
# [49] "Zfx"       "Nisch"     "Hax1"      "S100b"     "Hdgfl2"    "Ccdc85b"   "Spink10"   "Slc39a12" 
# [57] "Slc30a4"   "Pex26"     "Hsd11b1"   "Dgkg"      "Kmt5c"     "Srcin1"    "Myo6"      "Opn3"     
# [65] "Rora"      "Clpp"      "Ube2e1"    "Cdh13"     "Tceal6"    "Nudt18"    "Atp13a2"   "Rps4l"    
# [73] "Crym"      "Thy1"      "Ankrd39"   "Pycr2"     "Casp8"     "Selenon"   "Tsc22d2"   "Hnrnpa2b1"
# [81] "Amz2"      "Tubg2"     "Pkia"      "Lsm8"      "Ptdss1"    "Nxnl2"     "Anapc4"    "Ssx2ip"   
# [89] "Cyth3"     "Trmt44"    "Dbpht2"    "H4c8"      "Smg6"      "Hsph1"     "Ndrg2"     "Tpst1"    
# [97] "Lamtor4"   "Sfr1"      "Top2b"     "Ptges3"    "Scrg1"     "Ndfip1"    "Pdcd7"     "Prdm6"    
# [105] "Rab5c"     "Hikeshi"   "Ranbp1"    "Lpp"       "Lmo2"      "Skap2"     "Dmap1"     "Dnase1l1" 
# [113] "Reep5"     "Smim3" 



#What if we look at the correlation between datasets again, 
#but only using the top genes found in all 7 datasets?
sum(metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>6, na.rm=TRUE)
#[1] 361

TopGenesInAll11Datasets<-row.names(metaOutputFDR)[metaOutputFDR[,7]<0.05 & 
                                                    metaOutputFDR[,6]>6 & 
                                                    is.na(metaOutputFDR[,6])==FALSE]

cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges$x%in%TopGenesInAll11Datasets),-1]), use="pairwise.complete.obs")

# GSE109315_StressResilient_Vs_Control
# GSE109315_StressResilient_Vs_Control                             1.00000000
# GSE109315_StressSusceptible_Vs_Control                          -0.01256176
# GSE114224_Left_CSDS_vs_Control                                   0.19657085
# GSE114224_Right_CSDS_vs_Control                                  0.20559859
# GSE81672_StressResilient_Vs_Control                              0.29842824
# GSE81672_StressSusceptible_Vs_Control                            0.44087285
# GSE85472_CSDS_vs_Control                                         0.15369108
# GSE109315_StressSusceptible_Vs_Control GSE114224_Left_CSDS_vs_Control
# GSE109315_StressResilient_Vs_Control                              -0.01256176                      0.1965709
# GSE109315_StressSusceptible_Vs_Control                             1.00000000                      0.1284864
# GSE114224_Left_CSDS_vs_Control                                     0.12848638                      1.0000000
# GSE114224_Right_CSDS_vs_Control                                    0.11788957                      0.9208452
# GSE81672_StressResilient_Vs_Control                                0.06646361                      0.4152583
# GSE81672_StressSusceptible_Vs_Control                             -0.14611229                      0.2142673
# GSE85472_CSDS_vs_Control                                           0.02023647                      0.2780679
# GSE114224_Right_CSDS_vs_Control GSE81672_StressResilient_Vs_Control
# GSE109315_StressResilient_Vs_Control                         0.2055986                          0.29842824
# GSE109315_StressSusceptible_Vs_Control                       0.1178896                          0.06646361
# GSE114224_Left_CSDS_vs_Control                               0.9208452                          0.41525834
# GSE114224_Right_CSDS_vs_Control                              1.0000000                          0.41267188
# GSE81672_StressResilient_Vs_Control                          0.4126719                          1.00000000
# GSE81672_StressSusceptible_Vs_Control                        0.2366054                          0.33841478
# GSE85472_CSDS_vs_Control                                     0.3204683                          0.19221045
# GSE81672_StressSusceptible_Vs_Control GSE85472_CSDS_vs_Control
# GSE109315_StressResilient_Vs_Control                               0.4408728               0.15369108
# GSE109315_StressSusceptible_Vs_Control                            -0.1461123               0.02023647
# GSE114224_Left_CSDS_vs_Control                                     0.2142673               0.27806793
# GSE114224_Right_CSDS_vs_Control                                    0.2366054               0.32046827
# GSE81672_StressResilient_Vs_Control                                0.3384148               0.19221045
# GSE81672_StressSusceptible_Vs_Control                              1.0000000               0.16844829
# GSE85472_CSDS_vs_Control                                           0.1684483               1.00000000

pdf("Heatmap_CorMatrix_HPC_StressDatasets_GenesFDR05inAll.pdf", height=10, width=10)
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges$x%in%TopGenesInAll11Datasets),
                                               -1]), use="pairwise.complete.obs"))
dev.off()

CorMatrixTopGenes<-cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges
                                                                $x%in%TopGenesInAll11Datasets),-1]),
                       use="pairwise.complete.obs")

write.csv(CorMatrixTopGenes, "CorMatrix_PFC_StressDatasets_GenesFDR05inAll.csv")

#read the data for HC
Stress_HC_MetaAnalysis<-read.csv("metaOutputFDR_orderedByPval-HC.csv", header=TRUE, stringsAsFactors=FALSE)

str(Stress_HC_MetaAnalysis)
# 'data.frame':	17889 obs. of  8 variables:
# $ X                    : chr  "Slc7a12" "Rpl35" "Zfp277" "Klhdc7a" ...
# $ Log2FC_estimate      : num  0.0767 -0.0259 -0.078 0.1455 0.1164 ...
# $ SE                   : num  0.0054 0.00187 0.00582 0.01116 0.00907 ...
# $ pval                 : num  6.46e-46 7.31e-44 5.41e-41 6.84e-39 9.93e-38 ...
# $ CI_lb                : num  0.0662 -0.0296 -0.0894 0.1237 0.0986 ...
# $ CI_ub                : num  0.0873 -0.0223 -0.0666 0.1674 0.1342 ...
# $ Number_Of_Comparisons: int  7 7 10 10 11 11 11 7 9 11 ...
# $ FDR                  : num  1.16e-41 6.54e-40 3.23e-37 3.06e-35 3.55e-34 ...

metaOutputFDR_AsDF<-data.frame(X=row.names(metaOutputFDR), metaOutputFDR)

str(metaOutputFDR_AsDF)
# 'data.frame':	25334 obs. of  8 variables:
# $ X                    : chr  "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# $ Log2FC_estimate      : num  -0.07545 -0.00625 -0.23821 0.11246 0.0194 ...
# $ SE                   : num  0.23953 0.00829 0.33025 0.09458 0.00656 ...
# $ pval                 : num  0.75275 0.45111 0.47072 0.23441 0.00309 ...
# $ CI_lb                : num  -0.54492 -0.0225 -0.88549 -0.07291 0.00655 ...
# $ CI_ub                : num  0.394 0.01 0.4091 0.2978 0.0323 ...
# $ Number_Of_Comparisons: num  4 7 4 5 7 7 5 5 5 5 ...
# $ FDR                  : num  0.952 0.869 0.878 0.748 0.107 ...

# Comparing acute vs chronic in FC and PFC
Stress_PFC_vs_FC<-join(metaOutputFDR_AsDF,Stress_HC_MetaAnalysis, by="X", type="inner")

str(Stress_PFC_vs_HC)
# 'data.frame':	17657 obs. of  15 variables:
# $ X                    : chr  "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010K14Rik" ...
# $ Log2FC_estimate      : num  -0.07545 -0.00625 0.11246 0.0194 0.00716 ...
# $ SE                   : num  0.23953 0.00829 0.09458 0.00656 0.07153 ...
# $ pval                 : num  0.75275 0.45111 0.23441 0.00309 0.92024 ...
# $ CI_lb                : num  -0.54492 -0.0225 -0.07291 0.00655 -0.13303 ...
# $ CI_ub                : num  0.394 0.01 0.2978 0.0323 0.1474 ...
# $ Number_Of_Comparisons: num  4 7 5 7 5 5 5 5 5 7 ...
# $ FDR                  : num  0.952 0.869 0.748 0.107 0.985 ...
# $ Log2FC_estimate      : num  -0.03251 -0.000632 -0.088493 -0.065414 -0.009939 ...
# $ SE                   : num  0.11665 0.05085 0.03212 0.03927 0.00754 ...
# $ pval                 : num  0.78047 0.99008 0.00586 0.09578 0.18742 ...
# $ CI_lb                : num  -0.2611 -0.1003 -0.1514 -0.1424 -0.0247 ...
# $ CI_ub                : num  0.19612 0.09904 -0.02555 0.01156 0.00484 ...
# $ Number_Of_Comparisons: int  7 10 10 10 9 10 7 9 10 10 ...
# $ FDR                  : num  0.916 1 0.132 0.479 0.554 ...

plot(Stress_PFC_vs_HC[,2]~Stress_PFC_vs_HC[,9])

Stress_PFC_vs_HC[(abs(Stress_PFC_vs_HC[,2])>0.2)&(abs(Stress_PFC_vs_HC[,9])>0.2),]
# X Log2FC_estimate        SE       pval       CI_lb       CI_ub Number_Of_Comparisons
# 119   1700082M22Rik       0.2235603 0.1455888 0.12464641 -0.06178853  0.50890922                     4
# 433   4933432I03Rik      -0.2606715 0.2293358 0.25569022 -0.71016134  0.18881843                     4
# 484   6530403H02Rik      -0.2910372 0.3305137 0.37855592 -0.93883221  0.35675777                     4
# 522   9330199G10Rik       0.6311152 0.5914906 0.28597585 -0.52818506  1.79041548                     4
# 612   A730090N16Rik       0.2541312 0.2636561 0.33510908 -0.26262515  0.77088761                     4
# 1570       Arhgef16       0.2221442 0.2405331 0.35572095 -0.24929200  0.69358042                     5
# 2386  C430049B03Rik      -0.3277429 0.1494414 0.02829869 -0.62064261 -0.03484325                     4
# 2881         Cd209f      -0.6668075 0.5198319 0.19958436 -1.68565931  0.35204440                     5
# NA             <NA>              NA        NA         NA          NA          NA                    NA
# 6389        Gm11732       0.3972515 0.2396428 0.09738171 -0.07243975  0.86694271                     4
# 6395        Gm12108       0.2348505 0.2202047 0.28619327 -0.19674281  0.66644372                     4
# 6461        Gm20337      -0.3965294 0.3217454 0.21778767 -1.02713870  0.23407998                     4
# 9617         Ms4a14       0.5923883 0.3417952 0.08306590 -0.07751794  1.26229460                     5
# 15574      Tmem184a       0.5516624 0.3732330 0.13939057 -0.17986088  1.28318572                     7
# NA.1           <NA>              NA        NA         NA          NA          NA                    NA
# FDR Log2FC_estimate         SE         pval        CI_lb        CI_ub Number_Of_Comparisons
# 119   0.6287917       0.2062305 0.07022093 3.315319e-03  0.068599979  0.343860960                     8
# 433   0.7652491      -0.4192244 0.20881160 4.467855e-02 -0.828487663 -0.009961221                     8
# 484   0.8374974       0.2443910 0.15060793 1.046542e-01 -0.050795071  0.539577151                     7
# 522   0.7850582      -0.2511770 0.17950036 1.617194e-01 -0.602991287  0.100637208                     7
# 612   0.8168879      -0.2542777 0.16960089 1.338034e-01 -0.586689312  0.078133976                     7
# 1570  0.8256493       0.2434992 0.12166607 4.535215e-02  0.005038136  0.481960363                    10
# 2386  0.3496850      -0.2534091 0.32214732 4.315014e-01 -0.884806283  0.377988021                     7
# 2881  0.7150714       0.3057158 0.04618293 3.599974e-11  0.215198939  0.396232682                     9
# NA           NA              NA         NA           NA           NA           NA                    NA
# 6389  0.5835465      -0.2140771 0.29650854 4.702989e-01 -0.795223112  0.367069002                     7
# 6395  0.7850923       0.2598875 0.15708417 9.803660e-02 -0.047991791  0.567766854                     8
# 6461  0.7323189      -0.5955608 0.35320459 9.176424e-02 -1.287829087  0.096707461                     7
# 9617  0.5487030      -0.3227939 0.39877787 4.182517e-01 -1.104384194  0.458796314                     7
# 15574 0.6497606       0.2039885 0.08824029 2.079206e-02  0.031040739  0.376936301                     9
# NA.1         NA              NA         NA           NA           NA           NA                    NA
# FDR
# 119   8.958707e-02
# 433   3.826315e-01
# 484   4.856387e-01
# 522   5.334939e-01
# 612   5.103541e-01
# 1570  3.845046e-01
# 2386  7.193386e-01
# 2881  1.399999e-08
# NA              NA
# 6389  7.413143e-01
# 6395  4.820717e-01
# 6461  4.740483e-01
# 9617  7.105643e-01
# 15574 2.718927e-01
# NA.1            NA

Stress_PFC_vs_HC[(Stress_PFC_vs_HC[,4]<0.05)&(Stress_PFC_vs_HC[,11]<0.05),]

P05inBothPFCandHC<-Stress_PFC_vs_HC[(Stress_PFC_vs_HC[,4]<0.05)&(Stress_PFC_vs_HC[,11]<0.05),]
str(P05inBothPFCandHC)

plot(P05inBothPFCandHC[,2]~P05inBothPFCandHC[,9])
P05inBothPFCandHC[(abs(P05inBothPFCandHC[,2])>0.1)&(abs(P05inBothPFCandHC[,9])>0.1),]

#read the data for HC
Stress_FC_MetaAnalysis<-read.csv("metaOutputFDR_orderedByPval_wHDRFData-FC.csv", header=TRUE, stringsAsFactors=FALSE)

str(Stress_FC_MetaAnalysis)
# 'data.frame':	18523 obs. of  8 variables:
#   $ X                    : chr  "Slc12a6" "Thg1l" "B3gnt3" "Emc8" ...
# $ Log2FC_estimate      : num  -0.0812 0.103 0.1548 0.0837 -0.0805 ...
# $ SE                   : num  0.0114 0.0147 0.0231 0.0142 0.014 ...
# $ pval                 : num  1.33e-12 2.41e-12 1.98e-11 3.76e-09 7.91e-09 ...
# $ CI_lb                : num  -0.1036 0.0742 0.1095 0.0559 -0.1079 ...
# $ CI_ub                : num  -0.0587 0.1318 0.2 0.1116 -0.0532 ...
# $ Number_Of_Comparisons: int  16 18 18 11 18 18 18 18 18 18 ...
# $ FDR                  : num  2.23e-08 2.23e-08 1.22e-07 1.74e-05 2.93e-05 6.40e-05 7.69e-05 7.69e-05 8.44e-05 8.63e-05 ...

metaOutputFDR_AsDF<-data.frame(X=row.names(metaOutputFDR), metaOutputFDR)

str(metaOutputFDR_AsDF)
# 'data.frame':	25334 obs. of  8 variables:
# $ X                    : chr  "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# $ Log2FC_estimate      : num  -0.07545 -0.00625 -0.23821 0.11246 0.0194 ...
# $ SE                   : num  0.23953 0.00829 0.33025 0.09458 0.00656 ...
# $ pval                 : num  0.75275 0.45111 0.47072 0.23441 0.00309 ...
# $ CI_lb                : num  -0.54492 -0.0225 -0.88549 -0.07291 0.00655 ...
# $ CI_ub                : num  0.394 0.01 0.4091 0.2978 0.0323 ...
# $ Number_Of_Comparisons: num  4 7 4 5 7 7 5 5 5 5 ...
# $ FDR                  : num  0.952 0.869 0.878 0.748 0.107 ...

Stress_PFC_vs_FC<-join(metaOutputFDR_AsDF,Stress_FC_MetaAnalysis, by="X", type="inner")

str(Stress_PFC_vs_FC)
# 'data.frame':	18427 obs. of  15 variables:
#   $ X                    : chr  "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# $ Log2FC_estimate      : num  -0.07545 -0.00625 -0.23821 0.11246 0.0194 ...
# $ SE                   : num  0.23953 0.00829 0.33025 0.09458 0.00656 ...
# $ pval                 : num  0.75275 0.45111 0.47072 0.23441 0.00309 ...
# $ CI_lb                : num  -0.54492 -0.0225 -0.88549 -0.07291 0.00655 ...
# $ CI_ub                : num  0.394 0.01 0.4091 0.2978 0.0323 ...
# $ Number_Of_Comparisons: num  4 7 4 5 7 7 5 5 5 5 ...
# $ FDR                  : num  0.952 0.869 0.878 0.748 0.107 ...
# $ Log2FC_estimate      : num  0.027 -0.0382 -0.0734 0.0241 0.017 ...
# $ SE                   : num  0.0836 0.0336 0.0407 0.0573 0.0386 ...
# $ pval                 : num  0.7466 0.2558 0.0714 0.6746 0.6588 ...
# $ CI_lb                : num  -0.1369 -0.1041 -0.1533 -0.0883 -0.0585 ...
# $ CI_ub                : num  0.19092 0.02769 0.00641 0.13638 0.0926 ...
# $ Number_Of_Comparisons: int  16 18 15 18 18 18 14 18 16 14 ...
# $ FDR                  : num  0.948 0.739 0.5 0.925 0.92 ...

plot(Stress_PFC_vs_FC[,2]~Stress_PFC_vs_FC[,9])

Stress_PFC_vs_FC[(abs(Stress_PFC_vs_FC[,2])>0.2)&(abs(Stress_PFC_vs_FC[,9])>0.2),]
# X Log2FC_estimate        SE       pval       CI_lb       CI_ub Number_Of_Comparisons
# 119   1700082M22Rik       0.2235603 0.1455888 0.12464641 -0.06178853  0.50890922                     4
# 433   4933432I03Rik      -0.2606715 0.2293358 0.25569022 -0.71016134  0.18881843                     4
# 484   6530403H02Rik      -0.2910372 0.3305137 0.37855592 -0.93883221  0.35675777                     4
# 522   9330199G10Rik       0.6311152 0.5914906 0.28597585 -0.52818506  1.79041548                     4
# 612   A730090N16Rik       0.2541312 0.2636561 0.33510908 -0.26262515  0.77088761                     4
# 1570       Arhgef16       0.2221442 0.2405331 0.35572095 -0.24929200  0.69358042                     5
# 2386  C430049B03Rik      -0.3277429 0.1494414 0.02829869 -0.62064261 -0.03484325                     4
# 2881         Cd209f      -0.6668075 0.5198319 0.19958436 -1.68565931  0.35204440                     5
# NA             <NA>              NA        NA         NA          NA          NA                    NA
# 6389        Gm11732       0.3972515 0.2396428 0.09738171 -0.07243975  0.86694271                     4
# 6395        Gm12108       0.2348505 0.2202047 0.28619327 -0.19674281  0.66644372                     4
# 6461        Gm20337      -0.3965294 0.3217454 0.21778767 -1.02713870  0.23407998                     4
# 9617         Ms4a14       0.5923883 0.3417952 0.08306590 -0.07751794  1.26229460                     5
# 15574      Tmem184a       0.5516624 0.3732330 0.13939057 -0.17986088  1.28318572                     7
# NA.1           <NA>              NA        NA         NA          NA          NA                    NA
# FDR Log2FC_estimate         SE         pval        CI_lb        CI_ub Number_Of_Comparisons
# 119   0.6287917       0.2062305 0.07022093 3.315319e-03  0.068599979  0.343860960                     8
# 433   0.7652491      -0.4192244 0.20881160 4.467855e-02 -0.828487663 -0.009961221                     8
# 484   0.8374974       0.2443910 0.15060793 1.046542e-01 -0.050795071  0.539577151                     7
# 522   0.7850582      -0.2511770 0.17950036 1.617194e-01 -0.602991287  0.100637208                     7
# 612   0.8168879      -0.2542777 0.16960089 1.338034e-01 -0.586689312  0.078133976                     7
# 1570  0.8256493       0.2434992 0.12166607 4.535215e-02  0.005038136  0.481960363                    10
# 2386  0.3496850      -0.2534091 0.32214732 4.315014e-01 -0.884806283  0.377988021                     7
# 2881  0.7150714       0.3057158 0.04618293 3.599974e-11  0.215198939  0.396232682                     9
# NA           NA              NA         NA           NA           NA           NA                    NA
# 6389  0.5835465      -0.2140771 0.29650854 4.702989e-01 -0.795223112  0.367069002                     7
# 6395  0.7850923       0.2598875 0.15708417 9.803660e-02 -0.047991791  0.567766854                     8
# 6461  0.7323189      -0.5955608 0.35320459 9.176424e-02 -1.287829087  0.096707461                     7
# 9617  0.5487030      -0.3227939 0.39877787 4.182517e-01 -1.104384194  0.458796314                     7
# 15574 0.6497606       0.2039885 0.08824029 2.079206e-02  0.031040739  0.376936301                     9
# NA.1         NA              NA         NA           NA           NA           NA                    NA
# FDR
# 119   8.958707e-02
# 433   3.826315e-01
# 484   4.856387e-01
# 522   5.334939e-01
# 612   5.103541e-01
# 1570  3.845046e-01
# 2386  7.193386e-01
# 2881  1.399999e-08
# NA              NA
# 6389  7.413143e-01
# 6395  4.820717e-01
# 6461  4.740483e-01
# 9617  7.105643e-01
# 15574 2.718927e-01
# NA.1            NA

Stress_PFC_vs_FC[(Stress_PFC_vs_FC[,4]<0.05)&(Stress_PFC_vs_FC[,11]<0.05),]

P05inBothPFCandFC<-Stress_PFC_vs_FC[(Stress_PFC_vs_FC[,4]<0.05)&(Stress_PFC_vs_FC[,11]<0.05),]
str(P05inBothPFCandFC)
# 'data.frame':	255 obs. of  15 variables:
#   $ X                    : chr  "1700020G03Rik" "1700025G04Rik" "2700062C07Rik" "3110082I17Rik" ...
# $ Log2FC_estimate      : num  0.6874 -0.0142 0.02 0.0172 0.4888 ...
# $ SE                   : num  0.33057 0.00585 0.0077 0.00655 0.23741 ...
# $ pval                 : num  0.03758 0.01553 0.00957 0.00865 0.03949 ...
# $ CI_lb                : num  0.03948 -0.02564 0.00486 0.00436 0.02353 ...
# $ CI_ub                : num  1.3353 -0.00269 0.03505 0.03002 0.95416 ...
# $ Number_Of_Comparisons: num  4 7 7 7 5 4 4 4 7 7 ...
# $ FDR                  : num  0.396 0.262 0.204 0.194 0.406 ...
# $ Log2FC_estimate      : num  0.138 0.0338 0.0731 0.1045 0.2522 ...
# $ SE                   : num  0.0691 0.0155 0.0308 0.0329 0.1194 ...
# $ pval                 : num  0.04592 0.02897 0.01758 0.00146 0.03476 ...
# $ CI_lb                : num  0.0025 0.00347 0.01276 0.04014 0.01806 ...
# $ CI_ub                : num  0.2734 0.0642 0.1335 0.1689 0.4862 ...
# $ Number_Of_Comparisons: int  11 18 16 18 18 13 16 11 18 18 ...
# $ FDR                  : num  0.4309 0.3625 0.2948 0.0845 0.3865 ...

plot(P05inBothPFCandFC[,2]~P05inBothPFCandFC[,9])
P05inBothPFCandHC[(abs(P05inBothPFCandFC[,2])>0.1)&(abs(P05inBothPFCandFC[,9])>0.1),]
