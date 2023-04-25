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
setwd("~/Documents/2022 Fall/Research/CSDS and CUMS/Gemma Data")

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

ListOfDEResults<-list(DEResults_GSE109315, DEResults_GSE114224_Left, DEResults_GSE114224_Right, DEResults_GSE81672, DEResults_GSE85472,DEResults_GSE81587,DEResults_GSE102556.1,DEResults_GSE151807)



AligningDatasets(ListOfDEResults)
# [1] "MetaAnalysis_FoldChange_Dfs:"
# List of 8
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
# $ :'data.frame':	23318 obs. of  2 variables:
#   ..$ x                       : chr [1:23318] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ GSE81587_CUMS_vs_Control: num [1:23318] 1.536 -0.0853 -0.1828 -0.1766 0.0161 ...
# $ :'data.frame':	24789 obs. of  2 variables:
#   ..$ x                          : chr [1:24789] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ GSE102556.1_CUMS_vs_Control: num [1:24789] 0.3705 -0.00622 0.2275 -0.2992 -0.06806 ...
# $ :'data.frame':	18537 obs. of  2 variables:
#   ..$ x                        : chr [1:18537] "0610005C13Rik" "0610009B22Rik" "0610009L18Rik" "0610010K14Rik" ...
# ..$ GSE151807_CUMS_vs_Control: num [1:18537] -0.033 -0.1549 0.1254 0.036 -0.0529 ...
# NULL
# [1] "MetaAnalysis_FoldChanges:"
# 'data.frame':	35237 obs. of  11 variables:
#   $ x                                     : chr  "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# $ GSE109315_StressResilient_Vs_Control  : num  -0.084 -0.0858 0.0425 -0.408 0.2268 ...
# $ GSE109315_StressSusceptible_Vs_Control: num  0.1929 -0.0733 -0.0991 -0.0613 0.3685 ...
# $ GSE114224_Left_CSDS_vs_Control        : num  NA NA -0.0093 NA NA ...
# $ GSE114224_Right_CSDS_vs_Control       : num  NA NA 0.00277 NA NA ...
# $ GSE81672_StressResilient_Vs_Control   : num  -0.7165 NA -0.0636 0.6613 0.1336 ...
# $ GSE81672_StressSusceptible_Vs_Control : num  -1.667 NA -0.0998 -1.264 -0.2959 ...
# $ GSE85472_CSDS_vs_Control              : num  NA NA -0.035 NA 0.0365 ...
# $ GSE81587_CUMS_vs_Control              : num  1.536 NA -0.0853 -0.1828 -0.1766 ...
# $ GSE102556.1_CUMS_vs_Control           : num  0.3705 NA -0.00622 0.2275 -0.2992 ...
# $ GSE151807_CUMS_vs_Control             : num  -0.033 NA -0.155 NA 0.125 ...
# NULL
# [1] "MetaAnalysis_SV_Dfs:"
# List of 8
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
# 'data.frame':	35237 obs. of  11 variables:
#   $ x                                     : chr  "0610005C13Rik" "0610006L08Rik" "0610009B22Rik" "0610009E02Rik" ...
# $ GSE109315_StressResilient_Vs_Control  : num  0.1704 0.0202 0.0063 0.1143 0.0344 ...
# $ GSE109315_StressSusceptible_Vs_Control: num  0.09173 0.01089 0.00339 0.06156 0.01856 ...
# $ GSE114224_Left_CSDS_vs_Control        : num  NA NA 0.000191 NA NA ...
# $ GSE114224_Right_CSDS_vs_Control       : num  NA NA 0.000124 NA NA ...
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
# cor(as.matrix(MetaAnalysis_FoldChanges[,-1]), use="pairwise.complete.obs")
# GSE109315_StressResilient_Vs_Control GSE109315_StressSusceptible_Vs_Control
# GSE109315_StressResilient_Vs_Control                            1.000000000                            0.209603799
# GSE109315_StressSusceptible_Vs_Control                          0.209603799                            1.000000000
# GSE114224_Left_CSDS_vs_Control                                 -0.009182396                            0.089912682
# GSE114224_Right_CSDS_vs_Control                                 0.062960251                           -0.024486299
# GSE81672_StressResilient_Vs_Control                             0.041576271                            0.058634459
# GSE81672_StressSusceptible_Vs_Control                           0.053830157                           -0.022506798
# GSE85472_CSDS_vs_Control                                       -0.020812101                            0.088833896
# GSE81587_CUMS_vs_Control                                        0.001838177                            0.041118824
# GSE102556.1_CUMS_vs_Control                                    -0.005626461                            0.010055882
# GSE151807_CUMS_vs_Control                                       0.004870402                           -0.001522263
# GSE114224_Left_CSDS_vs_Control GSE114224_Right_CSDS_vs_Control
# GSE109315_StressResilient_Vs_Control                     -0.009182396                     0.062960251
# GSE109315_StressSusceptible_Vs_Control                    0.089912682                    -0.024486299
# GSE114224_Left_CSDS_vs_Control                            1.000000000                     0.219358022
# GSE114224_Right_CSDS_vs_Control                           0.219358022                     1.000000000
# GSE81672_StressResilient_Vs_Control                       0.099781605                     0.050370225
# GSE81672_StressSusceptible_Vs_Control                    -0.011179858                     0.067188740
# GSE85472_CSDS_vs_Control                                  0.094358420                    -0.057335481
# GSE81587_CUMS_vs_Control                                 -0.019236298                     0.013345454
# GSE102556.1_CUMS_vs_Control                               0.039552231                    -0.027798898
# GSE151807_CUMS_vs_Control                                 0.007359793                    -0.007076517
# GSE81672_StressResilient_Vs_Control GSE81672_StressSusceptible_Vs_Control
# GSE109315_StressResilient_Vs_Control                            0.04157627                           0.053830157
# GSE109315_StressSusceptible_Vs_Control                          0.05863446                          -0.022506798
# GSE114224_Left_CSDS_vs_Control                                  0.09978160                          -0.011179858
# GSE114224_Right_CSDS_vs_Control                                 0.05037023                           0.067188740
# GSE81672_StressResilient_Vs_Control                             1.00000000                           0.455978932
# GSE81672_StressSusceptible_Vs_Control                           0.45597893                           1.000000000
# GSE85472_CSDS_vs_Control                                        0.05527107                          -0.015546804
# GSE81587_CUMS_vs_Control                                       -0.01957326                          -0.002982487
# GSE102556.1_CUMS_vs_Control                                    -0.01466969                          -0.074005675
# GSE151807_CUMS_vs_Control                                       0.00992222                           0.052802083
# GSE85472_CSDS_vs_Control GSE81587_CUMS_vs_Control GSE102556.1_CUMS_vs_Control
# GSE109315_StressResilient_Vs_Control               -0.020812101              0.001838177                -0.005626461
# GSE109315_StressSusceptible_Vs_Control              0.088833896              0.041118824                 0.010055882
# GSE114224_Left_CSDS_vs_Control                      0.094358420             -0.019236298                 0.039552231
# GSE114224_Right_CSDS_vs_Control                    -0.057335481              0.013345454                -0.027798898
# GSE81672_StressResilient_Vs_Control                 0.055271070             -0.019573262                -0.014669693
# GSE81672_StressSusceptible_Vs_Control              -0.015546804             -0.002982487                -0.074005675
# GSE85472_CSDS_vs_Control                            1.000000000             -0.039247153                 0.005804135
# GSE81587_CUMS_vs_Control                           -0.039247153              1.000000000                 0.041525439
# GSE102556.1_CUMS_vs_Control                         0.005804135              0.041525439                 1.000000000
# GSE151807_CUMS_vs_Control                          -0.094564384              0.110097305                 0.034211134
# GSE151807_CUMS_vs_Control
# GSE109315_StressResilient_Vs_Control                 0.004870402
# GSE109315_StressSusceptible_Vs_Control              -0.001522263
# GSE114224_Left_CSDS_vs_Control                       0.007359793
# GSE114224_Right_CSDS_vs_Control                     -0.007076517
# GSE81672_StressResilient_Vs_Control                  0.009922220
# GSE81672_StressSusceptible_Vs_Control                0.052802083
# GSE85472_CSDS_vs_Control                            -0.094564384
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
# 0     1     2     3     4     5     6     7     8     9 
# 12163  1257  2349  1766  5184  1444  1392  2453  5366  1863 


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
NumberOfComparisons=11
CutOffForNAs=5
#I want at least 7 datasets (so 5 na's is too many)

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

# [1] "Table of # of NAs per Row (Gene):"
# MetaAnalysis_FoldChanges_NAsPerRow
# 0     1     2     3     4     5     6     7     8     9 
# 12163  1257  2349  1766  5184  1444  1392  2453  5366  1863 
# [1] "MetaAnalysis_FoldChanges_ForMeta:"
# 'data.frame':	22719 obs. of  11 variables:
#   $ x                                     : chr  "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# $ GSE109315_StressResilient_Vs_Control  : num  -0.084 0.0425 -0.408 0.2268 -0.0261 ...
# $ GSE109315_StressSusceptible_Vs_Control: num  0.1929 -0.0991 -0.0613 0.3685 0.0872 ...
# $ GSE114224_Left_CSDS_vs_Control        : num  NA -0.0093 NA NA 0.0218 ...
# $ GSE114224_Right_CSDS_vs_Control       : num  NA 0.00277 NA NA 0.02463 ...
# $ GSE81672_StressResilient_Vs_Control   : num  -0.7165 -0.0636 0.6613 0.1336 0.0292 ...
# $ GSE81672_StressSusceptible_Vs_Control : num  -1.667 -0.0998 -1.264 -0.2959 -0.0832 ...
# $ GSE85472_CSDS_vs_Control              : num  NA -0.035 NA 0.0365 -0.116 ...
# $ GSE81587_CUMS_vs_Control              : num  1.536 -0.0853 -0.1828 -0.1766 0.0161 ...
# $ GSE102556.1_CUMS_vs_Control           : num  0.3705 -0.00622 0.2275 -0.2992 -0.06806 ...
# $ GSE151807_CUMS_vs_Control             : num  -0.033 -0.155 NA 0.125 0.036 ...
# NULL
# There were 50 or more warnings (use warnings() to see the first 50)

str(metaOutput)
# num [1:22719, 1:6] 0.0167 -0.0403 -0.1559 0.061 0.0195 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:22719] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ : chr [1:6] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

head(metaOutput)
# Log2FC_estimate          SE        pval        CI_lb       CI_ub Number_Of_Comparisons
# 0610005C13Rik     0.016699109 0.122156587 0.891265966 -0.222723403 0.256121621                     8
# 0610009B22Rik    -0.040267582 0.021094106 0.056268619 -0.081611271 0.001076107                    11
# 0610009E02Rik    -0.155860725 0.203530042 0.443802860 -0.554772277 0.243050827                     7
# 0610009L18Rik     0.060979966 0.071008591 0.390468055 -0.078194314 0.200154246                     9
# 0610010K14Rik     0.019489738 0.006505182 0.002735172  0.006739815 0.032239661                    11
# 0610012G03Rik    -0.005078878 0.012221231 0.677718637 -0.029032051 0.018874295                    11

tail(metaOutput)
# Log2FC_estimate          SE       pval        CI_lb      CI_ub Number_Of_Comparisons
# Strbp    0.0009552082 0.006320727 0.87987856 -0.011433189 0.01334361                     9
# Syne1   -0.0005225234 0.005903905 0.92947554 -0.012093965 0.01104892                     9
# Ttc3     0.0148097236 0.015399964 0.33621408 -0.015373651 0.04499310                     8
# Usp9x    0.0017593789 0.004807233 0.71437571 -0.007662626 0.01118138                     9
# Cltc    -0.0134662132 0.016368860 0.41069417 -0.045548589 0.01861616                     7
# Macf1    0.0536989659 0.024213487 0.02657328  0.006241404 0.10115653                     7
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
# num [1:22719, 1:2] 0.89127 0.05627 0.4438 0.39047 0.00274 ...
# - attr(*, "dimnames")=List of 2
# ..$ : NULL
# ..$ : chr [1:2] "rawp" "BH"

#adjusted pvalue object is in same orientation as metaOutput so can simply be binded together
#We can double-check that it is in the same order by comparing the uncorrected p-values in the metaPvalAdj object to the pvalues in the metaOutput matrix

plot(metaPvalAdj[,1]~metaOutput[,3])
#Perfectly straight line - they are exactly the same! *phew*

metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])

str(metaOutputFDR)
# num [1:22719, 1:7] 0.0167 -0.0403 -0.1559 0.061 0.0195 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:22719] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
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
# num [1:19838, 1:7] 0.0167 -0.03858 -0.05958 0.0753 0.00285 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:19838] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following false discovery rate correction?"
# [1] 494
# [1] "What are the top results?"
# Log2FC_estimate          SE         pval       CI_lb       CI_ub Number_Of_Comparisons          FDR
# Col5a1       0.04374801 0.003085410 1.236047e-45  0.03770072  0.04979530                    11 2.452069e-41
# Nars        -0.03089716 0.003375091 5.462673e-20 -0.03751221 -0.02428210                    11 5.418426e-16
# Shisa7       0.03192035 0.003683342 4.470208e-18  0.02470114  0.03913957                    11 2.955999e-14
# Rgs10       -0.03190231 0.003830416 8.174991e-17 -0.03940978 -0.02439483                    11 4.054387e-13
# Tmem267      0.03030467 0.003933888 1.323993e-14  0.02259439  0.03801495                    11 5.253075e-11
# Prmt1        0.03296609 0.004337612 2.960050e-14  0.02446453  0.04146765                    11 9.786914e-11
############################

#10) Determine which are the top differentially expressed genes and create forest 
#plots to visualize the effect sizes for those top differentially expressed genes 
#across the different studies. 

row.names(metaOutputFDR_OrderbyPval)[c(1:100)]
# [1] "Col5a1"        "Nars"          "Rgs10"         "Shisa7"        "Tmem267"       "Prmt1"         "Tceal3"       
# [8] "Prr7"          "Stk32c"        "Usp32"         "Ramp3"         "Pou1f1"        "P2ry12"        "Tdrd7"        
# [15] "Zmat3"         "1810021B22Rik" "Rad23b"        "Actr1a"        "Cmpk2"         "Pcdhb20"       "Lor"          
# [22] "Klhl7"         "Acsbg1"        "Lrrtm3"        "Map10"         "Ptk2"          "Gmds"          "Med20"        
# [29] "Lonrf2"        "Pex26"         "Adam9"         "Rwdd4a"        "Sqor"          "Lrpap1"        "Gm26782"      
# [36] "Exosc7"        "Zfx"           "Opn3"          "Crym"          "Hax1"          "Spink10"       "Stard3nl"     
# [43] "Hsd11b1"       "Zfp397"        "Tpst1"         "Rora"          "Ccdc85b"       "Clpp"          "Nudt18"       
# [50] "Pycr2"         "Ankrd39"       "Amz2"          "Slc7a10"       "Slc39a12"      "Srcin1"        "Trappc6b"     
# [57] "Hnrnpa2b1"     "Casp8"         "Nxnl2"         "Tubg2"         "Hdgfl2"        "Kcnq2"         "Kmt5c"        
# [64] "Blvrb"         "Ndfip1"        "Rab5c"         "Ptdss1"        "Ptgds"         "Smg6"          "Ssx2ip"       
# [71] "Impdh2"        "Ndufa5"        "Htr1d"         "Dbpht2"        "Top2b"         "Snrnp27"       "Dnase1l1"     
# [78] "Eef2"          "Or13c7"        "Fam98a"        "Dmap1"         "Pdcd7"         "Lypd1"         "Fgf12"        
# [85] "Dis3"          "Hikeshi"       "Gm46920"       "Pld4"          "Aldoa"         "Cxcl13"        "Tmem128"      
# [92] "Cacybp"        "Gm46945"       "Rpp21"         "Smim3"         "Lpp"           "Nceh1"         "Ccdc102a"     
# [99] "Micos13"       "Mmp24"       

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
MakeForestPlots("Rgs10")
MakeForestPlots("Shisa7") 
MakeForestPlots("Tmem267") 
MakeForestPlots("Prmt1")
MakeForestPlots("Tceal3") # warning
MakeForestPlots("Prr7")# warning
MakeForestPlots("Stk32c")
# Warning message:
#   Studies with NAs omitted from model fitting. 



#In general, it seems to me like the results that have both a significant FDR 
#& large effect size tend to be the most convincing 

#Here's a summary of the distribution for the Log2FCs
summary(metaOutputFDR[,1])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -1.21067 -0.01378  0.00009 -0.00031  0.01191  1.58796       50 

# how many are statistically significant
sum(metaOutputFDR[,7]<0.05, na.rm = TRUE)
# [1] 334
# store significant data
metaOutputSig <- metaOutputFDR[metaOutputFDR[,7]<0.05,]
# remove na
metaOutputSigNoNa <- na.omit(metaOutputSig)
# count upregulated genes
metaOutputupregulated <- metaOutputSigNoNa[metaOutputSigNoNa[,1]>0,]
str(metaOutputupregulated)
# num [1:163, 1:7] 0.0218 0.0171 0.0128 0.0308 0.0203 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:163] "1700030K09Rik" "2310039H08Rik" "3110070M22Rik" "Actr1a" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
metaOutputdownregulated <- metaOutputSigNoNa[metaOutputSigNoNa[,1]<0,]
str(metaOutputdownregulated)
# num [1:171, 1:7] -0.104 -0.204 -0.272 -0.168 -0.377 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:171] "1700066M21Rik" "1810019D21Rik" "1810021B22Rik" "2010320M18Rik" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...



#Let's see how many results are both statistically significant 
#(using a more conservative FDR) and have a larger estimated Log2FC (>0.10 or <-0.10)
sum(metaOutputFDR[,7]<0.05 & abs(metaOutputFDR[,1])>0.1, na.rm=TRUE)
#[1] 44

#What are their gene symbols?
row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & 
                                       abs(metaOutputFDR_OrderbyPval[,1])>0.1 & 
                                       is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]
# [1] "1810021B22Rik" "Gm26782"       "Gm46920"       "Gm46945"       "Gm15631"       "2410018L13Rik" "Gm36163"      
# [8] "Gm30539"       "4930429B21Rik" "2010320M18Rik" "Tagln2"        "Gm38418"       "LOC108168680"  "Gm29994"      
# [15] "Gm30392"       "Gm46405"       "Gcnt4"         "Gm14199"       "Gm14124"       "Hmgb1"         "Cfap97d2"     
# [22] "Gm1043"        "Ckap2"         "Gm34140"       "Gm33931"       "1700066M21Rik" "Gm36603"       "Slc44a5"      
# [29] "Xlr3b"         "Gm34922"       "Gm31465"       "Gm30142"       "Cd99"          "Csmd2"         "Gm40564"      
# [36] "A330009N23Rik" "Mnd1"          "Gm11627"       "Gm33324"       "1810019D21Rik" "Gm12519"       "Gm14305"      
# [43] "Gm39805"       "Gm38561" 

MakeForestPlots("Gm26782")
# Warning message:
#   Studies with NAs omitted from model fitting. 
MakeForestPlots("Gm30539")
# Warning message:
#   Studies with NAs omitted from model fitting. 
MakeForestPlots("Gm15631") 
# Warning message:
#   Studies with NAs omitted from model fitting. 
MakeForestPlots("2410018L13Rik") 
# Warning message:
#   Studies with NAs omitted from model fitting. 

#They all look good

#Instead of effect sizes, that if I focus on genes represented in more than 5 datasets?

colnames(metaOutputFDR)
# [1] "Log2FC_estimate"       "SE"                    "pval"                  "CI_lb"                
# [5] "CI_ub"                 "Number_Of_Comparisons" "FDR" 

sum(metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>5, na.rm=TRUE)
#[1] 344

row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & 
                                       metaOutputFDR_OrderbyPval[,6]>5 & 
                                       is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]

#Adding in a p-value threshold too:
sum(metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>5 & metaOutputFDR[,3]<0.00001, na.rm=TRUE)
#[1] 84

row.names(metaOutputFDR_OrderbyPval)[metaOutputFDR_OrderbyPval[,7]<0.05 & 
                                       metaOutputFDR_OrderbyPval[,6]>5 & 
                                       metaOutputFDR_OrderbyPval[,3]<0.00001 & 
                                        is.na(metaOutputFDR_OrderbyPval[,1])==FALSE]
# [1] "Col5a1"        "Nars"          "Rgs10"         "Shisa7"        "Tmem267"       "Prmt1"         "Tceal3"       
# [8] "Prr7"          "Stk32c"        "Usp32"         "Ramp3"         "Pou1f1"        "P2ry12"        "Tdrd7"        
# [15] "Zmat3"         "1810021B22Rik" "Rad23b"        "Actr1a"        "Cmpk2"         "Pcdhb20"       "Lor"          
# [22] "Klhl7"         "Acsbg1"        "Lrrtm3"        "Map10"         "Ptk2"          "Gmds"          "Med20"        
# [29] "Lonrf2"        "Pex26"         "Adam9"         "Rwdd4a"        "Sqor"          "Lrpap1"        "Gm26782"      
# [36] "Exosc7"        "Zfx"           "Opn3"          "Crym"          "Hax1"          "Spink10"       "Stard3nl"     
# [43] "Hsd11b1"       "Zfp397"        "Tpst1"         "Rora"          "Ccdc85b"       "Clpp"          "Nudt18"       
# [50] "Pycr2"         "Ankrd39"       "Amz2"          "Slc7a10"       "Slc39a12"      "Srcin1"        "Trappc6b"     
# [57] "Hnrnpa2b1"     "Casp8"         "Nxnl2"         "Tubg2"         "Hdgfl2"        "Kcnq2"         "Kmt5c"        
# [64] "Blvrb"         "Ndfip1"        "Rab5c"         "Ptdss1"        "Ptgds"         "Smg6"          "Ssx2ip"       
# [71] "Impdh2"        "Ndufa5"        "Htr1d"         "Dbpht2"        "Top2b"         "Snrnp27"       "Dnase1l1"     
# [78] "Eef2"          "Or13c7"        "Fam98a"        "Dmap1"         "Pdcd7"         "Lypd1"         "Fgf12" 



#What if we look at the correlation between datasets again, 
#but only using the top genes found in all 7 datasets?
sum(metaOutputFDR[,7]<0.05 & metaOutputFDR[,6]>6, na.rm=TRUE)
#[1] 334

TopGenesInAll11Datasets<-row.names(metaOutputFDR)[metaOutputFDR[,7]<0.05 & 
                                                    metaOutputFDR[,6]>6 & 
                                                    is.na(metaOutputFDR[,6])==FALSE]

cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges$x%in%TopGenesInAll11Datasets),-1]), use="pairwise.complete.obs")
# GSE109315_StressResilient_Vs_Control GSE109315_StressSusceptible_Vs_Control
# GSE109315_StressResilient_Vs_Control                             1.00000000                              0.5375206
# GSE109315_StressSusceptible_Vs_Control                           0.53752056                              1.0000000
# GSE114224_Left_CSDS_vs_Control                                   0.13681265                              0.2019750
# GSE114224_Right_CSDS_vs_Control                                  0.12049250                              0.1967166
# GSE81672_StressResilient_Vs_Control                              0.46818413                              0.5896635
# GSE81672_StressSusceptible_Vs_Control                            0.52699504                              0.4097039
# GSE85472_CSDS_vs_Control                                         0.10861635                              0.2040945
# GSE81587_CUMS_vs_Control                                         0.19611815                              0.4029892
# GSE102556.1_CUMS_vs_Control                                      0.36607396                              0.4685374
# GSE151807_CUMS_vs_Control                                        0.07924778                              0.1567104
# GSE114224_Left_CSDS_vs_Control GSE114224_Right_CSDS_vs_Control
# GSE109315_StressResilient_Vs_Control                        0.1368127                      0.12049250
# GSE109315_StressSusceptible_Vs_Control                      0.2019750                      0.19671655
# GSE114224_Left_CSDS_vs_Control                              1.0000000                      0.91708130
# GSE114224_Right_CSDS_vs_Control                             0.9170813                      1.00000000
# GSE81672_StressResilient_Vs_Control                         0.3783022                      0.36447062
# GSE81672_StressSusceptible_Vs_Control                       0.1469268                      0.17146001
# GSE85472_CSDS_vs_Control                                    0.3053087                      0.35246259
# GSE81587_CUMS_vs_Control                                    0.1129144                      0.10399474
# GSE102556.1_CUMS_vs_Control                                 0.1139123                      0.08344762
# GSE151807_CUMS_vs_Control                                   0.2714816                      0.24452804
# GSE81672_StressResilient_Vs_Control GSE81672_StressSusceptible_Vs_Control
# GSE109315_StressResilient_Vs_Control                             0.4681841                             0.5269950
# GSE109315_StressSusceptible_Vs_Control                           0.5896635                             0.4097039
# GSE114224_Left_CSDS_vs_Control                                   0.3783022                             0.1469268
# GSE114224_Right_CSDS_vs_Control                                  0.3644706                             0.1714600
# GSE81672_StressResilient_Vs_Control                              1.0000000                             0.4809467
# GSE81672_StressSusceptible_Vs_Control                            0.4809467                             1.0000000
# GSE85472_CSDS_vs_Control                                         0.1883060                             0.1321393
# GSE81587_CUMS_vs_Control                                         0.3499361                             0.1458025
# GSE102556.1_CUMS_vs_Control                                      0.4277913                             0.3582670
# GSE151807_CUMS_vs_Control                                        0.2396167                             0.1146241
# GSE85472_CSDS_vs_Control GSE81587_CUMS_vs_Control GSE102556.1_CUMS_vs_Control
# GSE109315_StressResilient_Vs_Control                  0.1086164                0.1961182                  0.36607396
# GSE109315_StressSusceptible_Vs_Control                0.2040945                0.4029892                  0.46853741
# GSE114224_Left_CSDS_vs_Control                        0.3053087                0.1129144                  0.11391228
# GSE114224_Right_CSDS_vs_Control                       0.3524626                0.1039947                  0.08344762
# GSE81672_StressResilient_Vs_Control                   0.1883060                0.3499361                  0.42779135
# GSE81672_StressSusceptible_Vs_Control                 0.1321393                0.1458025                  0.35826704
# GSE85472_CSDS_vs_Control                              1.0000000                0.2111827                  0.07215520
# GSE81587_CUMS_vs_Control                              0.2111827                1.0000000                  0.19596813
# GSE102556.1_CUMS_vs_Control                           0.0721552                0.1959681                  1.00000000
# GSE151807_CUMS_vs_Control                             0.1887835                0.2149992                 -0.00974468
# GSE151807_CUMS_vs_Control
# GSE109315_StressResilient_Vs_Control                  0.07924778
# GSE109315_StressSusceptible_Vs_Control                0.15671041
# GSE114224_Left_CSDS_vs_Control                        0.27148158
# GSE114224_Right_CSDS_vs_Control                       0.24452804
# GSE81672_StressResilient_Vs_Control                   0.23961665
# GSE81672_StressSusceptible_Vs_Control                 0.11462414
# GSE85472_CSDS_vs_Control                              0.18878351
# GSE81587_CUMS_vs_Control                              0.21499920
# GSE102556.1_CUMS_vs_Control                          -0.00974468
# GSE151807_CUMS_vs_Control                             1.00000000

pdf("Heatmap_CorMatrix_PFC_StressDatasets_GenesFDR05inAll.pdf", height=10, width=10)
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges$x%in%TopGenesInAll11Datasets),
                                               -1]), use="pairwise.complete.obs"))
dev.off()

CorMatrixTopGenes<-cor(as.matrix(MetaAnalysis_FoldChanges[which(MetaAnalysis_FoldChanges
                                                                $x%in%TopGenesInAll11Datasets),-1]),
                       use="pairwise.complete.obs")

write.csv(CorMatrixTopGenes, "CorMatrix_PFC_StressDatasets_GenesFDR05inAll.csv")
