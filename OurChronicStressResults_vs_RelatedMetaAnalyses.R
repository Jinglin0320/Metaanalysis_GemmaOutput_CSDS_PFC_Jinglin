#Comparing Jinglin's Chronic Stress effects on the PFC to other related meta-analyses or large sample size experiments
#Megan Hagenauer
#2025-08-29

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2022_JinglinXiong_ChronicStress_FrontalCortex/Final Product_2025_Corrected/Reference Articles")

list.files()

DEGs_inOtherMeta<-read.csv("DEGs_vs_OtherMetaAnalyses.csv", header=TRUE, stringsAsFactors = FALSE)

colnames(DEGs_inOtherMeta)
str(DEGs_inOtherMeta)

#Some of the Stankiewicz_2022 columns may need a bit of cleaning - there seems to be blank spaces in them, e.g.,

table(DEGs_inOtherMeta$ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction)
# Down Down, Up       Up Up, Down 
# 66       26        6       29        6 

#in general, the entire table may suffer from that problem, since it was constructed in Excel
#Let's clean it up... or maybe it doesn't really matter for our current purposes? I double-checked the columns and it looks like they are identified as the correct data types.


###############

#Comparisons: Scatterplots & Correlations

transparent_grey <- rgb(red = 128, green = 128, blue = 128, alpha = 128, maxColorValue = 255)

###############

#Comparison with Toni's ELS Meta-Analysis Results:

pdf("Scatterplot_ELS_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(ELS_Duan_2025_Log2FC_estimate~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta, xlab="Chronic Stress: Log2FC", ylab="ELS: Log2FC")

points(ELS_Duan_2025_Log2FC_estimate~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta[DEGs_inOtherMeta$ELS_Duan_2025_pval<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(ELS_Duan_2025_Log2FC_estimate~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

#linear regression output for that relationship:
summary.lm(TrendLine)
rm(TrendLine)

# Call:
#   lm(formula = ELS_Duan_2025_Log2FC_estimate ~ ChronicStress_Xiong_Log2FC_estimate, 
#      data = DEGs_inOtherMeta)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.269439 -0.047357  0.002105  0.049250  0.264736 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          0.00797    0.01005   0.793     0.43    
# ChronicStress_Xiong_Log2FC_estimate  0.30907    0.05413   5.710 1.44e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07797 on 90 degrees of freedom
# (41 observations deleted due to missingness)
# Multiple R-squared:  0.2659,	Adjusted R-squared:  0.2578 
# F-statistic:  32.6 on 1 and 90 DF,  p-value: 1.437e-07

#Calculating the Pearson rank correlation:
cor.test(DEGs_inOtherMeta$ELS_Duan_2025_Log2FC_estimate, DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate, method="pearson")

# Pearson's product-moment correlation
# 
# data:  DEGs_inOtherMeta$ELS_Duan_2025_Log2FC_estimate and DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate
# t = 5.7099, df = 90, p-value = 1.437e-07
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3475687 0.6516650
# sample estimates:
#       cor 
# 0.5156761

#Calculating the Spearman rank correlation:
cor.test(DEGs_inOtherMeta$ELS_Duan_2025_Log2FC_estimate, DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate, method="spearman")

# Spearman's rank correlation rho
# 
# data:  DEGs_inOtherMeta$ELS_Duan_2025_Log2FC_estimate and DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate
# S = 84928, p-value = 0.0007948
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3455296 

sum(DEGs_inOtherMeta$ELS_Duan_2025_pval<0.05, na.rm=TRUE)
#[1] 11
sum(DEGs_inOtherMeta$ELS_Duan_2025_FDR<0.05, na.rm=TRUE)
#[1] 0

###################################

#Comparison with Cosette's Sleep Deprivation Results

pdf("Scatterplot_AcuteSleepDep_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(AcuteSleepDep_Rhoads_Log2FC_estimate~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta, xlab="Chronic Stress: Log2FC", ylab="Acute Sleep Deprivation: Log2FC")

points(AcuteSleepDep_Rhoads_Log2FC_estimate~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta[DEGs_inOtherMeta$AcuteSleepDep_Rhoads_pval<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(AcuteSleepDep_Rhoads_Log2FC_estimate~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

#linear regression output for that relationship:
summary.lm(TrendLine)
rm(TrendLine)

# Call:
#   lm(formula = AcuteSleepDep_Rhoads_Log2FC_estimate ~ ChronicStress_Xiong_Log2FC_estimate, 
#      data = DEGs_inOtherMeta)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.34856 -0.05181  0.01210  0.06481  0.20722 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         0.006886   0.009960   0.691    0.491    
# ChronicStress_Xiong_Log2FC_estimate 0.584996   0.055194  10.599   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09523 on 116 degrees of freedom
# (15 observations deleted due to missingness)
# Multiple R-squared:  0.492,	Adjusted R-squared:  0.4876 
# F-statistic: 112.3 on 1 and 116 DF,  p-value: < 2.2e-16

#Calculating the Pearson rank correlation:
cor.test(DEGs_inOtherMeta$AcuteSleepDep_Rhoads_Log2FC_estimate, DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate, method="pearson")

# Pearson's product-moment correlation
# 
# data:  DEGs_inOtherMeta$AcuteSleepDep_Rhoads_Log2FC_estimate and DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate
# t = 10.599, df = 116, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5962473 0.7829077
# sample estimates:
#       cor 
# 0.7014105 

#Calculating the Spearman rank correlation:
cor.test(DEGs_inOtherMeta$AcuteSleepDep_Rhoads_Log2FC_estimate, DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate, method="spearman")

# Spearman's rank correlation rho
# 
# data:  DEGs_inOtherMeta$AcuteSleepDep_Rhoads_Log2FC_estimate and DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate
# S = 172560, p-value = 4.271e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3698027


sum(DEGs_inOtherMeta$AcuteSleepDep_Rhoads_pval<0.05, na.rm=TRUE)
#[1] 17
sum(DEGs_inOtherMeta$AcuteSleepDep_Rhoads_FDR<0.05, na.rm=TRUE)
#[1] 1

#######################

#Comparison with Chronic Corticosterone:


pdf("Scatterplot_ChronicCort_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(ChronicDailyCorticosterone_Jaszczyk_2025_logFC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta, xlab="Chronic Stress: Log2FC", ylab="Chronic Corticosterone: Log2FC")

points(ChronicDailyCorticosterone_Jaszczyk_2025_logFC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta[DEGs_inOtherMeta$ChronicDailyCorticosterone_Jaszczyk_2025_PValue<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(ChronicDailyCorticosterone_Jaszczyk_2025_logFC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

#linear regression output for that relationship:
summary.lm(TrendLine)
rm(TrendLine)

# Call:
#   lm(formula = ChronicDailyCorticosterone_Jaszczyk_2025_logFC ~ 
#        ChronicStress_Xiong_Log2FC_estimate, data = DEGs_inOtherMeta)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.33486 -0.09452  0.02909  0.11416  1.04818 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                         -0.03991    0.02955  -1.351    0.180
# ChronicStress_Xiong_Log2FC_estimate  0.22024    0.16269   1.354    0.179
# 
# Residual standard error: 0.266 on 108 degrees of freedom
# (23 observations deleted due to missingness)
# Multiple R-squared:  0.01669,	Adjusted R-squared:  0.007581 
# F-statistic: 1.833 on 1 and 108 DF,  p-value: 0.1786

#Calculating the Pearson rank correlation:
cor.test(DEGs_inOtherMeta$ChronicDailyCorticosterone_Jaszczyk_2025_logFC, DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate, method="pearson")

# Pearson's product-moment correlation
# 
# data:  DEGs_inOtherMeta$ChronicDailyCorticosterone_Jaszczyk_2025_logFC and DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate
# t = 1.3538, df = 108, p-value = 0.1786
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.05950557  0.30894420
# sample estimates:
#       cor 
# 0.1291752

#Calculating the Spearman rank correlation:
cor.test(DEGs_inOtherMeta$ChronicDailyCorticosterone_Jaszczyk_2025_logFC, DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate, method="spearman")

# Spearman's rank correlation rho
# 
# data:  DEGs_inOtherMeta$ChronicDailyCorticosterone_Jaszczyk_2025_logFC and DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate
# S = 153319, p-value = 0.00103
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3087989 


sum(DEGs_inOtherMeta$ChronicDailyCorticosterone_Jaszczyk_2025_PValue<0.05, na.rm=TRUE)
#[1] 48
sum(DEGs_inOtherMeta$ChronicDailyCorticosterone_Jaszczyk_2025_FDR<0.05, na.rm=TRUE)
#[1] 29

##############################

#Comparison with 10 days CSDS Meta-Analysis:


pdf("Scatterplot_10DaysCSDS_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(CSDS_10Days_Reshetnikov_2022_log2FoldChange~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta, xlab="Chronic Stress: Log2FC", ylab="CSDS 10 Days: Log2FC")

points(CSDS_10Days_Reshetnikov_2022_log2FoldChange~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta[DEGs_inOtherMeta$CSDS_10Days_Reshetnikov_2022_CSDS_pvalue<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(CSDS_10Days_Reshetnikov_2022_log2FoldChange~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

#linear regression output for that relationship:
summary.lm(TrendLine)
rm(TrendLine)

# Call:
#   lm(formula = CSDS_10Days_Reshetnikov_2022_log2FoldChange ~ ChronicStress_Xiong_Log2FC_estimate, 
#      data = DEGs_inOtherMeta)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.24172 -0.03670 -0.02127  0.06779  0.22696 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                         0.0009958  0.0304633   0.033 0.974264
# ChronicStress_Xiong_Log2FC_estimate 0.6413225  0.1357337   4.725 0.000147
# 
# (Intercept)                            
# ChronicStress_Xiong_Log2FC_estimate ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1141 on 19 degrees of freedom
# (112 observations deleted due to missingness)
# Multiple R-squared:  0.5402,	Adjusted R-squared:  0.516 
# F-statistic: 22.32 on 1 and 19 DF,  p-value: 0.0001474


#Calculating the Pearson rank correlation:
cor.test(DEGs_inOtherMeta$CSDS_10Days_Reshetnikov_2022_log2FoldChange, DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate, method="pearson")

# Pearson's product-moment correlation
# 
# data:  DEGs_inOtherMeta$CSDS_10Days_Reshetnikov_2022_log2FoldChange and DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate
# t = 4.7249, df = 19, p-value = 0.0001474
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4442745 0.8856711
# sample estimates:
#       cor 
# 0.7349979 

#Calculating the Spearman rank correlation:
cor.test(DEGs_inOtherMeta$CSDS_10Days_Reshetnikov_2022_log2FoldChange, DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate, method="spearman")

# Spearman's rank correlation rho
# 
# data:  DEGs_inOtherMeta$CSDS_10Days_Reshetnikov_2022_log2FoldChange and DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate
# S = 338, p-value = 4.316e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7805195 

sum(DEGs_inOtherMeta$CSDS_10Days_Reshetnikov_2022_CSDS_pvalue<0.05, na.rm=TRUE)
#[1] 21
sum(DEGs_inOtherMeta$CSDS_10Days_Reshetnikov_2022_padj<0.05, na.rm=TRUE)
#[1] 0

##################################

#Comparison with 30 Days CSDS:


pdf("Scatterplot_30DaysCSDS_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(CSDS_30Days_Reshetnikov_2022_log2FoldChange~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta, xlab="Chronic Stress: Log2FC", ylab="CSDS 30 Days: Log2FC")

points(CSDS_30Days_Reshetnikov_2022_log2FoldChange~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta[DEGs_inOtherMeta$CSDS_30Days_Reshetnikov_2022pvalue<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(CSDS_30Days_Reshetnikov_2022_log2FoldChange~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_inOtherMeta) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

#linear regression output for that relationship:
summary.lm(TrendLine)
rm(TrendLine)

# Call:
#   lm(formula = CSDS_30Days_Reshetnikov_2022_log2FoldChange ~ ChronicStress_Xiong_Log2FC_estimate, 
#      data = DEGs_inOtherMeta)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.2357 -0.1797 -0.1517  0.1821  0.4425 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                         -0.02141    0.06958  -0.308    0.764
# ChronicStress_Xiong_Log2FC_estimate -0.03424    0.47525  -0.072    0.944
# 
# Residual standard error: 0.2346 on 12 degrees of freedom
# (119 observations deleted due to missingness)
# Multiple R-squared:  0.0004323,	Adjusted R-squared:  -0.08287 
# F-statistic: 0.005189 on 1 and 12 DF,  p-value: 0.9438

#Calculating the Pearson rank correlation:
cor.test(DEGs_inOtherMeta$CSDS_30Days_Reshetnikov_2022_log2FoldChange, DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate, method="pearson")

# data:  DEGs_inOtherMeta$CSDS_30Days_Reshetnikov_2022_log2FoldChange and DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate
# t = -0.072037, df = 12, p-value = 0.9438
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.5453544  0.5154750
# sample estimates:
#   cor 
# -0.02079084 

#Calculating the Spearman rank correlation:
cor.test(DEGs_inOtherMeta$CSDS_30Days_Reshetnikov_2022_log2FoldChange, DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate, method="spearman")

# Spearman's rank correlation rho
# 
# data:  DEGs_inOtherMeta$CSDS_30Days_Reshetnikov_2022_log2FoldChange and DEGs_inOtherMeta$ChronicStress_Xiong_Log2FC_estimate
# S = 424, p-value = 0.8201
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.06813187 

sum(DEGs_inOtherMeta$CSDS_30Days_Reshetnikov_2022pvalue<0.05, na.rm=TRUE)
#[1] 14
sum(DEGs_inOtherMeta$CSDS_30Days_Reshetnikov_2022_padj<0.05, na.rm=TRUE)
#[1] 4

####################

#Maybe we could make some boxplots to compare with the categorical variables?

#I'm not sure that the categories other than Down and Up are interesting to us - may need to condense down to those

table(DEGs_inOtherMeta$ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction)

pdf("Boxplot_AcuteCorticosterone_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=7)

boxplot(ChronicStress_Xiong_Log2FC_estimate~ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction, data=DEGs_inOtherMeta, col="darkgrey", ylab="Chronic Stress: Log2FC")

stripchart(ChronicStress_Xiong_Log2FC_estimate~ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction, data=DEGs_inOtherMeta, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')

dev.off()

#I should probably only compare the Ups and Downs:

pdf("Boxplot_AcuteCorticosterone_vs_ChronicStress_Xiong_Log2FC_DEGs_UpDown.pdf", height=5, width=4)

boxplot(ChronicStress_Xiong_Log2FC_estimate~ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction, data=DEGs_inOtherMeta[(DEGs_inOtherMeta$ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction=="Up"|DEGs_inOtherMeta$ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction=="Down"),], col="darkgrey", ylab="Chronic Stress: Log2FC")

stripchart(ChronicStress_Xiong_Log2FC_estimate~ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction, data=DEGs_inOtherMeta[(DEGs_inOtherMeta$ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction=="Up"|DEGs_inOtherMeta$ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction=="Down"),], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')

dev.off()

#defaults to a two-sided unpaired Welch's t-test:
t.test(ChronicStress_Xiong_Log2FC_estimate~ ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction, data=DEGs_inOtherMeta[(DEGs_inOtherMeta$ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction=="Up"|DEGs_inOtherMeta$ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction=="Down"),]) 

# Welch Two Sample t-test
# 
# data:  ChronicStress_Xiong_Log2FC_estimate by ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction
# t = -1.4506, df = 51.36, p-value = 0.153
# alternative hypothesis: true difference in means between group Down and group Up is not equal to 0
# 95 percent confidence interval:
#   -0.14505874  0.02334987
# sample estimates:
#   mean in group Down   mean in group Up 
# -0.1198534         -0.0589990 

###################

#Responsive to Acute Stress:


table(DEGs_inOtherMeta$Responsivity.to.stress_Stankiewicz_2022_AcuteStressResponsive)
# Down     Up   Up&Down 
# 120       1       9       3 

pdf("Boxplot_AcuteStress_vs_ChronicStress_Xiong_Log2FC_DEGs_UpDown.pdf", height=5, width=4)

boxplot(ChronicStress_Xiong_Log2FC_estimate~Responsivity.to.stress_Stankiewicz_2022_AcuteStressResponsive, data=DEGs_inOtherMeta[(DEGs_inOtherMeta$Responsivity.to.stress_Stankiewicz_2022_AcuteStressResponsive=="Up"|DEGs_inOtherMeta$Responsivity.to.stress_Stankiewicz_2022_AcuteStressResponsive=="Down"),], col="darkgrey", ylab="Chronic Stress: Log2FC")

stripchart(ChronicStress_Xiong_Log2FC_estimate~Responsivity.to.stress_Stankiewicz_2022_AcuteStressResponsive, data=DEGs_inOtherMeta[(DEGs_inOtherMeta$Responsivity.to.stress_Stankiewicz_2022_AcuteStressResponsive=="Up"|DEGs_inOtherMeta$Responsivity.to.stress_Stankiewicz_2022_AcuteStressResponsive=="Down"),], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')

dev.off()

#############

#Medium Stress Responsive:

table(DEGs_inOtherMeta$Responsivity.to.stress_Stankiewicz_2022_MediumStressResponsive)
#         Down      Up Up&Down 
# 121       4       4       4 

pdf("Boxplot_MediumStressResponsive_vs_ChronicStress_Xiong_Log2FC_DEGs_UpDown.pdf", height=5, width=4)

boxplot(ChronicStress_Xiong_Log2FC_estimate~Responsivity.to.stress_Stankiewicz_2022_MediumStressResponsive, data=DEGs_inOtherMeta[(DEGs_inOtherMeta$Responsivity.to.stress_Stankiewicz_2022_MediumStressResponsive=="Up"|DEGs_inOtherMeta$Responsivity.to.stress_Stankiewicz_2022_MediumStressResponsive=="Down"),], col="darkgrey", ylab="Chronic Stress: Log2FC")

stripchart(ChronicStress_Xiong_Log2FC_estimate~Responsivity.to.stress_Stankiewicz_2022_MediumStressResponsive, data=DEGs_inOtherMeta[(DEGs_inOtherMeta$Responsivity.to.stress_Stankiewicz_2022_MediumStressResponsive=="Up"|DEGs_inOtherMeta$Responsivity.to.stress_Stankiewicz_2022_MediumStressResponsive=="Down"),], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')

dev.off()

#defaults to a two-sided unpaired Welch's t-test:
t.test(ChronicStress_Xiong_Log2FC_estimate~ Responsivity.to.stress_Stankiewicz_2022_MediumStressResponsive, data=DEGs_inOtherMeta[(DEGs_inOtherMeta$Responsivity.to.stress_Stankiewicz_2022_MediumStressResponsive=="Up"|DEGs_inOtherMeta$Responsivity.to.stress_Stankiewicz_2022_MediumStressResponsive=="Down"),]) 

###########

#Prolonged Stress - probably too small to work with:

table(DEGs_inOtherMeta$Responsivity.to.stress_Stankiewicz_2022_ProlongedStressResponsive)
#         Down      Up    Up&Down 
# 117       2       1      13


######################

#Reading in a few human studies that might serve as an interesting point of comparison too:

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2022_JinglinXiong_ChronicStress_FrontalCortex/Final Product_2025_Corrected/Reference Articles/MetaAnalysis_Suicide_PFC_Piras2022_1-s2.0-S0924977X21016631-mmc1")

list.files()

Suicide<-read.csv("MetaAnalysis_Suicide_PFC_Piras2022_forR.csv", header=TRUE, stringsAsFactors = FALSE)

str(Suicide)
#'data.frame':	299 obs. of  8 variables:
#Note: Some of the genes in this study show up multiple times in this spreadsheet because there were 4 different meta-analyses

table(table(Suicide$Symbol))
# 1   2   3 
# 148  68   5

#It may make sense to reorganize the data.frame into a data.frame with each meta-analysis as its own columns

table(Suicide$Meta.Analysis)
# Meta-analysis A (OFC)      Meta-analysis B (DLPFC/ PFC) 
# 60                               109 
# Meta-analysis C (DLPFC) Meta-analysis D (OFC, DLPFC, PFC) 
# 126                                 4

Suicide_OFC<-Suicide[Suicide$Meta.Analysis=="Meta-analysis A (OFC)",]
colnames(Suicide_OFC)
#I need to rename most of the columns before joining
colnames(Suicide_OFC)[c(1,4:8)]<-paste("Suicide_OFC", colnames(Suicide_OFC)[c(1,4:8)], sep="_")

Suicide_PFC<-Suicide[Suicide$Meta.Analysis=="Meta-analysis B (DLPFC/ PFC)",]

colnames(Suicide_PFC)[c(1,4:8)]<-paste("Suicide_PFC", colnames(Suicide_PFC)[c(1,4:8)], sep="_")

Suicide_DLPFC<-Suicide[Suicide$Meta.Analysis=="Meta-analysis C (DLPFC)",]

colnames(Suicide_DLPFC)[c(1,4:8)]<-paste("Suicide_DlPFC", colnames(Suicide_DLPFC)[c(1,4:8)], sep="_")

Suicide_FC<-Suicide[Suicide$Meta.Analysis=="Meta-analysis D (OFC, DLPFC, PFC)",]

colnames(Suicide_FC)[c(1,4:8)]<-paste("Suicide_FC", colnames(Suicide_FC)[c(1,4:8)], sep="_")

library(plyr)

SuicideDF<-join_all(list(Suicide_OFC, Suicide_DLPFC, Suicide_PFC, Suicide_FC), by="Symbol", type="full")

str(SuicideDF)
#'data.frame':	221 obs. of  26 variables:

#Gandal 2018's meta-analyses are also really high quality

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2022_JinglinXiong_ChronicStress_FrontalCortex/Final Product_2025_Corrected/Reference Articles/Gandal_2018_PsychiatricDisorder_Cortex")

list.files()

Gandal_Psych_RNASeq<-read.csv("GandalFrontalCortex_RNASeqMeta_wSymbols_NoNA_NoMultiMapped.csv", header=TRUE, stringsAsFactors = FALSE) 

Gandal_Psych_Microarray<-read.csv("Gandal_2018_PsychTranscriptomeVsGenetics_Suppl_TableS1_MicroarrayMetaAnalysis.csv", header=TRUE, stringsAsFactors = FALSE) 

colnames(Gandal_Psych_Microarray)

colnames(Gandal_Psych_Microarray)[c(8:25)]<-paste("Gandal_Microarray", colnames(Gandal_Psych_Microarray)[c(8:25)], sep="_")

colnames(Gandal_Psych_RNASeq)
colnames(Gandal_Psych_RNASeq)[c(3:10)]<-paste("Gandal_RNASeq", colnames(Gandal_Psych_RNASeq)[c(3:10)], sep="_")
colnames(Gandal_Psych_RNASeq)[2]<-"ensembl_gene_id"

Gandal_Psych<-join(Gandal_Psych_Microarray, Gandal_Psych_RNASeq, by="ensembl_gene_id", type="full")

str(Gandal_Psych)
#'data.frame':	26992 obs. of  36 variables:


str(SuicideDF)
#includes Entrez.ID as an identifier

#Which Entrez ID column should I join by?
sum(SuicideDF$Entrez.ID %in% Gandal_Psych$entrezgene)
#[1] 217
sum(SuicideDF$Entrez.ID %in% Gandal_Psych$ENTREZID)
#[1] 213

#hmm... do they capture different information?
sum((SuicideDF$Entrez.ID %in% Gandal_Psych$entrezgene)==(SuicideDF$Entrez.ID %in% Gandal_Psych$ENTREZID))
#[1] 209
#Apparently slightly different. Sigh.
sum(is.na(Gandal_Psych$entrezgene))
#[1] 5552
sum(is.na(Gandal_Psych$ENTREZID))
#[1] 12783 #much bigger

Gandal_Psych$Entrez.ID<-Gandal_Psych$entrezgene
Gandal_Psych$Entrez.ID[is.na(Gandal_Psych$Entrez.ID)]<-Gandal_Psych$ENTREZID[is.na(Gandal_Psych$Entrez.ID)]

sum(is.na(Gandal_Psych$Entrez.ID))
#[1] 5087 #huh

sum(SuicideDF$Entrez.ID %in% Gandal_Psych$Entrez.ID)
#[1] 220
#All but one - pretty good.

Gandal_Psych_Suicide<-join(Gandal_Psych,SuicideDF, by="Entrez.ID", type="full")
str(Gandal_Psych_Suicide)

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2022_JinglinXiong_ChronicStress_FrontalCortex/Final Product_2025_Corrected/Reference Articles")

write.csv(Gandal_Psych_Suicide, "CombinedDF_Gandal_Psych_Piras_Suicide.csv")

#Now I need the ortholog info


list.files()
Orthologs<-read.csv("HCOP_Orthologs_Best.csv", header=TRUE, stringsAsFactors = FALSE)
colnames(Orthologs)

colnames(Gandal_Psych_Suicide)
#[37] "Entrez.ID"
colnames(Gandal_Psych_Suicide)[37]<-"human_entrez_gene"

sum(duplicated(Gandal_Psych_Suicide$human_entrez_gene))
#[1] 5413
sum(is.na(Gandal_Psych_Suicide$human_entrez_gene))
#[1] 5087
sum(is.na(Gandal_Psych_Suicide$hgnc_symbol))

sum((is.na(Gandal_Psych_Suicide$human_entrez_gene)==TRUE)&(is.na(Gandal_Psych_Suicide$hgnc_symbol)==TRUE))
#[1] 345
#Oh interesting.

sum(is.na(Orthologs$human_entrez_gene))
#[1] 5185
sum(is.na(Orthologs$human_symbol))
#[1] 5185
sum(is.na(Orthologs$human_ensembl_gene))
#[1] 5185
#Ok, let's trim out the NAs before joining because I think they are causing problems

Orthologs_NoNA<-Orthologs[is.na(Orthologs$human_ensembl_gene)==FALSE,]
str(Orthologs_NoNA)
#'data.frame':	17274 obs. of  18 variables:

write.csv(Orthologs_NoNA, "Orthologs_Best_noHumanNA.csv")

setwd("~/Library/CloudStorage/GoogleDrive-hagenaue@umich.edu/My Drive/BrainAlchemyProject/ProjectFolders/2022_JinglinXiong_ChronicStress_FrontalCortex/Final Product_2025_Corrected/Reference Articles")

Gandal_Psych_Suicide_wOrthologs<-join(Gandal_Psych_Suicide, Orthologs_NoNA, by="human_entrez_gene", type="left")

write.csv(Gandal_Psych_Suicide_wOrthologs, "Gandal_Psych_Suicide_wOrthologs.csv")

list.files()

DEGs_vs_OtherMetaAnalyses<-read.csv("DEGs_vs_OtherMetaAnalyses.csv", header=TRUE, stringsAsFactors = FALSE)

#I'm going to add the acute sleep deprivation results from the large cortical study too:
SleepDep_GSE114845<-read.delim("Limma_Results_SleepDeprivation_GSE114845.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE) 

str(SleepDep_GSE114845)
colnames(SleepDep_GSE114845)

str(DEGs_vs_OtherMetaAnalyses)
colnames(SleepDep_GSE114845)[13]<-"GeneSymbol"



DEGs_vs_OtherMetaAnalysesSleepDep<-join(DEGs_vs_OtherMetaAnalyses, SleepDep_GSE114845, by="GeneSymbol", type="left")

str(DEGs_vs_OtherMetaAnalysesSleepDep)
write.csv(DEGs_vs_OtherMetaAnalysesSleepDep, "DEGs_vs_OtherMetaAnalyses_wSleepDep_GSE114845.csv")

plot(DEGs_vs_OtherMetaAnalysesSleepDep$ChronicStress_Xiong_Log2FC_estimate~DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived)

summary.lm(lm(DEGs_vs_OtherMetaAnalysesSleepDep$ChronicStress_Xiong_Log2FC_estimate~DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived))

# Call:
#   lm(formula = DEGs_vs_OtherMetaAnalysesSleepDep$ChronicStress_Xiong_Log2FC_estimate ~ 
#        DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.29074 -0.09005 -0.01816  0.05190  0.77768 
# 
# Coefficients:
#   Estimate
# (Intercept)                                                    -0.07170
# DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived -0.16863
# Std. Error
# (Intercept)                                                       0.01477
# DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived    0.03093
# t value Pr(>|t|)
# (Intercept)                                                     -4.854 3.87e-06
# DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived  -5.452 2.92e-07
# 
# (Intercept)                                                    ***
#   DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.158 on 114 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.2068,	Adjusted R-squared:  0.1999 
# F-statistic: 29.73 on 1 and 114 DF,  p-value: 2.92e-07

#sanity check:
summary.lm(lm(DEGs_vs_OtherMetaAnalysesSleepDep$AcuteSleepDep_Rhoads_Log2FC_estimate~DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived))

# Call:
#   lm(formula = DEGs_vs_OtherMetaAnalysesSleepDep$AcuteSleepDep_Rhoads_Log2FC_estimate ~ 
#        DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.38378 -0.03548  0.01245  0.05157  0.48329 
# 
# Coefficients:
#   Estimate
# (Intercept)                                                    -0.03468
# DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived -0.17584
# Std. Error
# (Intercept)                                                       0.01021
# DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived    0.02095
# t value Pr(>|t|)
# (Intercept)                                                     -3.396 0.000954
# DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived  -8.395 1.89e-13
# 
# (Intercept)                                                    ***
#   DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1067 on 109 degrees of freedom
# (22 observations deleted due to missingness)
# Multiple R-squared:  0.3927,	Adjusted R-squared:  0.3871 
# F-statistic: 70.47 on 1 and 109 DF,  p-value: 1.888e-13

#Did I output GSE114845 with the log2fc reversed? I vaguely remember that. Let me double-check...
#Nope, it just seems to go in the opposite direction from Cosette's results for the chronic stress genes. Very interesting.

#I should make a nicer plot of it:


pdf("Scatterplot_GSE114845SleepDep_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(Coef.ConditionSleep.Deprived~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalysesSleepDep, xlab="Chronic Stress: Log2FC", ylab="Acute Sleep Deprivation: Log2FC")

points(Coef.ConditionSleep.Deprived~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalysesSleepDep[DEGs_vs_OtherMetaAnalysesSleepDep$P.value.ConditionSleep.Deprived<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(Coef.ConditionSleep.Deprived~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalysesSleepDep) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = Coef.ConditionSleep.Deprived ~ ChronicStress_Xiong_Log2FC_estimate, 
#      data = DEGs_vs_OtherMetaAnalysesSleepDep)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.93484 -0.27361 -0.02743  0.13618  2.19039 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -0.04366    0.04357  -1.002    0.318    
# ChronicStress_Xiong_Log2FC_estimate -1.22651    0.22496  -5.452 2.92e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4261 on 114 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.2068,	Adjusted R-squared:  0.1999 
# F-statistic: 29.73 on 1 and 114 DF,  p-value: 2.92e-07


cor.test(DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived,DEGs_vs_OtherMetaAnalysesSleepDep$ChronicStress_Xiong_Log2FC_estimate, method="pearson")
# data:  DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived and DEGs_vs_OtherMetaAnalysesSleepDep$ChronicStress_Xiong_Log2FC_estimate
# t = -5.4522, df = 114, p-value = 2.92e-07
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.5883175 -0.2970976
# sample estimates:
#   cor 
# -0.4547806

cor.test(DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived,DEGs_vs_OtherMetaAnalysesSleepDep$ChronicStress_Xiong_Log2FC_estimate, method="spearman")
# data:  DEGs_vs_OtherMetaAnalysesSleepDep$Coef.ConditionSleep.Deprived and DEGs_vs_OtherMetaAnalysesSleepDep$ChronicStress_Xiong_Log2FC_estimate
# S = 264702, p-value = 0.8513
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.01757583

sum(DEGs_vs_OtherMetaAnalysesSleepDep$P.value.ConditionSleep.Deprived<0.05, na.rm=TRUE)
#[1] 81

rm(TrendLine)

colnames(DEGs_vs_OtherMetaAnalysesSleepDep)
colnames(Gandal_Psych_Suicide_wOrthologs)
#[74] "mouse_symbol"
colnames(Gandal_Psych_Suicide_wOrthologs)[74]<-"GeneSymbol" 

DEGs_vs_OtherMetaAnalyses_HumanPsych<-join(DEGs_vs_OtherMetaAnalysesSleepDep, Gandal_Psych_Suicide_wOrthologs, by="GeneSymbol", type="left")

write.csv(DEGs_vs_OtherMetaAnalyses_HumanPsych, "DEGs_vs_OtherMetaAnalyses_HumanPsych.csv")

#How many suicide-related genes are in our DEG list?
DEGs_vs_OtherMetaAnalyses_HumanPsych$GeneSymbol[(is.na(DEGs_vs_OtherMetaAnalyses_HumanPsych$Suicide_OFC_Fisher_log2.FC..averaged.)==FALSE)|(is.na(DEGs_vs_OtherMetaAnalyses_HumanPsych$Suicide_DlPFC_Fisher_log2.FC..averaged.)==FALSE)|(is.na(DEGs_vs_OtherMetaAnalyses_HumanPsych$Suicide_FC_Fisher_log2.FC..averaged.)==FALSE)|(is.na(DEGs_vs_OtherMetaAnalyses_HumanPsych$Suicide_PFC_Fisher_log2.FC..averaged.)==FALSE) ]
#[1] "Aqp4"    "Csrp1"   "Elovl5"  "Kif1c"   "Mid1ip1" "Sox9"    "Tcf7l2" 

#How about PTSD?
table(DEGs_vs_OtherMetaAnalyses_HumanPsych$Responsivity.to.stress_Stankiewicz_2022_PTSD)
#       Down 
# 132    3 

DEGs_vs_OtherMetaAnalyses_HumanPsych$GeneSymbol[DEGs_vs_OtherMetaAnalyses_HumanPsych$Responsivity.to.stress_Stankiewicz_2022_PTSD=="Down"]
#[1] "Fa2h"  "Lpar1" "Mtus1"


#columns of interest for plotting:
Gandal_Microarray_MDD.beta_log2FC
Gandal_Microarray_SCZ.beta_log2FC
Gandal_Microarray_BD.beta_log2FC
Gandal_Microarray_AAD.beta_log2FC
Gandal_RNASeq_BD.logFC
Gandal_RNASeq_SCZ.logFC
#suicide looks like it has too few overlapping data points

#Let's make some scatterplots:

pdf("Scatterplot_Gandal_MDD_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(Gandal_Microarray_MDD.beta_log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych, xlab="Chronic Stress: Log2FC", ylab="MDD: Log2FC")

points(Gandal_Microarray_MDD.beta_log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych[DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_MDD.P.value<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(Gandal_Microarray_MDD.beta_log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

#linear regression output for that relationship:
summary.lm(TrendLine)
# Call:
#   lm(formula = Gandal_Microarray_MDD.beta_log2FC ~ ChronicStress_Xiong_Log2FC_estimate, 
#      data = DEGs_vs_OtherMetaAnalyses_HumanPsych)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.189635 -0.025028 -0.000018  0.029462  0.106454 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                         -0.004648   0.005392  -0.862   0.3908   
# ChronicStress_Xiong_Log2FC_estimate  0.097402   0.028977   3.361   0.0011 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04545 on 100 degrees of freedom
# (33 observations deleted due to missingness)
# Multiple R-squared:  0.1015,	Adjusted R-squared:  0.09253 
# F-statistic:  11.3 on 1 and 100 DF,  p-value: 0.0011

rm(TrendLine)

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_MDD.beta_log2FC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="pearson")

# Pearson's product-moment correlation
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_MDD.beta_log2FC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# t = 3.3614, df = 100, p-value = 0.0011
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1323430 0.4831547
# sample estimates:
#       cor 
# 0.3186175 

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_MDD.beta_log2FC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="spearman")

# Spearman's rank correlation rho
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_MDD.beta_log2FC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# S = 170064, p-value = 0.7018
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.03837737 

sum(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_MDD.P.value<0.05, na.rm=TRUE)
#[1] 20

pdf("Scatterplot_Gandal_MicroarraySchiz_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(Gandal_Microarray_SCZ.beta_log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych, xlab="Chronic Stress: Log2FC", ylab="Schizophrenia: Log2FC")

points(Gandal_Microarray_SCZ.beta_log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych[DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_SCZ.P.value<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(Gandal_Microarray_SCZ.beta_log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

#linear regression output for that relationship:
summary.lm(TrendLine)

# Call:
#   lm(formula = Gandal_Microarray_SCZ.beta_log2FC ~ ChronicStress_Xiong_Log2FC_estimate, 
#      data = DEGs_vs_OtherMetaAnalyses_HumanPsych)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.28391 -0.07386 -0.00475  0.05025  0.44634 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          0.04197    0.01616   2.598 0.011373 *  
#   ChronicStress_Xiong_Log2FC_estimate  0.30166    0.08059   3.743 0.000363 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1122 on 72 degrees of freedom
# (61 observations deleted due to missingness)
# Multiple R-squared:  0.1629,	Adjusted R-squared:  0.1513 
# F-statistic: 14.01 on 1 and 72 DF,  p-value: 0.0003627

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_SCZ.beta_log2FC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="pearson")

# Pearson's product-moment correlation
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_SCZ.beta_log2FC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# t = 3.7432, df = 72, p-value = 0.0003627
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1929013 0.5787351
# sample estimates:
#       cor 
# 0.4036099

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_SCZ.beta_log2FC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="spearman")

# Spearman's rank correlation rho
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_SCZ.beta_log2FC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# S = 59514, p-value = 0.314
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1186411 

sum(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_SCZ.P.value<0.05, na.rm=TRUE)

pdf("Scatterplot_Gandal_MicroarrayBP_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(Gandal_Microarray_BD.beta_log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych, xlab="Chronic Stress: Log2FC", ylab="BPD: Log2FC")

points(Gandal_Microarray_BD.beta_log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych[DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_BD.P.value<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(Gandal_Microarray_BD.beta_log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = Gandal_Microarray_BD.beta_log2FC ~ ChronicStress_Xiong_Log2FC_estimate, 
#      data = DEGs_vs_OtherMetaAnalyses_HumanPsych)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.26789 -0.05769 -0.00176  0.03571  0.24142 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                          0.02435    0.01310   1.858   0.0672 .
# ChronicStress_Xiong_Log2FC_estimate  0.12225    0.06537   1.870   0.0655 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09098 on 72 degrees of freedom
# (61 observations deleted due to missingness)
# Multiple R-squared:  0.04632,	Adjusted R-squared:  0.03307 
# F-statistic: 3.497 on 1 and 72 DF,  p-value: 0.06554

rm(TrendLine)

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_BD.beta_log2FC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="pearson")

# Pearson's product-moment correlation
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_BD.beta_log2FC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# t = 1.87, df = 72, p-value = 0.06554
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.01396338  0.42292233
# sample estimates:
#      cor 
# 0.215222 

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_BD.beta_log2FC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="spearman")

# Spearman's rank correlation rho
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_BD.beta_log2FC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# S = 72522, p-value = 0.5309
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.07400441 

sum(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_BD.P.value<0.05, na.rm=TRUE)

pdf("Scatterplot_Gandal_AlcoholAbuse_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(Gandal_Microarray_AAD.beta_log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych, xlab="Chronic Stress: Log2FC", ylab="AAD: Log2FC")

points(Gandal_Microarray_AAD.beta_log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych[DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_AAD.P.value<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(Gandal_Microarray_AAD.beta_log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = Gandal_Microarray_AAD.beta_log2FC ~ ChronicStress_Xiong_Log2FC_estimate, 
#      data = DEGs_vs_OtherMetaAnalyses_HumanPsych)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.64203 -0.07230 -0.00566  0.10405  0.47028 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                          0.02393    0.02385   1.004   0.3181  
# ChronicStress_Xiong_Log2FC_estimate  0.26802    0.12661   2.117   0.0369 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.192 on 96 degrees of freedom
# (37 observations deleted due to missingness)
# Multiple R-squared:  0.04459,	Adjusted R-squared:  0.03464 
# F-statistic: 4.481 on 1 and 96 DF,  p-value: 0.03686

rm(TrendLine)

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_AAD.beta_log2FC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="pearson")

# Pearson's product-moment correlation
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_AAD.beta_log2FC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# t = 2.1168, df = 96, p-value = 0.03686
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.01330979 0.39312142
# sample estimates:
#       cor 
# 0.2111729 

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_AAD.beta_log2FC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="spearman")

# Spearman's rank correlation rho
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_AAD.beta_log2FC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# S = 133396, p-value = 0.1417
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1495279

sum(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_Microarray_AAD.P.value<0.05, na.rm=TRUE)
#[1] 18


pdf("Scatterplot_Gandal_RNASeqBPD_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(Gandal_RNASeq_BD.logFC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych, xlab="Chronic Stress: Log2FC", ylab="BPD: Log2FC")

points(Gandal_RNASeq_BD.logFC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych[DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_BD.P.Value<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(Gandal_RNASeq_BD.logFC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = Gandal_RNASeq_BD.logFC ~ ChronicStress_Xiong_Log2FC_estimate, 
#      data = DEGs_vs_OtherMetaAnalyses_HumanPsych)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.32300 -0.04147  0.00721  0.05197  0.22473 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          0.01927    0.01236   1.560    0.123    
# ChronicStress_Xiong_Log2FC_estimate  0.41845    0.06506   6.432  7.6e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08964 on 83 degrees of freedom
# (50 observations deleted due to missingness)
# Multiple R-squared:  0.3327,	Adjusted R-squared:  0.3246 
# F-statistic: 41.37 on 1 and 83 DF,  p-value: 7.601e-09

rm(TrendLine)

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_BD.logFC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="pearson")

# Pearson's product-moment correlation
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_BD.logFC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# t = 6.4322, df = 83, p-value = 7.601e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4146030 0.7034213
# sample estimates:
#       cor 
# 0.5767632 

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_BD.logFC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="spearman")

# Spearman's rank correlation rho
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_BD.logFC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# S = 91350, p-value = 0.3279
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1073892 

sum(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_BD.P.Value<0.05, na.rm=TRUE)
#[1] 15

pdf("Scatterplot_Gandal_RNASeq_Schiz_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(Gandal_RNASeq_SCZ.logFC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych, xlab="Chronic Stress: Log2FC", ylab="Schizophrenia: Log2FC")

points(Gandal_RNASeq_SCZ.logFC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych[DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_SCZ.P.Value<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(Gandal_RNASeq_SCZ.logFC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = Gandal_RNASeq_SCZ.logFC ~ ChronicStress_Xiong_Log2FC_estimate, 
#      data = DEGs_vs_OtherMetaAnalyses_HumanPsych)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.37641 -0.05874  0.00266  0.06393  0.26367 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          0.02145    0.01515   1.416 0.160438    
# ChronicStress_Xiong_Log2FC_estimate  0.28859    0.07973   3.620 0.000506 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1099 on 83 degrees of freedom
# (50 observations deleted due to missingness)
# Multiple R-squared:  0.1363,	Adjusted R-squared:  0.1259 
# F-statistic:  13.1 on 1 and 83 DF,  p-value: 0.0005058

rm(TrendLine)

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_SCZ.logFC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="pearson")

# Pearson's product-moment correlation
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_SCZ.logFC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# t = 3.6198, df = 83, p-value = 0.0005058
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1694591 0.5398851
# sample estimates:
#       cor 
# 0.3692475

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_SCZ.logFC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="spearman")

# Spearman's rank correlation rho
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_SCZ.logFC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# S = 98992, p-value = 0.7663
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.03271512 

sum(DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_SCZ.P.Value<0.05, na.rm=TRUE)
#[1] 23
#I'm going to output a version of the table that's trimmed down and cleaned up to compare with the other stress findings:

colnames(DEGs_vs_OtherMetaAnalyses_HumanPsych)
# [1] "GeneSymbol"
# [2] "ChronicStress_Xiong_Log2FC_estimate" 
# [4] "ChronicStress_Xiong_pval" 
# [8] "ChronicStress_Xiong_FDR"
# [10] "Responsivity.to.stress_Stankiewicz_2022_AcuteStressResponsive"        
# [11] "Responsivity.to.stress_Stankiewicz_2022_MediumStressResponsive"       
# [12] "Responsivity.to.stress_Stankiewicz_2022_ProlongedStressResponsive"    
# [13] "Responsivity.to.stress_Stankiewicz_2022_PTSD"                         
# [14] "ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_Direction"        
# [15] "ResponsivityToAcuteCorticosterone_12h_Jaszczyk_2023_NumberSigTimepoints"
# [16] "ResponsivityToGlucocorticoids_Juszczak_Stankiewicz2018"
# [18] "ChronicDailyCorticosterone_Jaszczyk_2025_logFC"                         
# [20] "ChronicDailyCorticosterone_Jaszczyk_2025_PValue"                        
# [21] "ChronicDailyCorticosterone_Jaszczyk_2025_FDR" 
# [26] "ELS_Duan_2025_Log2FC_estimate"
# [28] "ELS_Duan_2025_pval"
# [32] "ELS_Duan_2025_FDR"
# [34] "AcuteSleepDep_Rhoads_Log2FC_estimate"
# [36] "AcuteSleepDep_Rhoads_pval"
# [40] "AcuteSleepDep_Rhoads_FDR"
# [41] "CSDS_10Days_Reshetnikov_2022_log2FoldChange"                            
# [42] "CSDS_10Days_Reshetnikov_2022_CSDS_pvalue"                               
# [43] "CSDS_10Days_Reshetnikov_2022_padj"                                      
# [44] "CSDS_30Days_Reshetnikov_2022_log2FoldChange"                            
# [45] "CSDS_30Days_Reshetnikov_2022pvalue"                                     
# [46] "CSDS_30Days_Reshetnikov_2022_padj"                                      
# [47] "ChronicStress_Gururajan2022_log2FoldChange"
# [50] "ChronicStress_Gururajan2022_pvalue"
# [55] "Coef.ConditionSleep.Deprived"
# [59] "P.value.ConditionSleep.Deprived"
# [61] "P.value.adj.ConditionSleep.Deprived" 
# [71] "Gandal_Microarray_ASD.beta_log2FC"                                      
# [72] "Gandal_Microarray_ASD.P.value"                                          
# [73] "Gandal_Microarray_ASD.FDR"                                              
# [74] "Gandal_Microarray_SCZ.beta_log2FC"                                      
# [75] "Gandal_Microarray_SCZ.P.value"                                          
# [76] "Gandal_Microarray_SCZ.FDR"                                              
# [77] "Gandal_Microarray_BD.beta_log2FC"                                       
# [78] "Gandal_Microarray_BD.P.value"                                           
# [79] "Gandal_Microarray_BD.FDR"                                               
# [80] "Gandal_Microarray_MDD.beta_log2FC"                                      
# [81] "Gandal_Microarray_MDD.P.value"                                          
# [82] "Gandal_Microarray_MDD.FDR"                                              
# [83] "Gandal_Microarray_AAD.beta_log2FC"                                      
# [84] "Gandal_Microarray_AAD.P.value"                                          
# [85] "Gandal_Microarray_AAD.FDR" 
# [90] "Gandal_RNASeq_SCZ.logFC"                                                
# [92] "Gandal_RNASeq_SCZ.P.Value"                                              
# [93] "Gandal_RNASeq_SCZ.adj.P.Val"
# [94] "Gandal_RNASeq_BD.logFC"                                                 
# [96] "Gandal_RNASeq_BD.P.Value"                                               
# [97] "Gandal_RNASeq_BD.adj.P.Val" 
# [98] "SYMBOL"
# [100] "human_entrez_gene"
# [103] "Suicide_OFC_GeneMeta_Z"                                                 
# [104] "Suicide_OFC_GeneMeta_FDR"                                               
# [105] "Suicide_OFC_Fisher_log2.FC..averaged."                                  
# [106] "Suicide_OFC_Fisher_p"                                                   
# [107] "Suicide_OFC_Fisher_FDR"                                                 
# [109] "Suicide_DlPFC_GeneMeta_Z"                                               
# [110] "Suicide_DlPFC_GeneMeta_FDR"                                             
# [111] "Suicide_DlPFC_Fisher_log2.FC..averaged."                                
# [112] "Suicide_DlPFC_Fisher_p"                                                 
# [113] "Suicide_DlPFC_Fisher_FDR"                                               
# [115] "Suicide_PFC_GeneMeta_Z"                                                 
# [116] "Suicide_PFC_GeneMeta_FDR"                                               
# [117] "Suicide_PFC_Fisher_log2.FC..averaged."                                  
# [118] "Suicide_PFC_Fisher_p"                                                   
# [119] "Suicide_PFC_Fisher_FDR"                                                 
# [121] "Suicide_FC_GeneMeta_Z"                                                  
# [122] "Suicide_FC_GeneMeta_FDR"                                                
# [123] "Suicide_FC_Fisher_log2.FC..averaged."                                   
# [124] "Suicide_FC_Fisher_p"                                                    
# [125] "Suicide_FC_Fisher_FDR"
# [126] "human_ensembl_gene" 

DEGs_vs_OtherMetaAnalyses_forTable<-DEGs_vs_OtherMetaAnalyses_HumanPsych[c(1:2,4,8,10:16,18,20:21,26,28,32,34,36,40,41:47,50,55,59,61,71:85,90,92:94,96:98,100,103:107,109:113,115:119,121:126)]

write.csv(DEGs_vs_OtherMetaAnalyses_forTable, "DEGs_vs_OtherMetaAnalyses_forTable.csv")

DEGs_vs_OtherMetaAnalyses_forTable_Pretty<-DEGs_vs_OtherMetaAnalyses_forTable

DEGs_vs_OtherMetaAnalyses_forTable_Pretty$ChronicDailyCorticosterone_Jaszczyk_2025_logFC[DEGs_vs_OtherMetaAnalyses_forTable_Pretty$ChronicDailyCorticosterone_Jaszczyk_2025_PValue>0.05]<-NA

DEGs_vs_OtherMetaAnalyses_forTable_Pretty$ELS_Duan_2025_Log2FC_estimate[DEGs_vs_OtherMetaAnalyses_forTable_Pretty$ELS_Duan_2025_pval>0.05]<-NA

DEGs_vs_OtherMetaAnalyses_forTable_Pretty$AcuteSleepDep_Rhoads_Log2FC_estimate[DEGs_vs_OtherMetaAnalyses_forTable_Pretty$AcuteSleepDep_Rhoads_pval>0.05]<-NA

DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Coef.ConditionSleep.Deprived[DEGs_vs_OtherMetaAnalyses_forTable_Pretty$P.value.ConditionSleep.Deprived>0.05]<-NA

DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_Microarray_ASD.beta_log2FC[DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_Microarray_ASD.P.value>0.05]<-NA

DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_Microarray_SCZ.beta_log2FC[DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_Microarray_SCZ.P.value>0.05]<-NA

DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_Microarray_BD.beta_log2FC[DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_Microarray_BD.P.value>0.05]<-NA

DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_Microarray_MDD.beta_log2FC[DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_Microarray_MDD.P.value>0.05]<-NA

DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_Microarray_AAD.beta_log2FC[DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_Microarray_AAD.P.value>0.05]<-NA

DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_RNASeq_SCZ.logFC[DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_RNASeq_SCZ.P.Value>0.05]<-NA

DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_RNASeq_BD.logFC[DEGs_vs_OtherMetaAnalyses_forTable_Pretty$Gandal_RNASeq_BD.P.Value>0.05]<-NA

write.csv(DEGs_vs_OtherMetaAnalyses_forTable_Pretty, "DEGs_vs_OtherMetaAnalyses_forTable_Pretty.csv")


#Out of curiosity, hypothetically chronic stress should be something that all of these psychiatric disorders have in common, so meta-analyzing their results should hypothetically improve signal-to-noise for chronic stress
#I'm not sure we have the info to easily run a meta-analysis, so I'm going to try a simple average

DEGs_vs_OtherMetaAnalyses_forTable$AveragePsych_Log2FC<-apply(cbind(DEGs_vs_OtherMetaAnalyses_forTable$Gandal_Microarray_SCZ.beta_log2FC, DEGs_vs_OtherMetaAnalyses_forTable$Gandal_Microarray_BD.beta_log2FC, DEGs_vs_OtherMetaAnalyses_forTable$Gandal_Microarray_MDD.beta_log2FC,
DEGs_vs_OtherMetaAnalyses_forTable$Gandal_Microarray_AAD.beta_log2FC, 
DEGs_vs_OtherMetaAnalyses_forTable$Gandal_RNASeq_SCZ.logFC, 
DEGs_vs_OtherMetaAnalyses_forTable$Gandal_RNASeq_BD.logFC,
DEGs_vs_OtherMetaAnalyses_forTable$Suicide_OFC_Fisher_log2.FC..averaged.,
DEGs_vs_OtherMetaAnalyses_forTable$Suicide_PFC_Fisher_log2.FC..averaged.), 1, function(y) mean(y, na.rm=TRUE))

hist(DEGs_vs_OtherMetaAnalyses_forTable$AveragePsych_Log2FC)


#I should add it to the original data frame too for ease:
DEGs_vs_OtherMetaAnalyses_HumanPsych$AveragePsych_Log2FC<-DEGs_vs_OtherMetaAnalyses_forTable$AveragePsych_Log2FC



pdf("Scatterplot_AveragePsych_Log2FC_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(AveragePsych_Log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych, xlab="Chronic Stress: Log2FC", ylab="Average Psychiatric: Log2FC")

# points(Gandal_RNASeq_SCZ.logFC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych[DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_SCZ.P.Value<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(AveragePsych_Log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = AveragePsych_Log2FC ~ ChronicStress_Xiong_Log2FC_estimate, 
#      data = DEGs_vs_OtherMetaAnalyses_HumanPsych)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.26425 -0.03333  0.00058  0.03936  0.22640 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          0.01132    0.00830   1.364    0.176    
# ChronicStress_Xiong_Log2FC_estimate  0.20907    0.04487   4.659 9.61e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07039 on 102 degrees of freedom
# (31 observations deleted due to missingness)
# Multiple R-squared:  0.1755,	Adjusted R-squared:  0.1674 
# F-statistic: 21.71 on 1 and 102 DF,  p-value: 9.614e-06

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$AveragePsych_Log2FC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="pearson")

# Pearson's product-moment correlation
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$AveragePsych_Log2FC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# t = 4.6592, df = 102, p-value = 9.614e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2461756 0.5658419
# sample estimates:
#       cor 
# 0.4189034 

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$AveragePsych_Log2FC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="spearman")

# Spearman's rank correlation rho
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$AveragePsych_Log2FC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# S = 166802, p-value = 0.2654
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1102007 

#who are the genes that are driving the relationship?

DEGs_vs_OtherMetaAnalyses_HumanPsych$GeneSymbol[which(DEGs_vs_OtherMetaAnalyses_HumanPsych$AveragePsych_Log2FC<(-0.2))]
#[1] "Arc"   "Dusp1" "Egr1"  "Fos" 
#The most extreme overlap is IEGs

DEGs_vs_OtherMetaAnalyses_HumanPsych$GeneSymbol[which((DEGs_vs_OtherMetaAnalyses_HumanPsych$AveragePsych_Log2FC<(-0.05)) & (DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate<(-0.10)))]
# [1] "Arc"     "Atf4"    "Btg2"    "Ccn1"    "Dusp1"   "Egr1"    "Elovl5" 
# [8] "Fos"     "Gdf11"   "Ier2"    "Siae"    "Slc39a8" "Tshr" 

DEGs_vs_OtherMetaAnalyses_HumanPsych$GeneSymbol[which((DEGs_vs_OtherMetaAnalyses_HumanPsych$AveragePsych_Log2FC>(0.05)) & (DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate>(0.10)))]
#[1] "Psrc1"




#I'd be curious in seeing which genes are showing up the most just across the animal chronic stress studies too
#Although note that we only have info for most of the genes from the ELS meta-analysis and chronic cort study, so those studies will dominate the results.

DEGs_vs_OtherMetaAnalyses_forTable$AverageChronicStress_Log2FC<-apply(cbind(DEGs_vs_OtherMetaAnalyses_forTable$ChronicDailyCorticosterone_Jaszczyk_2025_logFC, DEGs_vs_OtherMetaAnalyses_forTable$CSDS_10Days_Reshetnikov_2022_log2FoldChange, DEGs_vs_OtherMetaAnalyses_forTable$CSDS_30Days_Reshetnikov_2022_log2FoldChange,                                                  DEGs_vs_OtherMetaAnalyses_forTable$ChronicStress_Gururajan2022_log2FoldChange, 
                                                                    DEGs_vs_OtherMetaAnalyses_forTable$ELS_Duan_2025_Log2FC_estimate), 1, function(y) mean(y, na.rm=TRUE))

hist(DEGs_vs_OtherMetaAnalyses_forTable$AverageChronicStress_Log2FC)
#Yeah, one of those genes is much more extreme than the others...
#That's going to skew the stats - definitely need non-parametric insight


#I should add it to the original data frame too for ease:
DEGs_vs_OtherMetaAnalyses_HumanPsych$AverageChronicStress_Log2FC<-DEGs_vs_OtherMetaAnalyses_forTable$AverageChronicStress_Log2FC



pdf("Scatterplot_AverageChronicStress_Log2FC_vs_ChronicStress_Xiong_Log2FC_DEGs.pdf", height=5, width=5)

plot(AverageChronicStress_Log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych, xlab="Chronic Stress: Log2FC", ylab="Average Chronic Stress: Log2FC")

# points(Gandal_RNASeq_SCZ.logFC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych[DEGs_vs_OtherMetaAnalyses_HumanPsych$Gandal_RNASeq_SCZ.P.Value<0.05,], pch=16, col=transparent_grey)

TrendLine<-lm(AverageChronicStress_Log2FC~ChronicStress_Xiong_Log2FC_estimate, data=DEGs_vs_OtherMetaAnalyses_HumanPsych) #Fitting a trendline
abline(TrendLine, col="black", lwd=3)

dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = AverageChronicStress_Log2FC ~ ChronicStress_Xiong_Log2FC_estimate, 
#      data = DEGs_vs_OtherMetaAnalyses_HumanPsych)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.79626 -0.06760  0.00927  0.06920  0.49350 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -0.01470    0.01579  -0.931 0.353660    
# ChronicStress_Xiong_Log2FC_estimate  0.32463    0.08746   3.712 0.000324 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1446 on 111 degrees of freedom
# (22 observations deleted due to missingness)
# Multiple R-squared:  0.1104,	Adjusted R-squared:  0.1024 
# F-statistic: 13.78 on 1 and 111 DF,  p-value: 0.0003235

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$AverageChronicStress_Log2FC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="pearson")

# Pearson's product-moment correlation
# 
# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$AverageChronicStress_Log2FC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# t = 3.7117, df = 111, p-value = 0.0003235
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1572052 0.4871146
# sample estimates:
#       cor 
# 0.3322853 

cor.test(DEGs_vs_OtherMetaAnalyses_HumanPsych$AverageChronicStress_Log2FC, DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate, method="spearman")

#Spearman's rank correlation rho

# data:  DEGs_vs_OtherMetaAnalyses_HumanPsych$AverageChronicStress_Log2FC and DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate
# S = 129521, p-value = 2.703e-07
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.4613702

#Interesting. So unlike the human comparison, the pattern isn't driven just by a handful of genes.

DEGs_vs_OtherMetaAnalyses_HumanPsych$GeneSymbol[which((DEGs_vs_OtherMetaAnalyses_HumanPsych$AverageChronicStress_Log2FC>(0.10)) & (DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate>(0.10)))]
#[1] "Chrd"    "Jsrp1"   "Myo1h"   "Nrros"   "Nsl1"    "Slc24a4" "Xdh" 

DEGs_vs_OtherMetaAnalyses_HumanPsych$GeneSymbol[which((DEGs_vs_OtherMetaAnalyses_HumanPsych$AverageChronicStress_Log2FC<(-0.10)) & (DEGs_vs_OtherMetaAnalyses_HumanPsych$ChronicStress_Xiong_Log2FC_estimate<(-0.10)))]
# [1] "9430085M18Rik" "Aqp4"          "BC034090"      "Bfsp2"         "Ccn1"         
# [6] "Fa2h"          "Fam163a"       "Fam163a"       "Fbxo32"        "Gna12"        
# [11] "Lpar1"         "Mid1ip1"       "Nipal4"        "Olfml1"        "Prr18"        
# [16] "Rnd2"          "S1pr5"         "Slc39a8"       "Sntb1"         "Sox9"         
# [21] "Tek"           "Unc5b"         "Zfp979"

#Interesting. There is actually *no overlap* in this gene list with the list of genes driving the positive correlation with the Gandal human neuropsych results.
#... but quite a bit of overlap with the suicide and PTSD results. 