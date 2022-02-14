# Shiekhah AL Binali, Kaimrc training problem solution.
# install required packages like MAFtools, Biobase,g3viz,..etc
# 1- Download the dataset linked to above and load the MAF files 
library(maftools)
library(Biobase)
library(BiocGenerics)
library(BiocManager)
library(BiocParallel)
library(BiocVersion)
library(g3viz)
library(dplyr)
library(MAFDash)
# Read MAF files in Rstudio
patient0.maf<- read.maf("mafs/Patient-0.somatic.snvs.maf")
# code reference: https://rdrr.io/bioc/maftools/man/read.maf.html
#view the data
patient0.maf@data
# Apply the pevious two code on all MAF files.
patient1.maf<- read.maf("mafs/Patient-1.somatic.snvs.maf")
patient2.maf<- read.maf("mafs/Patient-2.somatic.snvs.maf")
patient3.maf<- read.maf("mafs/Patient-3.somatic.snvs.maf")
patient4.maf<- read.maf("mafs/Patient-4.somatic.snvs.maf")
patient5.maf<- read.maf("mafs/Patient-5.somatic.snvs.maf")
patient6.maf<- read.maf("mafs/Patient-6.somatic.snvs.maf")
patient7.maf<- read.maf("mafs/Patient-7.somatic.snvs.maf")
patient8.maf<- read.maf("mafs/Patient-8.somatic.snvs.maf")
patient9.maf<- read.maf("mafs/Patient-9.somatic.snvs.maf")
patient10.maf<- read.maf("mafs/Patient-10.somatic.snvs.maf")
patient11.maf<- read.maf("mafs/Patient-11.somatic.snvs.maf")
patient12.maf<- read.maf("mafs/Patient-12.somatic.snvs.maf")
patient13.maf<- read.maf("mafs/Patient-13.somatic.snvs.maf")
patient14.maf<- read.maf("mafs/Patient-14.somatic.snvs.maf")
patient15.maf<- read.maf("mafs/Patient-15.somatic.snvs.maf")
patient16.maf<- read.maf("mafs/Patient-16.somatic.snvs.maf")
patient17.maf<- read.maf("mafs/Patient-17.somatic.snvs.maf")
patient18.maf<- read.maf("mafs/Patient-18.somatic.snvs.maf")
patient19.maf<- read.maf("mafs/Patient-19.somatic.snvs.maf")
patient20.maf<- read.maf("mafs/Patient-20.somatic.snvs.maf")
patient21.maf<- read.maf("mafs/Patient-21.somatic.snvs.maf")
patient22.maf<- read.maf("mafs/Patient-22.somatic.snvs.maf")
patient23.maf<- read.maf("mafs/Patient-23.somatic.snvs.maf")
patient24.maf<- read.maf("mafs/Patient-24.somatic.snvs.maf")
patient25.maf<- read.maf("mafs/Patient-25.somatic.snvs.maf")
patient26.maf<- read.maf("mafs/Patient-26.somatic.snvs.maf")
patient27.maf<- read.maf("mafs/Patient-27.somatic.snvs.maf")
patient28.maf<- read.maf("mafs/Patient-28.somatic.snvs.maf")
patient29.maf<- read.maf("mafs/Patient-29.somatic.snvs.maf")
patient30.maf<- read.maf("mafs/Patient-30.somatic.snvs.maf")
patient31.maf<- read.maf("mafs/Patient-31.somatic.snvs.maf")
patient32.maf<- read.maf("mafs/Patient-32.somatic.snvs.maf")
patient33.maf<- read.maf("mafs/Patient-33.somatic.snvs.maf")
patient34.maf<- read.maf("mafs/Patient-34.somatic.snvs.maf")
patient35.maf<- read.maf("mafs/Patient-35.somatic.snvs.maf")
patient36.maf<- read.maf("mafs/Patient-36.somatic.snvs.maf")
patient37.maf<- read.maf("mafs/Patient-37.somatic.snvs.maf")
patient38.maf<- read.maf("mafs/Patient-38.somatic.snvs.maf")
patient39.maf<- read.maf("mafs/Patient-39.somatic.snvs.maf")
patient40.maf<- read.maf("mafs/Patient-40.somatic.snvs.maf")
patient41.maf<- read.maf("mafs/Patient-41.somatic.snvs.maf")
patient42.maf<- read.maf("mafs/Patient-42.somatic.snvs.maf")
patient43.maf<- read.maf("mafs/Patient-43.somatic.snvs.maf")
patient44.maf<- read.maf("mafs/Patient-44.somatic.snvs.maf")
patient45.maf<- read.maf("mafs/Patient-45.somatic.snvs.maf")
patient46.maf<- read.maf("mafs/Patient-46.somatic.snvs.maf")
patient47.maf<- read.maf("mafs/Patient-47.somatic.snvs.maf")
patient48.maf<- read.maf("mafs/Patient-48.somatic.snvs.maf")
patient49.maf<- read.maf("mafs/Patient-49.somatic.snvs.maf")


# 2- Subset for mutations that are not of the Variant Classification "Silent".
sub_patient0 <- subsetMaf(maf = patient0.maf,query = "Variant_Classification !='Silent'")
View(sub_patient0)
#https://rdrr.io/bioc/maftools/man/subsetMaf.html
# Subset for all patients.

sub_patient1 <- subsetMaf(maf = patient1.maf,query = "Variant_Classification !='Silent'")
sub_patient2 <- subsetMaf(maf = patient2.maf,query = "Variant_Classification !='Silent'")
View(sub_patient2)
sub_patient3 <- subsetMaf(maf = patient3.maf,query = "Variant_Classification !='Silent'")
sub_patient4 <- subsetMaf(maf = patient4.maf,query = "Variant_Classification !='Silent'")
sub_patient5 <- subsetMaf(maf = patient5.maf,query = "Variant_Classification !='Silent'")
sub_patient6 <- subsetMaf(maf = patient6.maf,query = "Variant_Classification !='Silent'")
sub_patient7 <- subsetMaf(maf = patient7.maf,query = "Variant_Classification !='Silent'")
sub_patient8 <- subsetMaf(maf = patient8.maf,query = "Variant_Classification !='Silent'")
sub_patient9 <- subsetMaf(maf = patient9.maf,query = "Variant_Classification !='Silent'")
sub_patient10 <- subsetMaf(maf = patient10.maf,query = "Variant_Classification !='Silent'")
sub_patient11 <- subsetMaf(maf = patient11.maf,query = "Variant_Classification !='Silent'")
sub_patient12 <- subsetMaf(maf = patient12.maf,query = "Variant_Classification !='Silent'")
sub_patient13 <- subsetMaf(maf = patient13.maf,query = "Variant_Classification !='Silent'")
sub_patient14 <- subsetMaf(maf = patient14.maf,query = "Variant_Classification !='Silent'")
sub_patient15 <- subsetMaf(maf = patient15.maf,query = "Variant_Classification !='Silent'")
sub_patient16 <- subsetMaf(maf = patient16.maf,query = "Variant_Classification !='Silent'")
sub_patient17 <- subsetMaf(maf = patient17.maf,query = "Variant_Classification !='Silent'")
sub_patient18 <- subsetMaf(maf = patient18.maf,query = "Variant_Classification !='Silent'")
sub_patient19 <- subsetMaf(maf = patient19.maf,query = "Variant_Classification !='Silent'")
sub_patient20 <- subsetMaf(maf = patient20.maf,query = "Variant_Classification !='Silent'")
sub_patient21 <- subsetMaf(maf = patient21.maf,query = "Variant_Classification !='Silent'")
sub_patient22 <- subsetMaf(maf = patient22.maf,query = "Variant_Classification !='Silent'")
sub_patient23 <- subsetMaf(maf = patient23.maf,query = "Variant_Classification !='Silent'")
sub_patient24 <- subsetMaf(maf = patient24.maf,query = "Variant_Classification !='Silent'")
sub_patient25 <- subsetMaf(maf = patient25.maf,query = "Variant_Classification !='Silent'")
sub_patient26 <- subsetMaf(maf = patient26.maf,query = "Variant_Classification !='Silent'")
sub_patient27 <- subsetMaf(maf = patient27.maf,query = "Variant_Classification !='Silent'")
sub_patient28 <- subsetMaf(maf = patient28.maf,query = "Variant_Classification !='Silent'")
sub_patient29 <- subsetMaf(maf = patient29.maf,query = "Variant_Classification !='Silent'")
sub_patient30 <- subsetMaf(maf = patient30.maf,query = "Variant_Classification !='Silent'")
sub_patient31 <- subsetMaf(maf = patient31.maf,query = "Variant_Classification !='Silent'")
sub_patient32 <- subsetMaf(maf = patient32.maf,query = "Variant_Classification !='Silent'")
sub_patient33 <- subsetMaf(maf = patient33.maf,query = "Variant_Classification !='Silent'")
sub_patient34 <- subsetMaf(maf = patient34.maf,query = "Variant_Classification !='Silent'")
sub_patient35 <- subsetMaf(maf = patient35.maf,query = "Variant_Classification !='Silent'")
sub_patient36 <- subsetMaf(maf = patient36.maf,query = "Variant_Classification !='Silent'")
sub_patient37 <- subsetMaf(maf = patient37.maf,query = "Variant_Classification !='Silent'")
sub_patient38 <- subsetMaf(maf = patient38.maf,query = "Variant_Classification !='Silent'")
sub_patient39 <- subsetMaf(maf = patient39.maf,query = "Variant_Classification !='Silent'")
sub_patient40 <- subsetMaf(maf = patient40.maf,query = "Variant_Classification !='Silent'")
sub_patient41 <- subsetMaf(maf = patient41.maf,query = "Variant_Classification !='Silent'")
sub_patient42 <- subsetMaf(maf = patient42.maf,query = "Variant_Classification !='Silent'")
sub_patient43 <- subsetMaf(maf = patient43.maf,query = "Variant_Classification !='Silent'")
sub_patient44 <- subsetMaf(maf = patient44.maf,query = "Variant_Classification !='Silent'")
sub_patient45 <- subsetMaf(maf = patient45.maf,query = "Variant_Classification !='Silent'")
sub_patient46 <- subsetMaf(maf = patient46.maf,query = "Variant_Classification !='Silent'")
sub_patient47 <- subsetMaf(maf = patient47.maf,query = "Variant_Classification !='Silent'")
sub_patient48 <- subsetMaf(maf = patient48.maf,query = "Variant_Classification !='Silent'")
sub_patient49 <- subsetMaf(maf = patient49.maf,query = "Variant_Classification !='Silent'")
# 3- Find the 15 most common mutations.
# I couldn't find a hint about whether the common mutation for each patient or all of them so I will do both to cover all possibilities.
# First for each file. 
#plot summary for the patient file to observe the information. 
plotmafSummary(maf = sub_patient0, rmOutlier = TRUE, dashboard = TRUE)
# plot the top 15 frequent mutations. 
oncoplot(sub_patient0, top = 15, removeNonMutated = TRUE)
# https://bioconductor.statistik.tu-dortmund.de/packages/3.5/bioc/vignettes/maftools/inst/doc/maftools.html
plotmafSummary(maf = sub_patient1, rmOutlier = TRUE, dashboard = TRUE)
oncoplot(sub_patient1, top = 15, removeNonMutated = TRUE)

plotmafSummary(maf = sub_patient33, rmOutlier = TRUE, dashboard = TRUE)
oncoplot(sub_patient1, top = 15, removeNonMutated = TRUE)

plotmafSummary(maf = sub_patient45, rmOutlier = TRUE, dashboard = TRUE)
oncoplot(sub_patient1, top = 15, removeNonMutated = TRUE)

# The results for frequent mutations for one patient does not make sense, thus I will focus only on all patients.

# Merging patient files into one big file.
all_patientsMAF <- maftools::merge_mafs(list(sub_patient0,sub_patient1,sub_patient2,sub_patient3,sub_patient4,sub_patient5,sub_patient6,
                                             sub_patient7,sub_patient8,sub_patient9, sub_patient10, sub_patient11,sub_patient12,sub_patient13,
                                             sub_patient14,sub_patient15, sub_patient16, sub_patient17, sub_patient18,sub_patient19,sub_patient20,
                                             sub_patient21,sub_patient22,sub_patient23,sub_patient24,sub_patient25,sub_patient26,sub_patient27,sub_patient28,
                                             sub_patient29,sub_patient30,sub_patient31,sub_patient32,sub_patient33,sub_patient34,sub_patient35,sub_patient36,
                                             sub_patient37,sub_patient38,sub_patient39,sub_patient40,sub_patient41,sub_patient42,sub_patient43,sub_patient44,
                                             sub_patient45,sub_patient46,sub_patient47,sub_patient48,sub_patient49), verbose = TRUE)

# https://rdrr.io/github/PoisonAlien/maftools/man/merge_mafs.html and https://github.com/PoisonAlien/maftools/issues/624
plotmafSummary(maf = all_patientsMAF, rmOutlier = TRUE, dashboard = TRUE)
oncoplot(all_patientsMAF, top = 15, removeNonMutated = TRUE)
# https://rdrr.io/bioc/maftools/man/plotmafSummary.html
# https://rdrr.io/bioc/maftools/man/oncoplot.html
 
# 4- Perform a statistical test to explore if any mutated genes are enriched in patients who either responded or not.
# first transfer text file into MAF in case of needed merging. 
#Group the patient files into response/ non-response group, then apply the statistical analysis.   
sample_info <- read.table(file="sample-information.tsv",  header = TRUE)
# subsetting based on response value.
responded_patients <- subset(sample_info, sample_info$Response== 'Responder', select=c(Patient_ID, Tumor_Sample_Barcode,Response))
nonresponded_patients <- subset(sample_info, sample_info$Response== 'Non-Responder', select=c(Patient_ID, Tumor_Sample_Barcode,Response))
# https://www.statmethods.net/management/subset.html
respondMAF  <- maftools::merge_mafs(list(sub_patient1,sub_patient2,sub_patient4,sub_patient7,sub_patient8,sub_patient11,sub_patient12,sub_patient13,
                                             sub_patient14, sub_patient17,sub_patient19,sub_patient24,sub_patient27,sub_patient29,sub_patient30,sub_patient31,
                                         sub_patient33,sub_patient34,sub_patient36,sub_patient37,sub_patient38,sub_patient40,sub_patient42,sub_patient44,
                                             sub_patient49), verbose = TRUE)

nonrespondMAF  <- maftools::merge_mafs(list(sub_patient0,sub_patient3,sub_patient5,sub_patient6, sub_patient9, sub_patient10,sub_patient15, sub_patient16,
                                            sub_patient18,sub_patient20,sub_patient21,sub_patient22,sub_patient23,sub_patient25,sub_patient26,sub_patient28,
                                            sub_patient32,sub_patient35,sub_patient39,sub_patient41,sub_patient43,sub_patient45,sub_patient46,sub_patient47,sub_patient48), verbose = TRUE)

#Setting the minimum mutation to 5 since it the default and there weren't any specification regarding it.
Comparision <- mafCompare(m1 = respondMAF, m2 = nonrespondMAF, m1Name = 'Responder', m2Name = 'Non-Responder',minMut = 5)
print(Comparision)

# After observing the results, Responder group has more mutations than Non responded group. 

# 5- Create a scatter plot of genes with the number of mutated patients on the x-axis and your results from question 4 on the y-axis.
#(The scatter plot did not work, due to time shortage I put the barplot) 

coBarplot(m1 = respondMAF, m2 = nonrespondMAF, m1Name = "Responder ", m2Name = "Non-Responder")
#https://rdrr.io/bioc/maftools/man/coBarplot.html
# 6- How many samples are wild-type versus mutant with respect to the most significantly enriched gene from Question 4?
## First, I need to create 2 groups mutant and wild type. After checking gene summary, I found out all Mutant sample has Missense_Mutation and that why I specified it.

mutatedsample <- subsetMaf(maf = all_patientsMAF,query = "Variant_Classification =='Missense_Mutation'")
wildsample   <-  subsetMaf(maf = all_patientsMAF,query = "Variant_Classification !='Missense_Mutation'")

# plot both groups using the number of nonsynonymous mutations per megabase for each patient in sample information
# Due to time shortage, I couldn't correct the code properly. 
oncoplot( maf = mutatedsample, additionalFeature= "sample_info$Nonsynonymous_mutations_per_Mb" , showTitle= TRUE)
oncoplot( maf = wildsample, additionalFeature= "sample_info$Nonsynonymous_mutations_per_Mb", showTitle= TRUE)

# http://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/oncoplots.html#01_Including_TransitionTransversions_into_oncoplot
