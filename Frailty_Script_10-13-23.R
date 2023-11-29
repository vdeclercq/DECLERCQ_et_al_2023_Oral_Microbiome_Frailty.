
##METADATA ANALYSIS - Descriptive analysis of covariates - sex, AGE, BMI, ETC.
#CHECK FOR NORMALITY
shapiro.test(Frailty_taxa_metadata$AGE)


#TABLES FOR MEDIAN, MEAN, QUARTILES - CONTINUOUS VARIABLES
library(PMCMR)

quantile(Frailty_taxa_metadata$AGE, prob=c(0.25, 0.50, 0.75))
aggregate(AGE ~ Sex, data = Frailty_taxa_metadata, summary)
wilcox.test(Frailty_taxa_metadata$AGE ~ Frailty_taxa_metadata$Sex)


#TABLES FOR COUNTS, PERCENTAGES - CATEGORICAL VARIABLES
library(gmodels)
library(reporttools)

Smoker_Table <- CrossTable(Frailty_taxa_metadata$smoker , 
                        Frailty_taxa_metadata$Sex, 
                        prop.r=TRUE, chisq = TRUE)


##TAXA ANALYSIS - FRAILTY
#preprep data for analysis.
#create taxa table (remove any extra columns & set row names)
library(readr)
taxa_table <- read_csv()

taxa_table <- as.data.frame(taxa_table[,-c(1:2)])
row.names(taxa_table) <- row.names(Frailty_taxa_metadata)
View(taxa_table)

#flip table so taxa are rows, samples columns
rare_filter_table <- t(taxa_table)
View(rare_filter_table)

#create metadata table with only variables of interest
#(filter metadata to match samples in the rare filtered taxa table).
filtered_metadata <- Frailty_taxa_metadata[colnames(rare_filter_table),]
View(filtered_metadata)


#Corncob differential abundance (DA) analysis
library(corncob)
library(phyloseq)
#create the object containing the metadata and taxa to be tested.
otu_tab <- otu_table(rare_filter_table,taxa_are_rows = TRUE)
sam_data <- sample_data(filtered_metadata)
phylo <- phyloseq(otu_tab, sam_data)

#runs corncob DA analysis unadjusted, returns plot and DA results
results <- differentialTest(formula = ~Frailty,
                            formula_null = ~1,
                            phi.formula = ~Frailty,
                            phi.formula_null = ~Frailty,
                            data = phylo,
                            fdr_cutoff = 0.1,
                            test = "Wald")
plot(results)
results$p_fdr
results$significant_taxa
results$significant_models

results3 <- differentialTest(formula = ~Frailty + Sex + smoker + HEIGHT + WEIGHT + VEG + Medications,
                             formula_null = ~1 + Sex + smoker + HEIGHT + WEIGHT + VEG + Medications,
                             phi.formula = ~Frailty + Sex + smoker + HEIGHT + WEIGHT + VEG + Medications,
                             phi.formula_null = ~Frailty + Sex + smoker + HEIGHT + WEIGHT + VEG + Medications,
                             data = phylo,
                             fdr_cutoff = 0.1,
                             test = "Wald")

#runs corncob DA analysis Adjusted, returns plot and DA results
results3 <- differentialTest(formula = ~Frailty + smoker + HEIGHT + WEIGHT + VEG + Medications,
                             formula_null = ~1 + smoker + HEIGHT + WEIGHT + VEG + Medications,
                             phi.formula = ~Frailty + smoker + HEIGHT + WEIGHT + VEG + Medications,
                             phi.formula_null = ~Frailty + smoker + HEIGHT + WEIGHT + VEG + Medications,
                             data = phylo,
                             fdr_cutoff = 0.1,
                             test = "Wald")

plot(results3)
results3$p_fdr
results3$significant_taxa
results3$significant_models


#ALDEx2 DA analysis
#install.packAGEs("BiocManAGEr")
#BiocManAGEr::install("ALDEx2")
library(ALDEx2)
#create model containing the variables to be tested - Unadjusted
matrixmodel <- model.matrix(~Frailty, filtered_metadata) 
View(matrixmodel)
#generate Monte Carlo samples of the Dirichlet distribution, performs centred log-ratio transformation.
CLR <- aldex.clr(rare_filter_table, matrixmodel, mc.samples = 128)
#calculates the expected values for each coefficient of a glm model on the data returned by aldex.clr
GLM <- aldex.glm(CLR, matrixmodel)
View(GLM)

#inclusion of covariates
matrixmodel3 <- model.matrix(~Frailty + Sex + smoker + HEIGHT + WEIGHT + VEG + Medications, filtered_metadata) 
CLR3 <- aldex.clr(rare_filter_table, matrixmodel3, mc.samples = 128)
GLM3 <- aldex.glm(CLR3, matrixmodel3)
View(GLM3)


#Maaslin2 DA analysis
#BiocManAGEr::install("Maaslin2")
library(Maaslin2)

#finds multivariable associations between taxa and metadata based on GLM 
#unadjusted Maaslin2 model
results <- Maaslin2(rare_filter_table, filtered_metadata, "test", 
                    transform = "AST", fixed_effects = c("Frailty" ))
View(results[[1]])

#adjusted model
results3 <- Maaslin2(rare_filter_table, filtered_metadata, "test", 
                     transform = "AST", fixed_effects = c("Sex", "smoker", "HEIGHT", "WEIGHT",
                                                          "VEG", "Medications", "Frailty" ))
View(results3[[1]])


#ANCOM2 DA analysis
library(tidyr)
library(nlme)
library(dplyr)
library(ggplot2)
library(compositions)
library(exactRankTests)
deps = c("exactRankTests", "nlme", "dplyr", "ggplot2", "compositions")

for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(dep)
  }
  library(dep, character.only = TRUE)
}
source("~/Desktop/ANCOM2.R")

#preprocessing step 
preprocess <- feature_table_pre_process(feature_table = rare_filter_table, 
                      meta_data = filtered_metadata, sample_var = "OTU_ID", 
                      out_cut = 0.05,zero_cut = 0.9, lib_cut = 1000)

#run main ANCOM function with preprocessed data, adjusting p-values for multiple comparisons
# without covariates, use NULL
rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
             preprocess$structure_zeros,"Frailty", "BH", 0.1, NULL)
View(rez[[1]])

#with covariates
rez3 <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
              preprocess$structure_zeros,"Frailty", "BH", 0.1, "Sex + smoker + HEIGHT + WEIGHT + VEG + Medications")
View(rez3[[1]])


# Relative Abundance
RA <- sweep(rare_filter_table, 2,colSums(rare_filter_table), "/")
colSums(RA)
flipRA <-data.frame(t(RA))
colnames(flipRA)

identical(rownames(flipRA), rownames(filtered_metadata))

#Scatter plot - RA against Frailty/Age as continuous variables
library(ggplot2)
library(ggpubr) 

DF <- data.frame(flipRA[ , 89],filtered_metadata$Frailty)
ggplot(DF, aes(x=filtered_metadata.Frailty, y=flipRA[ , 89])) +
  geom_point() +   theme_bw() +
  geom_smooth(method=lm, level=0.95, col='red', size=3) + 
  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20)) + 
  ylim(-0.0, 0.00002) +
  labs(title="Lachnospiraceae_.G.2", y= "Relative Abundance", x = "Frailty Index") +  
  theme(plot.title = element_text(hjust = 0.5, size = 24)) 
  
cor.test(filtered_metadata$Frailty, flipRA[ , 10],
         method = "pearson")
cor.test(filtered_metadata$Frailty, flipRA[ , 40],
         method = "spearman", exact=FALSE)


DF <- data.frame(flipRA[ , 19],filtered_metadata$AGE)
ggplot(DF, aes(x=filtered_metadata.AGE, y=flipRA[ , 19])) +
  geom_point() +   theme_bw() +
  geom_smooth(method=lm, level=0.95, col='red', size=3) + 
  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20)) + 
  ylim(-0.001, 0.05) +
  labs(title="Capnocytophaga", y= "Relative Abundance", x = "Age (years)") +  
  theme(plot.title = element_text(hjust = 0.5, size = 24))

cor.test(filtered_metadata$AGE, flipRA[ , 24],
         method = "pearson")
cor.test(filtered_metadata$AGE, flipRA[ , 19],
         method = "spearman", exact=FALSE)



# Beta Diversity
#use of previously generated Bray Curtis dissimilarity matrix or Weighted UniFrac distance matrix.
# Bray Curtis
bray_curtis_distance <-read.table("FILE_NAME",
                      sep="\t", header=TRUE, row.names = 1, check.names = FALSE)

#match metadata samples to bray-curtis samples 
intersect(rownames(Frailty_taxa_metadata), 
          rownames(bray_curtis_distance))
samples_to_keep <- intersect(rownames(Frailty_taxa_metadata),
                             rownames(bray_curtis_distance))
filtered_bray_curtis_distance <- bray_curtis_distance[samples_to_keep, samples_to_keep]
dim(filtered_bray_curtis_distance)
filtered_metadata <- Frailty_taxa_metadata[samples_to_keep,]
rownames(filtered_metadata) <- filtered_metadata[,1]
View(filtered_metadata)

#adonis test - permutational ANOVA of dissimilarities
library(vegan)

set.seed(23)
adonis2(filtered_bray_curtis_distance ~ filtered_metadata$Frailty,
        permutations = 10000, by="margin")

set.seed(23)
adonis2(filtered_bray_curtis_distance ~ filtered_metadata$Frailty
        + filtered_metadata$Sex 
        + filtered_metadata$smoker 
        + filtered_metadata$HEIGHT 
        + filtered_metadata$WEIGHT 
        + filtered_metadata$VEG 
        + filtered_metadata$Medications,
        permutations = 10000, by="margin")


#plotting Bray Curtis dissimilarity
library(vegan)
library(devtools)
library(ggord)
library(ggplot2)
bray_curtis_pcoa <- cmdscale(filtered_bray_curtis_distance, k=2, eig = TRUE)
barplot(bray_curtis_pcoa$eig[1:10])
component1 <- bray_curtis_pcoa$eig[1]/sum(bray_curtis_pcoa$eig)
component2 <- bray_curtis_pcoa$eig[2]/sum(bray_curtis_pcoa$eig)
component1*100
component2*100

plot_data <- data.frame(pc1=bray_curtis_pcoa$points[ ,1],
                        pc2=bray_curtis_pcoa$points[ ,2],
                        Frailty=filtered_metadata$Frailty)
bray_plot <- ggplot(plot_data, aes(x=pc1, y=pc2, colour=Frailty))+ 
  geom_point()+ 
  scale_color_gradient(low = "#1AFF1A" , high= "#4B0092")+
  theme_bw() + xlab("PC1(28.1%)") + ylab("PC2(7.8%)") + theme(legend.position = c(0.09,0.82)) +
  theme(text = element_text(size = 26)) + ylim(-0.35, 0.5) +
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
bray_plot


# Weighted Unifrac
weighted_unifrac <-read.table("FILE_NAME",
                      sep="\t", header=TRUE, row.names = 1, check.names = FALSE)
intersect(rownames(Frailty_taxa_metadata), rownames(weighted_unifrac))
samples_to_keep <- intersect(rownames(Frailty_taxa_metadata), 
                             rownames(weighted_unifrac))
filtered_weighted_unifrac <- weighted_unifrac[samples_to_keep, samples_to_keep]
dim(filtered_weighted_unifrac)
unifrac__metadata <- Frailty_taxa_metadata[samples_to_keep,]

set.seed(23)
adonis2(filtered_weighted_unifrac ~ unifrac__metadata$Frailty, 
        permutations = 10000, by="margin")

set.seed(23)
adonis2(filtered_weighted_unifrac ~ unifrac__metadata$Frailty
        + filtered_metadata$Sex 
        + filtered_metadata$smoker 
        + filtered_metadata$HEIGHT 
        + filtered_metadata$WEIGHT 
        + filtered_metadata$VEG 
        + filtered_metadata$Medications,
        permutations = 10000, by="margin")


#plotting weighted unifrac
weighted_unifrac_pcoa <- cmdscale(filtered_weighted_unifrac, k=2, eig = TRUE)
barplot(weighted_unifrac_pcoa$eig[1:10])
component1 <- weighted_unifrac_pcoa$eig[1]/sum(weighted_unifrac_pcoa$eig)
component2 <- weighted_unifrac_pcoa$eig[2]/sum(weighted_unifrac_pcoa$eig)
component1*100
component2*100

plot_data <- data.frame(pc1=weighted_unifrac_pcoa$points[ ,1],
                        pc2=weighted_unifrac_pcoa$points[ ,2],
                        Frailty=filtered_metadata$Frailty)
unifrac_plot <- ggplot(plot_data, aes(x=pc1, y=pc2, colour=Frailty))+ 
  geom_point()+ 
  scale_color_gradient(low = "#1AFF1A" , high= "#4B0092")+
  theme_bw() + xlab("PC1(52.0%)") + ylab("PC2(19.7%)") + theme(legend.position = c(0.082,0.82)) +
  theme(text = element_text(size = 26)) + ylim(-0.4, 0.35)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
unifrac_plot



##FI components associated with beta diversity
set.seed(23)
adonis2(filtered_weighted_unifrac ~ unifrac__metadata$FI_BMR_kcal, 
        permutations = 10000, by="margin")



##Alph Diversity
library(ggplot2)
library(PMCMR)

#Shannon Diversity Index - taxa diversity in a community
shannon <-read.table("FILE_NAME", sep="\t", 
                     header=TRUE, row.names = 1, check.names = FALSE)
View(shannon)
intersect(rownames(Frailty_taxa_metadata), rownames(shannon))
samples_to_keep1 <- intersect(rownames(Frailty_taxa_metadata), rownames(shannon))
filtered_shannon <- data.frame(shannon[samples_to_keep1,])
filter_metadata <- Frailty_taxa_metadata[samples_to_keep1,]
row.names(filter_metadata) <- filter_metadata$OTU_ID

DF <- data.frame(filtered_shannon$shannon.samples_to_keep1...,filter_metadata$Frailty)
ggplot(DF, aes(x=filter_metadata$Frailty, y=filtered_shannon$shannon.samples_to_keep1...)) +
  geom_point() +   theme_bw() +
  geom_smooth(method=lm, level=0.95, col='red', size=3) + 
  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20)) + 
  ylim(0.2, 6) +
  labs(title="Shannon Diversity", y= "Alpha Diversity", x = "Frailty Index") +  
  theme(plot.title = element_text(hjust = 0.5, size = 24))
cor.test(filter_metadata$Frailty, filtered_shannon$shannon.samples_to_keep1...,
         method = "pearson")
cor.test(filter_metadata$Frailty, filtered_shannon$shannon.samples_to_keep1...,
         method = "spearman", exact = FALSE)


#richness - number of different taxa in a community
richness <-read.table("FILE_NAME", sep="\t", 
                      header=TRUE, row.names = 1, check.names = FALSE)
intersect(rownames(Frailty_taxa_metadata), rownames(richness))
samples_to_keep1 <- intersect(rownames(Frailty_taxa_metadata), rownames(richness))
filtered_richness <- data.frame(richness[samples_to_keep1,])
filter_metadata <- Frailty_taxa_metadata[samples_to_keep1,]

DF <- data.frame(filtered_richness$richness.samples_to_keep1...,filter_metadata$Frailty)
ggplot(DF, aes(x=filter_metadata$Frailty, y=filtered_richness$richness.samples_to_keep1...)) +
  geom_point() +   theme_bw() +
  geom_smooth(method=lm, level=0.95, col='red', size=3) + 
  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20)) + 
  ylim(0.2, 170) +
  labs(title="Observed ASVs", y= "Alpha Diversity", x = "Frailty Index") +  
  theme(plot.title = element_text(hjust = 0.5, size = 24))
cor.test(filter_metadata$Frailty, filtered_richness$richness.samples_to_keep1...,
         method = "pearson")
cor.test(filter_metadata$Frailty, filtered_richness$richness.samples_to_keep1...,
         method = "spearman", exact = FALSE)


#Faith's Phylogenetic diversity - uses phylogentic distances to calculate diversity of community
faiths <-read.table("FILE_NAME", sep="\t", 
                    header=TRUE, row.names = 1, check.names = FALSE)
intersect(rownames(Frailty_taxa_metadata), rownames(faiths))
samples_to_keep1 <- intersect(rownames(Frailty_taxa_metadata), rownames(faiths))
filtered_faiths <- data.frame(faiths[samples_to_keep1,])
filter_metadata <- Frailty_taxa_metadata[samples_to_keep1,]

DF <- data.frame(filtered_faiths$faiths.samples_to_keep1...,filter_metadata$Frailty)
ggplot(DF, aes(x=filter_metadata$Frailty, y=filtered_faiths$faiths.samples_to_keep1...)) +
  geom_point() +   theme_bw() +
  geom_smooth(method=lm, level=0.95, col='red', size=3) + 
  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20)) + 
  ylim(1, 16) +
  labs(title="Faith's Phylogenetic diversity", y= "Alpha Diversity", x = "Frailty Index") +  
  theme(plot.title = element_text(hjust = 0.5, size = 24))
cor.test(filter_metadata$Frailty, filtered_faiths$faiths.samples_to_keep1...,
         method = "pearson")
cor.test(filter_metadata$Frailty, filtered_faiths$faiths.samples_to_keep1...,
         method = "spearman", exact = FALSE)


#Simpson_Evenness
Simpsons_e <-read.table("FILE_NAME", sep="\t", 
                      header=TRUE, row.names = 1, check.names = FALSE)
intersect(rownames(Frailty_taxa_metadata), rownames(Simpsons_e))
samples_to_keep1 <- intersect(rownames(Frailty_taxa_metadata), rownames(Simpsons_e))
filtered_Simpsons_e <- data.frame(Simpsons_e[samples_to_keep1,])
filter_metadata <- Frailty_taxa_metadata[samples_to_keep1,]

DF <- data.frame(filtered_Simpsons_e$Simpsons_e.samples_to_keep1...,filter_metadata$Frailty)
ggplot(DF, aes(x=filter_metadata$Frailty, y=filtered_Simpsons_e$Simpsons_e.samples_to_keep1...)) +
  geom_point() +   theme_bw() +
  geom_smooth(method=lm, level=0.95, col='red', size=3) + 
  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20)) + 
  ylim(0, 0.35) +
  labs(title="Simpson's Evenness", y= "Alpha Diversity", x = "Frailty Index") +  
  theme(plot.title = element_text(hjust = 0.5, size = 24))
cor.test(filter_metadata$Frailty, filtered_Simpsons_e$Simpsons_e.samples_to_keep1...,
         method = "pearson")
cor.test(filter_metadata$Frailty, filtered_Simpsons_e$Simpsons_e.samples_to_keep1...,
         method = "spearman", exact = FALSE)
