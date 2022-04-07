#Load packages
library(dplyr)
library(magrittr)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(phyloseq)

#Load data
library(HMP2Data)

#16S rRNA data
data("momspi16S_mtx") #OTU matrix
data("momspi16S_tax") #Taxonomy matrix
data("momspi16S_samp") #Sample matrix
#Note relationships between matrices
all.equal(colnames(momspi16S_mtx), rownames(momspi16S_samp))
all.equal(rownames(momspi16S_mtx), rownames(momspi16S_tax))

#Cytokines
data("momspiCyto_mtx") #Cytokine matrix
data("momspiCyto_samp") #Sample matrix
#Note relationchip between matrices
all.equal(colnames(momspi16S_mtx), rownames(momspi16S_samp))

#Manually construct MOMS-PI 16S rRNA phyloseq object
momspi16S <- phyloseq(otu_table(momspi16S_mtx, taxa_are_rows = TRUE),
                      tax_table(momspi16S_tax),
                      sample_data(momspi16S_samp))

#Manually construct MOMS-PI Cytokines SummarizedExperiment object
momspiCytokines <- SummarizedExperiment(momspiCyto_mtx,
                                        colData = momspiCyto_samp,
                                        rowData = data.frame(cytokine = 
                                                               rownames(momspiCyto_mtx)))

#Order samples by visit number
momspi16S_samp <- momspi16S_samp[with(momspi16S_samp, order(subject_id, sample_body_site, visit_number))]
momspiCyto_samp <- momspiCyto_samp[with(momspiCyto_samp, order(subject_id, sample_body_site, visit_number))]

#Select data collected at the same visit
combined_samp <- merge(momspi16S_samp, momspiCyto_samp,
                       by = c("subject_id", "sample_body_site",
                              "project_name", "study_full_name",
                              "subject_gender", "subject_race",
                              "visit_number"))

#Select data from the first visit only
combined_samp <- combined_samp[combined_samp$visit_number == 1, ]

#Select 16S data for those samples
combined_16S_phyloseq <- subset_samples(momspi16S,
                                        file_name %in% combined_samp$file_name.x)

#Remove OTUs not observed in any sample for this subset
combined_16S_phyloseq %<>% taxa_sums() %>% 
  is_greater_than(0) %>% 
  prune_taxa(combined_16S_phyloseq)
#Extract OTU table from phyloseq model
combined_16S_mtx <- otu_table(combined_16S_phyloseq)

#Use sample data to subsample cytokine data
combined_Cyto_mtx <- momspiCyto_mtx[, colnames(momspiCyto_mtx) %in% combined_samp$file_name.y]

#Ensure that samples are in rows and variables are in columns
combined_16S_mtx %<>% t
combined_Cyto_mtx %<>% t

#I think it would be difficult to automate all of these steps. In what format can I expect the data to be?

#Convert taxa to proportions: divide each column by its sum
#Could you skip this step and scale in the next? Would a flexible scaling argument in the function be useful?
combined_16S_mtx <- combined_16S_mtx/apply(combined_16S_mtx, 1, sum)
#Scale
taxa_mtx <- scale(combined_16S_mtx, center = TRUE, scale = FALSE)
#Compute trace
taxa_tr <- sum(taxa_mtx*taxa_mtx)/(dim(taxa_mtx)[1]-1)
#Normalize magnitude by dividing each value of X by sqrt(trace)
taxa_mtx <- taxa_mtx/sqrt(taxa_tr)

#Scale and normalize cytokines data
cyto_mtx <- scale(combined_Cyto_mtx, center = TRUE, scale = TRUE)
cyto_tr <- sum(cyto_mtx*cyto_mtx)/(dim(cyto_mtx)[1]-1)
cyto_mtx <- cyto_mtx/sqrt(cyto_tr)

#Get PCA objects
taxa_pca <- ade4::dudi.pca(taxa_mtx, scannf = FALSE, nf = 61,
                           center = FALSE, scale = FALSE)
cyto_pca <- ade4::dudi.pca(cyto_mtx, scannf = FALSE, nf = 61,
                          center = FALSE, scale = FALSE)

#COIA
coin <- ade4::coinertia(taxa_pca, cyto_pca, scannf = FALSE, nf = 2)
#Get RV coefficient
coin$RV
