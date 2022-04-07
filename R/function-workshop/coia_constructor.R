#Load packages
library(dplyr)
library(magrittr)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(phyloseq)
library(tidySummarizedExperiment)

#Load data
library(HMP2Data)


#Define functions


prep <- function(data, scale = c(TRUE, FALSE)){
  if(missing(scale)){scale <- TRUE}
  #Center and scale
  scaled <- scale(data, center = TRUE, scale = scale)
  #Fast trace computation
  tr <- sum(scaled * scaled)/(dim(scaled)[1]-1)
  #Normalize magnitude
  scaled/sqrt(tr)
}

COIA_from_2 <- function(data1, data2, data1_type = c("cov.mat", "corr.mat"),
                        data2_type = c("cov.mat", "corr.mat")){
  #Handle missing types
  if(missing(data1_type)){data1_type <- "corr.mat"}
  if(missing(data2_type)){data2_type <- "corr.mat"}
  #Set args
  s1 <- if(data1_type == "corr.mat") TRUE else FALSE
  s2 <- if(data2_type == "corr.mat") TRUE else FALSE
  #Preprocess data
  data1 <- prep(data1, scale = s1)
  data2 <- prep(data2, scale = s2)
  data1_pca <- dudi.pca(data1, scannf = FALSE, nf = 61, center = FALSE, 
                              scale = FALSE)
  data2_pca <- dudi.pca(data2, scannf = FALSE, nf = 61, center = FALSE,
                              scale = FALSE)
  coinertia(data1_pca, data2_pca, scannf = FALSE, nf = 2)
}


#Extending COIA to more than two tables?
COIA <- function(data, types, 
                 option = c("inertia", "lambda1", "uniform", "internal")){
  #Checks for:
  #Data is list of datasets as matrices or dataframes
  if(!(class(data) == "list")) stop("Data must be list of matrices or dataframes.")
  #Set types to all corr.mat if not specified
  if(missing(types)){types <- rep("corr.mat", length(data))}
  #If types are set but number of types specified does not equal number of datasets
  if(length(data) != length(types)) stop("Need type specified for each dataset.")
  #Only possible inputs for types should be corr.mat or cov.mat
  if(!all(unique(types) %in% c("corr.mat", "cov.mat"))) stop("corr.mat and cov.mat are the only known types.")
  
  #Create T/F vector to use in `scale` based on desired matrix type
  internal_test_fun <- function(a){if(a == "corr.mat") TRUE else FALSE}
  s <- sapply(types, internal_test_fun, USE.NAMES = FALSE)
  
  #Run PCA on every dataset
  pca_res <- list()
  for(i in 1:length(data)){
    data[[i]] <- prep(data[[i]], scale = s[i])
    pca_res[[i]] <- dudi.pca(data[[i]], scannf = FALSE, nf = 30, 
                             center = FALSE, scale = FALSE)
  }
  
  #Combine PCA results into `ktab` object
  tabby_kat <- ktab.list.dudi(pca_res)
  
  #Run multiple COIA on `ktab` object
  mcoa(tabby_kat, option = option, scannf = FALSE, nf = length(data)+1)
}

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
momspiCyto <- SummarizedExperiment(momspiCyto_mtx,
                                        colData = momspiCyto_samp,
                                        rowData = data.frame(cytokine = 
                                                               rownames(momspiCyto_mtx)))
#Construct MOMS-PI 16S rRNA SummarizedExperiment object
momspi16S <- SummarizedExperiment(momspi16S_mtx,
                                     colData = momspi16S_samp,
                                     rowData = momspi16S_tax)

momspicy <- MultiAssayExperiment(experiments = list(phy16S = momspi16S,
                                                    cyto = momspiCyto))

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

coia <- COIA_from_2(combined_16S_mtx, combined_Cyto_mtx)
