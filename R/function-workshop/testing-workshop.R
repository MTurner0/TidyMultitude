library(phyloseq)
library(magrittr)
library(testthat)

data("momspiCyto_mtx")
data("momspiCyto_samp")
momspiCyto <- SummarizedExperiment(assays = list(cyto_conc = momspiCyto_mtx),
                                   colData = momspiCyto_samp,
                                   rowData = data.frame(cytokine = 
                                                          rownames(momspiCyto_mtx)))

data("momspi16S_mtx")
data("momspi16S_samp")
data("momspi16S_tax")
momspi16S <- SummarizedExperiment(assays = list(counts = momspi16S_mtx),
                                  colData = momspi16S_samp,
                                  rowData = momspi16S_tax)

momspi_data <- MultiAssayExperiment(experiments = list(phy16S = momspi16S,
                                                       cyto = momspiCyto))

momspi16S_phyloseq <- phyloseq(otu_table(momspi16S_mtx, taxa_are_rows = TRUE), 
                      sample_data(momspi16S_samp), tax_table(momspi16S_tax))

#Combine 16S and cytokines data
#order both sets by visit number within a subject
momspi16S_samp <- momspi16S_samp[
  with(momspi16S_samp, order(subject_id, sample_body_site, visit_number)),
] 

momspiCyto_samp <- momspiCyto_samp[
  with(momspiCyto_samp, order(subject_id, sample_body_site, visit_number)),
]

#select data collected at the same visit
combined_samp <- merge(momspi16S_samp, momspiCyto_samp, 
                       by = c("subject_id", "sample_body_site", 
                              "project_name", "study_full_name",
                              "subject_gender", "subject_race",
                              "visit_number"))

#select data from first visit only
combined_samp <- combined_samp[combined_samp$visit_number ==  1,]

table(combined_samp$sample_body_site)#all vaginal samples

#select 16S data for those samples
combined_16S_phyloseq <- subset_samples(momspi16S_phyloseq, file_name %in% combined_samp$file_name.x)

#get rid of otus that are not observed in any sample for this subset
combined_16S_phyloseq %<>%
  taxa_sums() %>%
  is_greater_than(0) %>%
  prune_taxa(combined_16S_phyloseq)
target1_phy16S <- (otu_table(combined_16S_phyloseq))@.Data

target1_cyto <- momspiCyto_mtx[, colnames(momspiCyto_mtx) %in% combined_samp$file_name.y ]

checkpoint1 <- momspi_data %>% 
  # Order both sets by visit number within a subject
  arrange_colData(subject_id, sample_body_site, visit_number) %>% 
  # Use metadata to select data collected at
  # (1) the same visit (across the two experiments)
  intersect_colData(by = c("subject_id", "sample_body_site",
                           "project_name", "study_full_name", "subject_gender",
                           "subject_race", "visit_number")) %>% 
  # and (2) the first visit
  filter_colData(visit_number == 1) %>% 
  # Remove OTUs not observed in any subject
  trim_empty_rows(phy16S) %>% 
  assays()

checkpoint1_phy16S <- checkpoint1[["phy16S"]]
checkpoint1_cyto <- checkpoint1[["cyto"]]

expect_equal(apply(checkpoint1_phy16S, 1, mean), apply(target1_phy16S, 1, mean))

expect_equal(apply(checkpoint1_cyto, 1, mean), apply(target1_cyto, 1, mean))

## PASSES CHECKPOINT 1

combined_16S_mtx <- t(target1_phy16S)
combined_Cyto_mtx <- t(target1_cyto)

combined_16S_mtx <- combined_16S_mtx/apply(combined_16S_mtx, 1, sum)

taxa_mtx <- scale(combined_16S_mtx, center = TRUE, scale = FALSE)
taxa_tr <- sum(taxa_mtx*taxa_mtx)/(dim(taxa_mtx)[1]-1)
taxa_mtx <- taxa_mtx/sqrt(taxa_tr)
target2_phy16S <- t(taxa_mtx)

cyto_mtx <- scale(combined_Cyto_mtx, center = TRUE, scale = TRUE)
cyto_tr <- sum(cyto_mtx*cyto_mtx)/(dim(cyto_mtx)[1]-1)
cyto_mtx <- cyto_mtx/sqrt(cyto_tr)
target2_cyto <- t(cyto_mtx)

# Define my preprocessing functions
prep <- function(data, scale = c(TRUE, FALSE)){
  # Note that the data is transposed so that these steps are performed rowwise
  scaled <- scale(t(data), center = TRUE, scale = scale) # Center and scale
  tr <- sum(scaled * scaled)/(dim(scaled)[1]-1) # Fast trace computation
  processed <- scaled/sqrt(tr) # Normalize magnitude
  return(t(processed)) # Rotate back
}

expect_equal(prep(target1_cyto, scale = TRUE), target2_cyto)

(t(target1_phy16S)/colSums(target1_phy16S)) %>% 
  t() %>% 
  prep(., scale = FALSE) %>% 
  expect_equal(., target2_phy16S)

checkpoint2 <- momspi_data %>% 
  # Order both sets by visit number within a subject
  arrange_colData(subject_id, sample_body_site, visit_number) %>% 
  # Use metadata to select data collected at
  # (1) the same visit (across the two experiments)
  intersect_colData(by = c("subject_id", "sample_body_site",
                           "project_name", "study_full_name", "subject_gender",
                           "subject_race", "visit_number")) %>% 
  # and (2) the first visit
  filter_colData(visit_number == 1) %>% 
  # Remove OTUs not observed in any subject
  trim_empty_rows(phy16S) %>% 
  # For phy16S data: Convert taxa to proportions
  # Then, by row -- center (but do not scale) and normalize magnitude
  transmute(phy16S,
            cov = (t(counts)/colSums(counts))%>% 
              t() %>% 
              prep(., scale = FALSE)) %>% 
  # For cytokines data: by row -- center, 
  transmute(cyto,
            corr = prep(cyto_conc, scale = TRUE)) %>% 
  assays()

checkpoint2_phy16S <- checkpoint2[["phy16S"]]
checkpoint2_cyto <- checkpoint2[["cyto"]]

expect_equal(apply(checkpoint2_phy16S, 1, mean), apply(target2_phy16S, 1, mean))
expect_equal(apply(checkpoint2_cyto, 1, mean), apply(target2_cyto, 1, mean))

## PASSES CHECKPOINT 2

target3_phy16S <- ade4::dudi.pca(taxa_mtx, scannf=FALSE, nf =61,
                           center = FALSE, scale = FALSE)
target3_cyto <- ade4::dudi.pca(cyto_mtx, scannf=FALSE, nf =61,
                           center = FALSE, scale = FALSE)

checkpoint3 <- checkpoint2 %>% 
  # Ensure that samples are in rows, variables are in columns
  lapply(., t) %>% 
  # Perform PCA on both experiments
  lapply(., function(x) {
    ade4::dudi.pca(x, scannf = FALSE, nf = 61, 
                   center = FALSE, scale = FALSE)
  })

checkpoint3_phy16S <- checkpoint3[["phy16S"]]
checkpoint3_cyto <- checkpoint3[["cyto"]]

expect_equal(checkpoint3_phy16S$eig, target3_phy16S$eig)
expect_equal(checkpoint3_cyto$eig, target3_cyto$eig)

## PASSES CHECKPOINT 3

target_coin <- ade4::coinertia(target3_phy16S, target3_cyto, 
                               scannf = FALSE, nf = 2)

checkpoint_coin <- ade4::coinertia(checkpoint3_phy16S, checkpoint3_cyto,
                                   scannf = FALSE, nf = 2)

expect_equal(target_coin$RV, checkpoint_coin$RV)