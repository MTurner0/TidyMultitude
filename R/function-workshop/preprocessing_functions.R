#Order columns








# First batch will work with one dataset at a time
# Need more flexible normalization methods than recipes::step_normalize()


prep <- function(data, scale = c(TRUE, FALSE)){
  if(missing(scale)){scale <- TRUE}
  #Center and scale
  scaled <- scale(data, center = TRUE, scale = scale)
  #Fast trace computation
  tr <- sum(scaled * scaled)/(dim(scaled)[1]-1)
  #Normalize magnitude
  scaled/sqrt(tr)
}

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



#Normalize magnitude
normalize <- function(data){
  #Fast trace computation
  tr <- sum(data * data)/(dim(data)[1]-1)
  normed <- data/sqrt(tr)
}