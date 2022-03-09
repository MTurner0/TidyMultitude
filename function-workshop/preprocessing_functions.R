# First batch will work with one dataset at a time
# Need more flexible normalization methods than recipes::step_normalize()

#Normalize magnitude
magnorm <- function(data){
  #Fast trace computation
  tr <- sum(data * data)/(dim(data)[1]-1)
  normed <- data/sqrt(tr)
}