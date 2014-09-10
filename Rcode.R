library(R.matlab)
library(plyr)

path <- "/media/yet/DATA/seizure/Dog_"
Dogid <- 5

seq_data <- NULL
state <- c("_interictal", "_preictal", "_test")

#   str(data)
#   str(data[[1]])
#   str(data[[1]][[1]])# data matrix
#   str(data[[1]][[2]]) # data.length.sec
#   str(data[[1]][[3]])# sampling.frequency
#   str(data[[1]][[4]])# chanels
#   str(data[[1]][[5]])# sequence
#   attr(data[[1]], "dimnames")
#   attr(data, "header")
i <- 1
j <- 1
pm1 <- proc.time()
for(i in 1:50){
  segmentid <- sprintf("%04d", i)
  filename <- paste("Dog_", Dogid, state[j], "_segment_", segmentid, ".mat", sep = "")
  pathname <- file.path(paste(path, Dogid, sep = ""), filename)
  data <- readMat(pathname)
  seq_data[i] <- as.vector(data[[1]][[5]])
    
}
proc.time() -pm1
seq_data


pm1 <- proc.time()
seq_data <- laply(1:50, function(i){
  segmentid <- sprintf("%04d", i)
  filename <- paste("Dog_", Dogid, state[j], "_segment_", segmentid, ".mat", sep = "")
  pathname <- file.path(paste(path, Dogid, sep = ""), filename)
  data <- readMat(pathname)
  return(as.vector(data[[1]][[5]]))
}
)
proc.time() -pm1

### extract matrix data ####

pm1 <- proc.time()
matrix_data <- llply(1:50, function(i){
  segmentid <- sprintf("%04d", i)
  filename <- paste("Dog_", Dogid, state[j], "_segment_", segmentid, ".mat", sep = "")
  pathname <- file.path(paste(path, Dogid, sep = ""), filename)
  data <- readMat(pathname)
  return(summary(data[[1]][[1]]))
}
)
proc.time() -pm1
str(matrix_data)
#gc()
#gcinfo(FALSE)
