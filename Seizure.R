require(R.matlab)
require(plyr)
require(ggplot2)
require(reshape2)
library(tseries)
#-- Data is stored in 

get.seizure.file = function(dog.type,dog.count,seizure.type,file.count){
   #data.dir = '/run/user/1000/gvfs/smb-share:server=researchfiles.iastate.edu,share=las$/STAT/KaggleDataComp/data/raw_data/'
   #data.dir = '/media/yet/DATA/seizure'
  data.dir <- "K:/data/raw_data/"
   data.file.gen = c('_interictal_segment_', '_preictal_segment_','_test_segment_')
   data.folder = paste(dog.type,dog.count,sep='_')
   files = list.files(paste(data.dir,data.folder,sep=''))
   files_interictal = files[grepl('interictal',files)]
   files_preictal = files[grepl('preictal',files)]
   files_test = files[grepl('test',files)]

   if(seizure.type == 'preictal'){get.files = files_preictal}
   if(seizure.type == 'interictal'){get.files = files_interictal}
   if(seizure.type == 'test'){get.files = files_test}

   if(dog.type == 'Dog'){
      data.file = paste(data.dir,data.folder,get.files[file.count],sep='/')
   }
   return(readMat(data.file))
}


# Dog2

#note <- read.table("/run/user/1000/gvfs/smb-share:server=researchfiles.iastate.edu,share=las$/STAT/KaggleDataComp/1 - important notes/workflow.txt")

# Dog2: preictal 42 interictal 500 test 1000
# get output of pacf and acf for 42 preictals  segments of Dog 2
# Dog2 has 16 channels

# preictal ############
n_segment <- c(42, 500, 1000)
seizure <- c("preictal", "interictal", "test")
type <- c("Dog", "Patient")
i_type <- c(1:2)
i_dog <- c(1,2,3,4,5)
i_patient <- c(1,2)
dog2_out <- data.frame(Date=as.Date(character()),
                 File=character(), 
                 User=character(), 
                 stringsAsFactors=FALSE) 
pm1 <- proc.time()
for(i_seizure in 1:3){ # i_seizure <- 1
for(j in 1:n_segment[i_seizure]){
  d <- get.seizure.file(type[1],i_dog[2],seizure[i_seizure],j) # j <- 1
  dat_mat <- d[[1]][[1]]
  n_channel <- nrow(dat_mat)
  res <- ldply(1:n_channel, function(i){ # i <- 1
    p <- 10
    des_out <- c(quantile(dat_mat[i,], probs = c(0, .25, .5, .75, 1)),
                 sd(dat_mat[i,]),
                 mean(dat_mat[i,]),
                 pacf(dat_mat[i,], plot = F)$acf[1:p])
    
    ar_fit <- arima(dat_mat[i,], order = c(p,0,0))
    ar_out <- c(ar_fit$coef, ar_fit$sigma2, ar_fit$aic)
    out <- c(dog = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
             des_out, ar_out)
    names(out)  <- c("dog", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                     "75_quantile", "100_quantile", 
                     "mean", "sd", paste("pacf_", 1:p, sep = ""), 
                     paste("ar", p,"_ar", 1:p, sep = ""), "intercept",  
                     "sigma2", "aic")
    print(paste(type[1],i_dog[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
    return(out)
    
  })
dog2_out <- rbind(dog2_out, res)
}
}

proc.time()-pm1
write.csv(dog2_out, file = "dog2_ar10.csv", row.names = FALSE)
