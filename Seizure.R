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
J <- 42
dog2_preictal_out <- data.frame(Date=as.Date(character()),
                 File=character(), 
                 User=character(), 
                 stringsAsFactors=FALSE) 
pm1 <- proc.time()
for(j in 1:J){
  d <- get.seizure.file('Dog',2,'preictal',j)
  dat_mat <- d[[1]][[1]]
  n_channel <- nrow(dat_mat)
  res <- ldply(1:n_channel, function(i){
    quantile_out <- quantile(dat_mat[i,], probs = c(0, .25, .5, .75, 1))
    sd_out <- sd(dat_mat[i,])
    mean_out <- mean(dat_mat[i,])
    pacf_out <- pacf(dat_mat[i,], plot = F)$acf[1:3]
    
    ar_out <- arma(dat_mat[i,], order = c(3,0))
    summary_ar_out <- summary(ar_out)
    coef_ar_out <-summary_ar_out$coef[,1]  # coef of AR(3)
    sigma2_ar_out <- summary_ar_out$var
    aic_ar_out <- summary_ar_out$aic
    
    
    out <- c(dog = 2, seizure = "preictal", segment = j, channel = i, quantile_out, mean_out, sd_out, pacf_out, 
             coef_ar_out, sigma2_ar_out, aic_ar_out)
    names(out)  <- c("dog", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                     "75_quantile", "100_quantile", 
                     "mean", "sd", "pacf_1", "pacf_2", "pacf_3",
                     "ar3_ar1", "ar3_ar2", "ar3_ar3", "intercept",
                     "sigma2", "aic")
    return(out)
    
  })
dog2_preictal_out <- rbind(dog2_preictal_out, res)

}
proc.time()-pm1
write.csv(dog2_preictal_out, file = "dog2_preictal_ar3.csv", row.names = FALSE)



# interictal segment#############
J <- 2
dog2_interictal_out <- data.frame(Date=as.Date(character()),
                                File=character(), 
                                User=character(), 
                                stringsAsFactors=FALSE) 
pm1 <- proc.time()
for(j in 1:J){
  d <- get.seizure.file('Dog',2,'interictal',j)
  dat_mat <- d[[1]][[1]]
  n_channel <- nrow(dat_mat)
  res <- ldply(1:n_channel, function(i){
    quantile_out <- quantile(dat_mat[i,], probs = c(0, .25, .5, .75, 1))
    sd_out <- sd(dat_mat[i,])
    mean_out <- mean(dat_mat[i,])
    pacf_out <- pacf(dat_mat[i,], plot = F)$acf[1:3]
    
    ar_out <- arma(dat_mat[i,], order = c(3,0))
    summary_ar_out <- summary(ar_out)
    coef_ar_out <-summary_ar_out$coef[,1]  # coef of AR(3)
    sigma2_ar_out <- summary_ar_out$var
    aic_ar_out <- summary_ar_out$aic
    
    
    out <- c(dog = 2, seizure = "interictal", segment = j, channel = i, quantile_out, mean_out, sd_out, pacf_out, 
             coef_ar_out, sigma2_ar_out, aic_ar_out)
    names(out)  <- c("dog", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                     "75_quantile", "100_quantile", 
                     "mean", "sd", "pacf_1", "pacf_2", "pacf_3",
                     "ar3_ar1", "ar3_ar2", "ar3_ar3", "intercept",
                     "sigma2", "aic")
    return(out)
    
  })
  dog2_interictal_out <- rbind(dog2_interictal_out, res)
  
}
proc.time()-pm1
write.csv(dog2_interictal_out, file = "dog2_interictal_ar3.csv", row.names = FALSE)
