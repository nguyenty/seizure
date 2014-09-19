require(R.matlab)
require(plyr)
require(ggplot2)
require(reshape2)
#install.packages("FitAR")
library(FitAR)
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

##### seizure = 1#####preictal
dog2_preictal_1_42_ar6_ar7 <- data.frame(Date=as.Date(character()),
                 File=character(), 
                 User=character(), 
                 stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 1){ # i_seizure <- 1
for(j in 1:n_segment[i_seizure]){
  pm2 <- proc.time()
  d <- get.seizure.file(type[1],i_dog[2],seizure[i_seizure],j) # j <- 1
  dat_mat <- d[[1]][[1]]
  n_channel <- nrow(dat_mat)
  res <- ldply(1:n_channel, function(i){ # i <- 1
    
    Dog.stats <- c(quantile(dat_mat[i,], probs = c(0, .25, .5, .75, 1)),
                 sd(dat_mat[i,]),
                 mean(dat_mat[i,]),
                 pacf(dat_mat[i,], plot = F)$acf[1:7])
    
    d.fit7 <- FitAR(dat_mat[i,], 7)
    d.fit6 <- FitAR(dat_mat[i,], 6)  
    AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*6 - 2*d.fit6$loglikelihood)
    AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*7 - 2*d.fit7$loglikelihood)
    
    out <- c(dog = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
             Dog.stats,AIC6.stats,AIC7.stats)
    names(out)  <- c("dog", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                     "75_quantile", "100_quantile", 
                     "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                     paste("ar", 6,"_ar", 1:6, sep = ""),
                     "ar6_sigma2", "ar6_aic", "ar7_mu",
                     paste("ar", 7,"_ar", 1:7, sep = ""),
                     "ar7_sigma2", "ar7_aic" )
    print(paste(type[1],i_dog[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
    return(out)
    
  })
  proc.time()-pm2
  
  dog2_preictal_1_42_ar6_ar7 <- rbind(dog2_preictal_1_42_ar6_ar7, res)
}
}

proc.time()-pm1
write.csv(dog2_preictal_1_42_ar6_ar7, file = "dog2_preictal_1_42_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(dog2_preictal_1_42_ar6_ar7, file = "K:/Yet/output/dog2_preictal_1_42_ar6_ar7.csv", row.names = FALSE)

#### seizure = 2 interictal 1- 250######
dog2_out_interictal_1_250 <- data.frame(Date=as.Date(character()),
                       File=character(), 
                       User=character(), 
                       stringsAsFactors=FALSE) 
pm1 <- proc.time()
for(i_seizure in 2){ # i_seizure <- 1
  for(j in 1:250){
    pm2 <- proc.time()
    d <- get.seizure.file(type[1],i_dog[2],seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    res <- ldply(1:n_channel, function(i){ # i <- 1
      
      Dog.stats <- c(quantile(dat_mat[i,], probs = c(0, .25, .5, .75, 1)),
                     sd(dat_mat[i,]),
                     mean(dat_mat[i,]),
                     pacf(dat_mat[i,], plot = F)$acf[1:7])
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*6 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*7 - 2*d.fit7$loglikelihood)
      
      out <- c(dog = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("dog", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[1],i_dog[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    dog2_out_interictal_1_250  <- rbind(dog2_out_interictal_1_250 , res)
  }
}

proc.time()-pm1
write.csv(dog2_out_interictal_1_250 , file = "dog2_interictal_1_250_ar6_ar7.csv", row.names = FALSE)
write.csv(dog2_out_interictal_1_250 , file = "K:/Yet/output/dog2_interictal_1_250_ar6_ar7.csv", row.names = FALSE)



#### seizure = 2 interictal 251- 500######
dog2_out_interictal_251_500 <- data.frame(Date=as.Date(character()),
                                        File=character(), 
                                        User=character(), 
                                        stringsAsFactors=FALSE) 
pm1 <- proc.time()
for(i_seizure in 2){ # i_seizure <- 1
  for(j in 251:500){
    pm2 <- proc.time()
    d <- get.seizure.file(type[1],i_dog[2],seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    res <- ldply(1:n_channel, function(i){ # i <- 1
      
      Dog.stats <- c(quantile(dat_mat[i,], probs = c(0, .25, .5, .75, 1)),
                     sd(dat_mat[i,]),
                     mean(dat_mat[i,]),
                     pacf(dat_mat[i,], plot = F)$acf[1:7])
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*6 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*7 - 2*d.fit7$loglikelihood)
      
      out <- c(dog = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("dog", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[1],i_dog[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    dog2_out_interictal_251_500  <- rbind(dog2_out_interictal_251_500 , res)
  }
}

proc.time()-pm1
write.csv(dog2_out_interictal_251_500 , file = "dog2_interictal_251_500_ar6_ar7.csv", row.names = FALSE)
write.csv(dog2_out_interictal_251_500 , file = "K:/Yet/output/dog2_interictal_251_500_ar6_ar7.csv", row.names = FALSE)

#### seizure = 3 test 1- 250######
dog2_out_test_1_250 <- data.frame(Date=as.Date(character()),
                                          File=character(), 
                                          User=character(), 
                                          stringsAsFactors=FALSE) 
pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 1:250){
    pm2 <- proc.time()
    d <- get.seizure.file(type[1],i_dog[2],seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    res <- ldply(1:n_channel, function(i){ # i <- 1
      
      Dog.stats <- c(quantile(dat_mat[i,], probs = c(0, .25, .5, .75, 1)),
                     sd(dat_mat[i,]),
                     mean(dat_mat[i,]),
                     pacf(dat_mat[i,], plot = F)$acf[1:7])
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*6 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*7 - 2*d.fit7$loglikelihood)
      
      out <- c(dog = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("dog", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[1],i_dog[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    dog2_out_test_1_250  <- rbind(dog2_out_test_1_250, res)
  }
}

proc.time()-pm1
write.csv(dog2_out_test_1_250, file = "dog2_test_1_250_ar6_ar7.csv", row.names = FALSE)
write.csv(dog2_out_test_1_250, file = "K:/Yet/output/dog2_test_1_250_ar6_ar7.csv", row.names = FALSE)

#### seizure = 3 test 251- 500######
dog2_out_test_251_500 <- data.frame(Date=as.Date(character()),
                                  File=character(), 
                                  User=character(), 
                                  stringsAsFactors=FALSE) 
pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 251:500){
    pm2 <- proc.time()
    d <- get.seizure.file(type[1],i_dog[2],seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    res <- ldply(1:n_channel, function(i){ # i <- 1
      
      Dog.stats <- c(quantile(dat_mat[i,], probs = c(0, .25, .5, .75, 1)),
                     sd(dat_mat[i,]),
                     mean(dat_mat[i,]),
                     pacf(dat_mat[i,], plot = F)$acf[1:7])
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*6 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*7 - 2*d.fit7$loglikelihood)
      
      out <- c(dog = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("dog", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[1],i_dog[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    dog2_out_test_251_500  <- rbind(dog2_out_test_251_500, res)
  }
}

proc.time()-pm1
write.csv(dog2_out_test_251_500, file = "dog2_test_251_500_ar6_ar7.csv", row.names = FALSE)
write.csv(dog2_out_test_251_500, file = "K:/Yet/output/dog2_test_251_500_ar6_ar7.csv", row.names = FALSE)

#### seizure = 3 test 501- 750######
dog2_out_test_501_750 <- data.frame(Date=as.Date(character()),
                                    File=character(), 
                                    User=character(), 
                                    stringsAsFactors=FALSE) 
pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 501:750){
    pm2 <- proc.time()
    d <- get.seizure.file(type[1],i_dog[2],seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    res <- ldply(1:n_channel, function(i){ # i <- 1
      
      Dog.stats <- c(quantile(dat_mat[i,], probs = c(0, .25, .5, .75, 1)),
                     sd(dat_mat[i,]),
                     mean(dat_mat[i,]),
                     pacf(dat_mat[i,], plot = F)$acf[1:7])
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*6 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*7 - 2*d.fit7$loglikelihood)
      
      out <- c(dog = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("dog", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[1],i_dog[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    dog2_out_test_501_750  <- rbind(dog2_out_test_501_750, res)
  }
}

proc.time()-pm1
write.csv(dog2_out_test_501_750, file = "dog2_test_501_750_ar6_ar7.csv", row.names = FALSE)
write.csv(dog2_out_test_501_750, file = "K:/Yet/output/dog2_test_501_750_ar6_ar7.csv", row.names = FALSE)

#### seizure = 3 test 751- 1000######
dog2_out_test_751_1000 <- data.frame(Date=as.Date(character()),
                                    File=character(), 
                                    User=character(), 
                                    stringsAsFactors=FALSE) 
pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 751:1000){
    pm2 <- proc.time()
    d <- get.seizure.file(type[1],i_dog[2],seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    res <- ldply(1:n_channel, function(i){ # i <- 1
      
      Dog.stats <- c(quantile(dat_mat[i,], probs = c(0, .25, .5, .75, 1)),
                     sd(dat_mat[i,]),
                     mean(dat_mat[i,]),
                     pacf(dat_mat[i,], plot = F)$acf[1:7])
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*6 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*7 - 2*d.fit7$loglikelihood)
      
      out <- c(dog = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("dog", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[1],i_dog[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    dog2_out_test_751_1000  <- rbind(dog2_out_test_751_1000, res)
  }
}

proc.time()-pm1
write.csv(dog2_out_test_751_1000, file = "dog2_test_751_1000_ar6_ar7.csv", row.names = FALSE)
write.csv(dog2_out_test_751_1000, file = "K:/Yet/output/dog2_test_751_1000_ar6_ar7.csv", row.names = FALSE)
