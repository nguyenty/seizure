require(R.matlab)
require(plyr)
require(ggplot2)
require(reshape2)
library(reshape)
#install.packages("FitAR")
library(FitAR)
library(randomForest)
#-- Data is stored in 
# dog.type <- "Patient"
# dog.count <- 2
# seizure.type = "preictal"
# file.count <- 1
# data_dir <- "K:/data/raw_data/Patient_2/Patient_2_test_segment_0002.mat"
# dat <-readMat(data_dir)
# str(dat)
# 
# output <- read.table("K:/Yihua/Dog5_test2_AR6_AR7.txt", heade = T)
# str(output)
# summary(output$AR7.aic)
# 

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
   data.file = paste(data.dir,data.folder,get.files[file.count],sep='/')
   
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


n_segment <- c(18, 42, 150)
seizure <- c("preictal", "interictal", "test")
type <- c("Dog", "Patient")
i_type <- c(1:2)
i_dog <- c(1,2,3,4,5)
i_patient <- c(1,2)

#####Patien2_ seizure = 1 preictal 1-18#####
## Patient 2 interictal 42 preictal 18 test 150###
patient2_preictal_1_18_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                         File=character(), 
                                         User=character(), 
                                         stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 1){ # i_seizure <- 1
  for(j in 1:n_segment[i_seizure]){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],i_patient[2],seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    res <- ldply(1:n_channel, function(i){ # i <- 5
      
      Dog.stats <- c(quantile(dat_mat[i,], probs = c(0, .25, .5, .75, 1)),
                     sd(dat_mat[i,]),
                     mean(dat_mat[i,]),
                     pacf(dat_mat[i,], plot = F)$acf[1:7])
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*6 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*7 - 2*d.fit7$loglikelihood)
      
      out <- c(patient = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("patient", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],i_patient[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    patient2_preictal_1_18_ar6_ar7 <- rbind(patient2_preictal_1_18_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(patient2_preictal_1_18_ar6_ar7, file = "patient2_preictal_1_18_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(patient2_preictal_1_18_ar6_ar7, file = "K:/Yet/output/patient2_preictal_1_18_ar6_ar7.csv", row.names = FALSE)

#####Patien2_ seizure = 2 interictal 1- 21#####
## Patient 2 interictal 42 preictal 18 test 150###
patient2_interictal_1_21_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                             File=character(), 
                                             User=character(), 
                                             stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 2){ # i_seizure <- 1
  for(j in 1:21){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],i_patient[2],seizure[i_seizure],j) # j <- 1
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
      
      out <- c(patient = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("patient", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],i_patient[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    patient2_interictal_1_21_ar6_ar7 <- rbind(patient2_interictal_1_21_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(patient2_interictal_1_21_ar6_ar7, file = "patient2_interictal_1_21_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(patient2_interictal_1_21_ar6_ar7, file = "K:/Yet/output/patient2_interictal_1_21_ar6_ar7.csv", row.names = FALSE)


#####Patien2_ seizure = 2 interictal 22- 42#####
## Patient 2 interictal 42 preictal 18 test 150###
patient2_interictal_22_42_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                               File=character(), 
                                               User=character(), 
                                               stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 2){ # i_seizure <- 1
  for(j in 22:42){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],i_patient[2],seizure[i_seizure],j) # j <- 1
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
      
      out <- c(patient = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("patient", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],i_patient[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    patient2_interictal_22_42_ar6_ar7 <- rbind(patient2_interictal_22_42_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(patient2_interictal_22_42_ar6_ar7, file = "patient2_interictal_22_42_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(patient2_interictal_22_42_ar6_ar7, file = "K:/Yet/output/patient2_interictal_22_42_ar6_ar7.csv", row.names = FALSE)


#####Patien2_ seizure = 2 test 1-30#####
## Patient 2 interictal 42 preictal 18 test 150###
patient2_test_1_30_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                                File=character(), 
                                                User=character(), 
                                                stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 1:30){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],i_patient[2],seizure[i_seizure],j) # j <- 1
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
      
      out <- c(patient = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("patient", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],i_patient[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    patient2_test_1_30_ar6_ar7 <- rbind(patient2_test_1_30_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(patient2_test_1_30_ar6_ar7, file = "patient2_test_1_30_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(patient2_test_1_30_ar6_ar7, file = "K:/Yet/output/patient2_test_1_30_ar6_ar7.csv", row.names = FALSE)

#####Patien2_ seizure = 2 test 31-60#####
## Patient 2 interictal 42 preictal 18 test 150###
patient2_test_31_60_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                         File=character(), 
                                         User=character(), 
                                         stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 31:60){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],i_patient[2],seizure[i_seizure],j) # j <- 1
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
      
      out <- c(patient = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("patient", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],i_patient[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    patient2_test_31_60_ar6_ar7 <- rbind(patient2_test_31_60_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(patient2_test_31_60_ar6_ar7, file = "patient2_test_31_60_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(patient2_test_31_60_ar6_ar7, file = "K:/Yet/output/patient2_test_31_60_ar6_ar7.csv", row.names = FALSE)


#####Patien2_ seizure = 2 test 61-90#####
## Patient 2 interictal 42 preictal 18 test 150###
patient2_test_61_90_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                          File=character(), 
                                          User=character(), 
                                          stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 61:90){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],i_patient[2],seizure[i_seizure],j) # j <- 1
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
      
      out <- c(patient = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("patient", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],i_patient[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    patient2_test_61_90_ar6_ar7 <- rbind(patient2_test_61_90_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(patient2_test_61_90_ar6_ar7, file = "patient2_test_61_90_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(patient2_test_61_90_ar6_ar7, file = "K:/Yet/output/patient2_test_61_90_ar6_ar7.csv", row.names = FALSE)


#####Patien2_ seizure = 2 test 91-120#####
## Patient 2 interictal 42 preictal 18 test 150###
patient2_test_91_120_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                          File=character(), 
                                          User=character(), 
                                          stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 91:120){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],i_patient[2],seizure[i_seizure],j) # j <- 1
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
      
      out <- c(patient = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("patient", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],i_patient[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    patient2_test_91_120_ar6_ar7 <- rbind(patient2_test_91_120_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(patient2_test_91_120_ar6_ar7, file = "patient2_test_91_120_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(patient2_test_91_120_ar6_ar7, file = "K:/Yet/output/patient2_test_91_120_ar6_ar7.csv", row.names = FALSE)

#####Patien2_ seizure = 2 test 121-150#####
## Patient 2 interictal 42 preictal 18 test 150###
patient2_test_121_150_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                           File=character(), 
                                           User=character(), 
                                           stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 121:150){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],i_patient[2],seizure[i_seizure],j) # j <- 1
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
      
      out <- c(patient = 2, seizure = seizure[i_seizure], segment = j, channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("patient", "seizure", "segment", "channel",  "0_quantile","25_quantile", "50_quantile",
                       "75_quantile", "100_quantile", 
                       "mean", "sd", paste("pacf_", 1:7, sep = ""), "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],i_patient[2],seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    patient2_test_121_150_ar6_ar7 <- rbind(patient2_test_121_150_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(patient2_test_121_150_ar6_ar7, file = "patient2_test_121_150_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(patient2_test_121_150_ar6_ar7, file = "K:/Yet/output/patient2_test_121_150_ar6_ar7.csv", row.names = FALSE)


## Prediction Dog2####

dog2_preictal_1_42_ar6_ar7 <- read.csv("dog2_preictal_1_42_ar6_ar7.csv")
str(dog2_preictal_1_42_ar6_ar7)
dog2_interictal_1_250_ar6_ar7 <- read.csv("dog2_interictal_1_250_ar6_ar7.csv")
dog2_interictal_251_500_ar6_ar7 <- read.csv("dog2_interictal_251_500_ar6_ar7.csv")
dog2_test_1_250_ar6_ar7 <- read.csv("dog2_test_1_250_ar6_ar7.csv")
dog2_test_251_500_ar6_ar7 <- read.csv("dog2_test_251_500_ar6_ar7.csv")
dog2_test_501_750_ar6_ar7<- read.csv("dog2_test_501_750_ar6_ar7.csv")
dog2_test_751_1000_ar6_ar7 <- read.csv("dog2_test_751_1000_ar6_ar7.csv")
d1          <- rbind(dog2_preictal_1_42_ar6_ar7, 
                   dog2_interictal_1_250_ar6_ar7,
                   dog2_interictal_251_500_ar6_ar7,
                   dog2_test_1_250_ar6_ar7,
                   dog2_test_251_500_ar6_ar7,
                   dog2_test_501_750_ar6_ar7,
                   dog2_test_751_1000_ar6_ar7)


d1$ar6_aic <- d1$ar6_aic +2
d1$ar7_aic <- d1$ar7_aic +2
colnames(d1)[1:4] <- c("Dog", "Seizure", "File", "Channel")
colnames(d1)

## transpose so that all rows channels for a single file are on the same row
d <-cbind(d1[,1:4], d1[,10:11],max = d1[,9], min = d1[, 5], 
           median = d1[, 7], d1[,19:37])
write.csv(d, file ="dog2_AR6_AR7.csv")
write.csv(d, file = "K:/Yet/output/dog2_AR6_AR7.csv", row.names = FALSE)
d <- read.csv("K:/Yet/output/dog2_AR6_AR7.csv")
#first melt
melt.d = melt(d,1:4,5:ncol(d))

#then create a new name column for each channel*variable combination
melt.d$varname = paste(melt.d$variable,'.Ch',melt.d$Channel,sep='')


#then recast the data so that it works again
cast.d = dcast(melt.d,Dog+Seizure+File ~ varname,value.var = 'value')


test.d = cast.d[which(cast.d$Seizure == 'test'),]
train.d = cast.d[which(cast.d$Seizure != 'test'),]
dim(train.d)


train.d$outcome = as.factor(1*!(train.d$Seizure == 'interictal'))
# train.d$outcome = (1*!(train.d$Seizure == 'interictal'))
rf.d = train.d[,-c(1:3)]
dim(rf.d)
# Dog 2 has 542 rows and 387 columns. That's a pretty weird shape

samp = sample(1:nrow(rf.d))[1:20]
rf.d1 = rf.d[samp,]
rf.d2 = rf.d[-samp,]


# randomForest
# leave out test set first in rf.d1

# balance data with resampling preictal instances in training set
imb <- c(sum(rf.d2$outcome==1), sum(rf.d2$outcome==0))
index <- c(which(rf.d2$outcome==0), 
           sample(which(rf.d2$outcome==1), imb[2], replace=TRUE)) 
rf.d3 = rf.d[index,]

rf.fit1000<-randomForest(outcome ~ ., data=rf.d3, mtry=19, replace=TRUE, 
                         sampsize = 500, ntree=1000, importance=TRUE, 
                         na.action=na.omit) # 

# predict(rf.fit1000, rf.d1)
table(rf.d1$outcome,as.numeric(predict(rf.fit1000, rf.d1))-1>0.5)
table(rf.d1$outcome,predict(rf.fit1000, rf.d1, type = "prob")[,2]>0.5)
table(rf.d$outcome,as.numeric(predict(rf.fit1000, rf.d))-1>0.5)
predict(rf.fit1000, test.d, type = 'prob')


pred.dog2 = data.frame(cbind(test.d[,1:3], predict(rf.fit1000,test.d,type='prob')))

names(pred.dog2)[4:5] = c("prob.0","prob.1")

write.csv(file='K:/submissions/sept22/Yet_prediction_Dog2.csv',pred.dog2[,c(1:3,5)],row.names = FALSE)


###Try CV 1 #####
cv1_dog2 <- laply(1:nrow(rf.d), function(samp){#samp = 1

rf.d1 = rf.d[samp,]
rf.d2 = rf.d[-samp,]


# randomForest
# leave out test set first in rf.d1

# balance data with resampling preictal instances in training set
imb <- c(sum(rf.d2$outcome==1), sum(rf.d2$outcome==0))
index <- c(which(rf.d2$outcome==0), 
           sample(which(rf.d2$outcome==1), imb[2], replace=TRUE)) 
rf.d3 = rf.d[index,]

rf.fit1000<-randomForest(outcome ~ ., data=rf.d3, mtry=19, replace=TRUE, 
                         sampsize = 500, ntree=1000, importance=TRUE, 
                         na.action=na.omit) # 

 predict(rf.fit1000, rf.d1,type='prob')[,2]
#table(rf.d1$outcome,as.numeric(predict(rf.fit1000, rf.d1))-1>0.5)
}
)


table(rf.d$outcome,as.numeric(predict(rf.fit1000, rf.d))-1>0.5)
predict(rf.fit1000, test.d, type = 'prob')


pred.dog2 = data.frame(cbind(test.d[,1:3], predict(rf.fit1000,test.d,type='prob')))

names(pred.dog2)[4:5] = c("prob.0","prob.1")

write.csv(file='K:/submissions/sept22/Yet_prediction_Dog2.csv',pred.dog2[,c(1:3,5)],row.names = FALSE)

####################



##


## Prediction Patient2####

patient2_preictal_1_18_ar6_ar7 <- read.csv("patient2_preictal_1_18_ar6_ar7.csv")

patient2_interictal_1_21_ar6_ar7 <- read.csv("patient2_interictal_1_21_ar6_ar7.csv")
patient2_interictal_22_42_ar6_ar7 <- read.csv("patient2_interictal_22_42_ar6_ar7.csv")
patient2_test_1_30_ar6_ar7 <- read.csv("patient2_test_1_30_ar6_ar7.csv")
patient2_test_31_60_ar6_ar7 <- read.csv("patient2_test_31_60_ar6_ar7.csv")
patient2_test_61_90_ar6_ar7<- read.csv("patient2_test_61_90_ar6_ar7.csv")
patient2_test_91_120_ar6_ar7 <- read.csv("patient2_test_91_120_ar6_ar7.csv")
patient2_test_121_150_ar6_ar7 <- read.csv("patient2_test_121_150_ar6_ar7.csv")

d1          <- rbind(patient2_preictal_1_18_ar6_ar7, 
                     patient2_interictal_1_21_ar6_ar7 ,
                     patient2_interictal_22_42_ar6_ar7,
                     patient2_test_1_30_ar6_ar7 ,
                     patient2_test_31_60_ar6_ar7,
                     patient2_test_61_90_ar6_ar7,
                     patient2_test_91_120_ar6_ar7,
                     patient2_test_121_150_ar6_ar7
                     )


d1$ar6_aic <- d1$ar6_aic +2
d1$ar7_aic <- d1$ar7_aic +2
colnames(d1)[1:4] <- c("Dog", "Seizure", "File", "Channel")
colnames(d1)

## transpose so that all rows channels for a single file are on the same row
d <-cbind(d1[,1:4], d1[,10:11],max = d1[,9], min = d1[, 5], 
          median = d1[, 7], d1[,19:37])
write.csv(d, file ="patient_2_AR6_AR7.csv")
write.csv(d, file = "K:/Yet/output/patient2_AR6_AR7.csv", row.names = FALSE)

#first melt
melt.d = melt(d,1:4,5:ncol(d))

#then create a new name column for each channel*variable combination
melt.d$varname = paste(melt.d$variable,'.Ch',melt.d$Channel,sep='')


#then recast the data so that it works again
cast.d = dcast(melt.d,Dog+Seizure+File ~ varname,value.var = 'value')


test.d = cast.d[which(cast.d$Seizure == 'test'),]
train.d = cast.d[which(cast.d$Seizure != 'test'),]
dim(train.d)

train.d$outcome = as.factor(1*!(train.d$Seizure == 'interictal'))
# train.d$outcome = (1*!(train.d$Seizure == 'interictal'))
rf.d = train.d[,-c(1:3)]
dim(rf.d)
# Dog 2 has 542 rows and 387 columns. That's a pretty weird shape

samp = sample(1:nrow(rf.d))[1:10]
rf.d1 = rf.d[samp,]
rf.d2 = rf.d[-samp,]
dim(rf.d2)

# randomForest
# leave out test set first in rf.d1

# balance data with resampling preictal instances in training set
imb <- c(sum(rf.d2$outcome==1), sum(rf.d2$outcome==0))
imb
index1 <- c(which(rf.d2$outcome==0), 
           sample(which(rf.d2$outcome==1), imb[2], replace=TRUE)) 
index <- sample(index1,1000, replace = TRUE )
rf.d3 = rf.d[index,]

rf.fit1000<-randomForest(outcome ~ ., data=rf.d3, mtry=19, replace=TRUE, 
                         sampsize = 500, ntree=1000, importance=TRUE, 
                         na.action=na.omit) # 

# predict(rf.fit1000, rf.d1)
table(rf.d1$outcome,as.numeric(predict(rf.fit1000, rf.d1))-1>0.5)
table(rf.d$outcome,as.numeric(predict(rf.fit1000, rf.d))-1>0.5)
predict(rf.fit1000, test.d, type = 'prob')


pred.patient2 = data.frame(cbind(test.d[,1:3], predict(rf.fit1000,test.d,type='prob')))

names(pred.patient2)[4:5] = c("prob.0","prob.1")

write.csv(file='K:/submissions/sept22/Yet_prediction_Patient2.csv',pred.patient2[,c(1:3,5)],row.names = FALSE)





n_segment <- c(18, 50, 195)
seizure <- c("preictal", "interictal", "test")
type <- c("Dog", "Patient")
i_type <- c(1:2)
i_dog <- c(1,2,3,4,5)
i_patient <- c(1,2)

#####Patient1_ seizure = 1 preictal 1-18#####
## Patient 1 interictal 50 preictal 18 test 195###
Patient1_preictal_1_18_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                             File=character(), 
                                             User=character(), 
                                             stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 1){ # i_seizure <- 1
  for(j in 1:n_segment[i_seizure]){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],1,seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    #str(d)
    res <- ldply(1:n_channel, function(i){ # i <- 5
      
      Dog.stats <- c(mean(dat_mat[i,]),
                     sd(dat_mat[i,]),
                     max(dat_mat[i,]),
                     min(dat_mat[i,]),
                     median(dat_mat[i,]))
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*7 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*8 - 2*d.fit7$loglikelihood)
      
      out <- c(Patient = 1, Seizure = seizure[i_seizure], File = j, Channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("Patient", "Seizure", "File", "Channel",  
                       "mean", "sd", "max","min", "median", "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],1,seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    Patient1_preictal_1_18_ar6_ar7 <- rbind(Patient1_preictal_1_18_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(Patient1_preictal_1_18_ar6_ar7, file = "Patient1_preictal_1_18_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(Patient1_preictal_1_18_ar6_ar7, file = "K:/Yet/output/Patient1_preictal_1_18_ar6_ar7.csv", row.names = FALSE)

#####Patient1_ seizure = 2 interictal 1- 25#####
## Patient 1 interictal 50 preictal 18 test 195###
Patient1_interictal_1_25_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                               File=character(), 
                                               User=character(), 
                                               stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 2){ # i_seizure <- 1
  for(j in 1:25){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],1,seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    #str(d)
    res <- ldply(1:n_channel, function(i){ # i <- 5
      
      Dog.stats <- c(mean(dat_mat[i,]),
                     sd(dat_mat[i,]),
                     max(dat_mat[i,]),
                     min(dat_mat[i,]),
                     median(dat_mat[i,]))
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*7 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*8 - 2*d.fit7$loglikelihood)
      
      out <- c(Patient = 1, Seizure = seizure[i_seizure], File = j, Channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("Patient", "Seizure", "File", "Channel",  
                       "mean", "sd", "max","min", "median", "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],1,seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    Patient1_interictal_1_25_ar6_ar7 <- rbind(Patient1_interictal_1_25_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(Patient1_interictal_1_25_ar6_ar7, file = "Patient1_interictal_1_25_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(Patient1_interictal_1_25_ar6_ar7, file = "K:/Yet/output/Patient1_interictal_1_25_ar6_ar7.csv", row.names = FALSE)


#####Patient1_ seizure = 2 interictal 25- 50#####
## Patient 1 interictal 50 preictal 18 test 195###
Patient1_interictal_25_50_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                                File=character(), 
                                                User=character(), 
                                                stringsAsFactors=FALSE) 



n_segment <- c(18, 50, 195)
seizure <- c("preictal", "interictal", "test")
type <- c("Dog", "Patient")
i_type <- c(1:2)
i_dog <- c(1,2,3,4,5)
i_patient <- c(1,2)


pm1 <- proc.time()
for(i_seizure in 2){ # i_seizure <- 1
  for(j in 25:50){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],1,seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    #str(d)
    res <- ldply(1:n_channel, function(i){ # i <- 5
      
      Dog.stats <- c(mean(dat_mat[i,]),
                     sd(dat_mat[i,]),
                     max(dat_mat[i,]),
                     min(dat_mat[i,]),
                     median(dat_mat[i,]))
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*7 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*8 - 2*d.fit7$loglikelihood)
      
      out <- c(Patient = 1, Seizure = seizure[i_seizure], File = j, Channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("Patient", "Seizure", "File", "Channel",  
                       "mean", "sd", "max","min", "median", "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],1,seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    Patient1_interictal_25_50_ar6_ar7 <- rbind(Patient1_interictal_25_50_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(Patient1_interictal_25_50_ar6_ar7, file = "Patient1_interictal_25_50_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(Patient1_interictal_25_50_ar6_ar7, file = "K:/Yet/output/Patient1_interictal_25_50_ar6_ar7.csv", row.names = FALSE)


#####Patient1_ seizure = 2 test 1-30#####
## Patient 1 interictal 50 preictal 18 test 195###



n_segment <- c(18, 50, 195)
seizure <- c("preictal", "interictal", "test")
type <- c("Dog", "Patient")
i_type <- c(1:2)
i_dog <- c(1,2,3,4,5)
i_patient <- c(1,2)

Patient1_test_1_40_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                         File=character(), 
                                         User=character(), 
                                         stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 1:40){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],1,seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    #str(d)
    res <- ldply(1:n_channel, function(i){ # i <- 5
      
      Dog.stats <- c(mean(dat_mat[i,]),
                     sd(dat_mat[i,]),
                     max(dat_mat[i,]),
                     min(dat_mat[i,]),
                     median(dat_mat[i,]))
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*7 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*8 - 2*d.fit7$loglikelihood)
      
      out <- c(Patient = 1, Seizure = seizure[i_seizure], File = j, Channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("Patient", "Seizure", "File", "Channel",  
                       "mean", "sd", "max","min", "median", "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],1,seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    Patient1_test_1_40_ar6_ar7 <- rbind(Patient1_test_1_40_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(Patient1_test_1_40_ar6_ar7, file = "Patient1_test_1_40_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(Patient1_test_1_40_ar6_ar7, file = "K:/Yet/output/Patient1_test_1_40_ar6_ar7.csv", row.names = FALSE)

#####Patient1_ seizure = 2 test 31-60#####
## Patient 1 interictal 50 preictal 18 test 195###



n_segment <- c(18, 50, 195)
seizure <- c("preictal", "interictal", "test")
type <- c("Dog", "Patient")
i_type <- c(1:2)
i_dog <- c(1,2,3,4,5)
i_patient <- c(1,2)

Patient1_test_41_80_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                          File=character(), 
                                          User=character(), 
                                          stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 41:80){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],1,seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    #str(d)
    res <- ldply(1:n_channel, function(i){ # i <- 5
      
      Dog.stats <- c(mean(dat_mat[i,]),
                     sd(dat_mat[i,]),
                     max(dat_mat[i,]),
                     min(dat_mat[i,]),
                     median(dat_mat[i,]))
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*7 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*8 - 2*d.fit7$loglikelihood)
      
      out <- c(Patient = 1, Seizure = seizure[i_seizure], File = j, Channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("Patient", "Seizure", "File", "Channel",  
                       "mean", "sd", "max","min", "median", "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],1,seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    Patient1_test_41_80_ar6_ar7 <- rbind(Patient1_test_41_80_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(Patient1_test_41_80_ar6_ar7, file = "Patient1_test_41_80_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(Patient1_test_41_80_ar6_ar7, file = "K:/Yet/output/Patient1_test_41_80_ar6_ar7.csv", row.names = FALSE)


#####Patient1_ seizure = 2 test 61-90#####
## Patient 1 interictal 50 preictal 18 test 195###



n_segment <- c(18, 50, 195)
seizure <- c("preictal", "interictal", "test")
type <- c("Dog", "Patient")
i_type <- c(1:2)
i_dog <- c(1,2,3,4,5)
i_patient <- c(1,2)

Patient1_test_81_120_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                          File=character(), 
                                          User=character(), 
                                          stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 81:120){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],1,seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    #str(d)
    res <- ldply(1:n_channel, function(i){ # i <- 5
      
      Dog.stats <- c(mean(dat_mat[i,]),
                     sd(dat_mat[i,]),
                     max(dat_mat[i,]),
                     min(dat_mat[i,]),
                     median(dat_mat[i,]))
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*7 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*8 - 2*d.fit7$loglikelihood)
      
      out <- c(Patient = 1, Seizure = seizure[i_seizure], File = j, Channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("Patient", "Seizure", "File", "Channel",  
                       "mean", "sd", "max","min", "median", "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],1,seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    Patient1_test_81_120_ar6_ar7 <- rbind(Patient1_test_81_120_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(Patient1_test_81_120_ar6_ar7, file = "Patient1_test_81_120_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(Patient1_test_81_120_ar6_ar7, file = "K:/Yet/output/Patient1_test_81_120_ar6_ar7.csv", row.names = FALSE)


#####Patient1_ seizure = 2 test 91-120#####
## Patient 1 interictal 50 preictal 18 test 195###



n_segment <- c(18, 50, 195)
seizure <- c("preictal", "interictal", "test")
type <- c("Dog", "Patient")
i_type <- c(1:2)
i_dog <- c(1,2,3,4,5)
i_patient <- c(1,2)

Patient1_test_121_160_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                           File=character(), 
                                           User=character(), 
                                           stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 121:160){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],1,seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    #str(d)
    res <- ldply(1:n_channel, function(i){ # i <- 5
      
      Dog.stats <- c(mean(dat_mat[i,]),
                     sd(dat_mat[i,]),
                     max(dat_mat[i,]),
                     min(dat_mat[i,]),
                     median(dat_mat[i,]))
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*7 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*8 - 2*d.fit7$loglikelihood)
      
      out <- c(Patient = 1, Seizure = seizure[i_seizure], File = j, Channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("Patient", "Seizure", "File", "Channel",  
                       "mean", "sd", "max","min", "median", "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],1,seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    Patient1_test_121_160_ar6_ar7 <- rbind(Patient1_test_121_160_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(Patient1_test_121_160_ar6_ar7, file = "Patient1_test_121_160_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(Patient1_test_121_160_ar6_ar7, file = "K:/Yet/output/Patient1_test_121_160_ar6_ar7.csv", row.names = FALSE)

#####Patient1_ seizure = 2 test 121-150#####
## Patient 1 interictal 50 preictal 18 test 195###



n_segment <- c(18, 50, 195)
seizure <- c("preictal", "interictal", "test")
type <- c("Dog", "Patient")
i_type <- c(1:2)
i_dog <- c(1,2,3,4,5)
i_patient <- c(1,2)

Patient1_test_161_195_ar6_ar7 <- data.frame(Date=as.Date(character()),
                                            File=character(), 
                                            User=character(), 
                                            stringsAsFactors=FALSE) 

pm1 <- proc.time()
for(i_seizure in 3){ # i_seizure <- 1
  for(j in 161:195){
    pm2 <- proc.time()
    d <- get.seizure.file(type[2],1,seizure[i_seizure],j) # j <- 1
    dat_mat <- d[[1]][[1]]
    n_channel <- nrow(dat_mat)
    #str(d)
    res <- ldply(1:n_channel, function(i){ # i <- 5
      
      Dog.stats <- c(mean(dat_mat[i,]),
                     sd(dat_mat[i,]),
                     max(dat_mat[i,]),
                     min(dat_mat[i,]),
                     median(dat_mat[i,]))
      
      d.fit7 <- FitAR(dat_mat[i,], 7)
      d.fit6 <- FitAR(dat_mat[i,], 6)  
      
      AIC6.stats = c(d.fit6$muHat, d.fit6$phiHat, d.fit6$sigsqHat, 2*7 - 2*d.fit6$loglikelihood)
      AIC7.stats = c(d.fit7$muHat, d.fit7$phiHat, d.fit7$sigsqHat, 2*8 - 2*d.fit7$loglikelihood)
      
      out <- c(Patient = 1, Seizure = seizure[i_seizure], File = j, Channel = i, 
               Dog.stats,AIC6.stats,AIC7.stats)
      names(out)  <- c("Patient", "Seizure", "File", "Channel",  
                       "mean", "sd", "max","min", "median", "ar6_mu",
                       paste("ar", 6,"_ar", 1:6, sep = ""),
                       "ar6_sigma2", "ar6_aic", "ar7_mu",
                       paste("ar", 7,"_ar", 1:7, sep = ""),
                       "ar7_sigma2", "ar7_aic" )
      print(paste(type[2],1,seizure[i_seizure],j,"channel", i,  sep = "_"))
      return(out)
      
    })
    proc.time()-pm2
    
    Patient1_test_161_195_ar6_ar7 <- rbind(Patient1_test_161_195_ar6_ar7, res)
  }
}

proc.time()-pm1
write.csv(Patient1_test_161_195_ar6_ar7, file = "Patient1_test_161_195_ar6_ar7.csv", row.names = FALSE)

#dog2_preictal_1_42_ar6_ar7 <- dog2_out1
write.csv(Patient1_test_161_195_ar6_ar7, file = "K:/Yet/output/Patient1_test_161_195_ar6_ar7.csv", row.names = FALSE)
