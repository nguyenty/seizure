require(R.matlab)
require(plyr)
require(ggplot2)
require(reshape2)

#install.packages("FitAR")
library(FitAR)
#install.packages("randomForest")
library(randomForest)

get.seizure.file = function(dog.type,dog.count,seizure.type,file.count){
#   data.dir = '/run/user/1000/gvfs/smb-share:server=researchfiles.iastate.edu,share=las$/STAT/KaggleDataComp/data/raw_data/'
  #    data.dir = '/media/yet/DATA/seizure'
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


dog.type = c("Dog", "Patient")
dog.count = c(5, 2) # 5 dogs and 2 patients
seizure.type = c( "interictal", "preictal", "test")
file.count =list(matrix(c(480,24,502,
                       500,42,1000,
                       1440,72,907,
                       804,97,990,
                       450,30,191),ncol=3,byrow=T),
                 matrix(c(50, 18, 195, 
                          42, 18, 150), 
                        ncol = 3, byrow = T))



## run to find the mean of data of each second and see how they are#####
## Features are pacf of each channel #####
pm1 <- proc.time()
for (i in 1:2){ # dog.type i <- 2
  for (j in 1:dog.count[i]){ # dog.count j <- 1
    d <- data.frame(Date=as.Date(character()),File=character(), User=character(), stringsAsFactors=FALSE) 
    for (k in 1:3){ # seizure.type k <- 1
      for (l in 1:file.count[[i]][j, k]){ # file.count l <- 1
        dat <- get.seizure.file(dog.type[i], j, seizure.type[k], l)
        dat_mat <- dat[[1]][[1]]
        n_channel <- nrow(dat_mat)
        res <- ldply(1:n_channel, function(m){ # m <- 1
          
          Dog.stats <- c(mean(dat_mat[m,]),
                         sd(dat_mat[m,]),
                         median(dat_mat[m,]),
                         max(dat_mat[m,]), min(dat_mat[m,]),
                         pacf(dat_mat[m,], plot = F)$acf[1:10])
          
          out <- c(j, seizure = seizure.type[k], segment = l, channel = m, 
                   Dog.stats)
          names(out)  <- c(dog.type[i], "Seizure", "File", "Channel", 
                           "mean", "sd", "median", "max","min", paste("pacf_", 1:10, sep = ""))
          print(paste(dog.type[i],j,seizure.type[k],"segment",l,"channel", m,  sep = "_"))
          return(out)
        })
      d <- rbind(d, res)  
      }
    }
    path1 <- paste0(dog.type[i], "_", j, "_pacf10.csv" )
    path2 <- paste0("K:/Yet/output/",dog.type[i], "_", j, "_pacf10.csv" )
    write.csv(d, file = path1, row.names = F)
    write.csv(d, file = path2, row.names = F)
    #first melt
    melt.d = melt(d,1:4,5:ncol(d))
    
    #then create a new name column for each channel*variable combination
    melt.d$varname = paste(melt.d$variable,'.Ch',melt.d$Channel,sep='')
    
    
    #then recast the data so that it works again
    cast.d = dcast(melt.d,dog.type[i]+Seizure+File ~ varname,value.var = 'value')
    
    
    test.d = cast.d[which(cast.d$Seizure == 'test'),]
    train.d = cast.d[which(cast.d$Seizure != 'test'),]
    
    
    train.d$outcome = as.factor(1*!(train.d$Seizure == 'interictal'))
    # train.d$outcome = (1*!(train.d$Seizure == 'interictal'))
    rf.d = train.d[,-c(1:3)]
    
    # randomForest
    # leave out test set first in rf.d1
    
    # balance data with resampling preictal instances in training set
    imb <- c(sum(rf.d$outcome==1), sum(rf.d$outcome==0))
    if(i ==1){
      index <- c(which(rf.d$outcome==0), 
                  sample(which(rf.d$outcome==1), imb[2], replace=TRUE)) 
    }else{
      index1 <- c(which(rf.d$outcome==0), 
                  sample(which(rf.d$outcome==1), imb[2], replace=TRUE)) 
      index <- sample(index1, 1000, replace = TRUE )
    }
    
    
    rf.d3 = rf.d[index,]
    rf.fit1000<-randomForest(outcome ~ ., data=rf.d3,  replace=TRUE, 
                             sampsize = 500, ntree=1000, importance=TRUE, 
                             na.action=na.omit) 
    pred.d = data.frame(cbind(test.d[,1:3], predict(rf.fit1000,test.d,type='prob')))
    
    names(pred.d)[4:5] = c("prob.0","prob.1")
    path3 <- paste0(dog.type[i], "_", j, "_pred.csv" )
    write.csv(pred.d[,c(1:3, 5)], file = path3, row.names = F)
    
  }
}

proc.time() -pm1
