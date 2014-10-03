require(plyr)
require(ggplot2)
require(reshape2)
#library(reshape)
#install.packages("FitAR")
#library(FitAR)
#install.packages("randomForest")
library(randomForest)

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

for(i in 1:2){ # i <- 1
  for(j in 1:dog.count[i]){ # j <- 1
    path1 <- paste0(dog.type[i], "_", j, "_pacf10.csv" )
    
    d <- read.csv(file= path1)
    dim(d)
    head(d)
    #first melt
    melt.d = melt(d,1:4,5:ncol(d))
    
    #then create a new name column for each channel*variable combination
    melt.d$varname = paste(melt.d$variable,'.Ch',melt.d$Channel,sep='')
    
    
    #then recast the data so that it works again
    if(i ==1){cast.d =  dcast(melt.d,Dog+Seizure+File ~ varname,value.var = 'value')
    }else{
      cast.d = dcast(melt.d,Patient+Seizure+File ~ varname,value.var = 'value')
    }
     
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
