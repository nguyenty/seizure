require(R.matlab)
require(plyr)
require(ggplot2)
require(reshape2)
library(reshape)
#install.packages("FitAR")
library(FitAR)
#install.packages("randomForest")
library(randomForest)

get.seizure.file = function(dog.type,dog.count,seizure.type,file.count){
  data.dir = '/run/user/1000/gvfs/smb-share:server=researchfiles.iastate.edu,share=las$/STAT/KaggleDataComp/data/raw_data/'
  #    data.dir = '/media/yet/DATA/seizure'
  #   data.dir <- "K:/data/raw_data/"
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

dog.type = "Dog"
dog.count = 1:5
seizure.type = c("preictal", "interictal", "test")
file.count =  matrix(c(480,24,502,500,42,1000,1440,72,907,804,97,990,450,30,191),ncol=3,byrow=T)
file.count
## run to find the mean of data of each second and see how they are#####
