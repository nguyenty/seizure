#--------------------------------------**--------------------------------------#
#  File Name:
#  Purpose:
#
#  Creation Date: 24-09-2014
#  Last Modified: Wed Sep 24 23:35:10 2014
#  Created By:
#
#--------------------------------------**--------------------------------------#
#
#  FORTRAN and C: 
#  source('~/R/shlib/C_FORTRAN.shlib.r')
#  .Fortran("subroutine name",as.integer(input1),as.double(input2), etc)
#

require('reshape2',quietly=TRUE) #for casting to convert data from long to wide

require("ggplot2", quietly = TRUE)
theme_set(theme_bw())

require('plyr',quietly=TRUE)
require('rpart',quietly=TRUE)
require('earth',quietly=TRUE)
require('fields',quietly=TRUE)
require('nnet',quietly=TRUE)
require('neuralnet',quietly=TRUE)
require('randomForest',quietly=TRUE)
require('ggmap',quietly=TRUE)
require('party')

#set working directory
win.wd = 'K:/'
mac.wd = '/Volumes/las$/STAT/KaggleDataComp/'
Kdrive = mac.wd 
setwd(Kdrive)

#read in the predictions template
sample = read.csv('./submissions/sample_submission.csv')
dogs2files = read.csv('./submissions/fileID2filename.csv')


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# PREDICTION 1: Yet's Dog 2, Ran's Dog 4, Patients at .49 #
dog1.1 = read.csv('./submissions/sept22/ian_prediction_Dog1.csv')
dog2.1 = read.csv('./submissions/sept22/Yet_prediction_Dog2.csv')
dog3.1 = read.csv('./submissions/sept22/Xiaojun_prediction_Dog3.csv')
dog4.1 = read.csv('./submissions/sept22/Ran_prediction_Dog4.csv')
dog5.1 = read.csv('./submissions/sept22/Yihua_prediction_Dog5.csv')
for(i in 1:5){eval(parse(text = paste("names(dog",i,".1) = c('Dog','Seizure','File','prob')",sep=''))) }
dogs = rbind(dog1.1,dog2.1,dog3.1,dog4.1,dog5.1)
dogs$type = 'Dog'

dog.pred1 = merge(dogs2files,dogs,by.x=c('type','count','segment','number'),by.y=c('type','Dog','Seizure','File'),all.x=TRUE)

names(dog.pred1)[5:6] = c('clip','preictal')
dog.pred1$preictal[is.na(dog.pred1$preictal)] = .49

head(dog.pred1)

qplot(count,preictal,geom='boxplot',data=dog.pred1,color=interaction(count,type))

write.csv(dog.pred1[,c('clip','preictal')],file='./submissions/sept24/submission1.csv',row.names=FALSE)
