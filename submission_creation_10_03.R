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
# win.wd = 'K:/'
# mac.wd = '/Volumes/las$/STAT/KaggleDataComp/'
# Kdrive =win.wd 
# setwd(Kdrive)

#read in the predictions template
sample = read.csv('sample_submission.csv')
dogs2files = read.csv('fileID2filename.csv')


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# PREDICTION 1: Yet's Dog 2, Ran's Dog 4, Patients at .49 #
dog1.1 = read.csv('Dog_1_pred.csv')
dog2.1 = read.csv('Dog_2_pred.csv')
dog3.1 = read.csv('Dog_3_pred.csv')
dog4.1 = read.csv('Dog_4_pred.csv')
dog5.1 = read.csv('Dog_5_pred.csv')
patient1.1 = read.csv('Patient_1_pred.csv')
patient2.1 = read.csv('Patient_2_pred.csv')
for(i in 1:5){eval(parse(text = paste("names(dog",i,".1) = c('Dog','Seizure','File','prob')",sep=''))) }
dogs = rbind(dog1.1,dog2.1,dog3.1,dog4.1,dog5.1)
dogs$type = 'Dog'

dog.pred1 = merge(dogs2files,dogs,by.x=c('type','count','segment','number'),by.y=c('type','Dog','Seizure','File'),all.x=TRUE)

names(dog.pred1)[5:6] = c('clip','preictal')
dog.pred1$preictal[is.na(dog.pred1$preictal)] = .49

dog.pred1[-c(1:3590),6] <- c(patient1.1[, 4], patient2.1[, 4])
qplot(count,preictal,geom='boxplot',data=dog.pred1,color=interaction(count,type))

write.csv(dog.pred1[,c('clip','preictal')],file='submission13.csv',row.names=FALSE)
