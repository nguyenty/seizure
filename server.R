#--------------------------------------**--------------------------------------#
#  File Name: server.R
#  Purpose:
#
#  Creation Date: 11-09-2014
#  Last Modified: Fri Sep 12 17:49:55 2014
#  Created By:
#
#--------------------------------------**--------------------------------------#
#
#  FORTRAN and C: 
#  source('~/R/shlib/C_FORTRAN.shlib.r')
#  .Fortran("subroutine name",as.integer(input1),as.double(input2), etc)
#
#install.packages("R.matlab")
require(R.matlab)
require(plyr)
require(ggplot2)
require(reshape2)
require(shiny)
require(MASS)

#-- Data is stored in 


get.seizure.file = function(dog.type,dog.count,seizure.type,file.count){
   data.dir = '/Volumes/KaggleDataComp/data/raw_data/'
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

d.mat1 = get.seizure.file('Dog',1,'preictal',2)
d.mat2 = get.seizure.file('Dog',4,'preictal',2)
d.mat3 = get.seizure.file('Dog',4,'interictal',105)

unlist(d.mat1[[1]][[4]])
unlist(d.mat2[[1]][[4]])
== 
unlist(d.mat3[[1]][[4]])

d.mat[[1]][[3]]

str(d.mat)

d = data.frame(t(d.mat[[1]][[1]]))
names(d) = paste('ch',1:16,sep='')
d$t = 600*(1:nrow(d)/nrow(d))

p1 = ggplot()+geom_line(data=d[1:1000,],aes(t,ch1))
p2 = ggplot()+geom_line(data=d,aes(t,abs(ch2)))
p3 = ggplot()+geom_line(data=d,aes(t,ch3))
p4 = ggplot()+geom_line(data=d,aes(t,ch4))
p5 = ggplot()+geom_line(data=d,aes(t,ch5))
p6 = ggplot()+geom_line(data=d,aes(t,ch6))
p7 = ggplot()+geom_line(data=d,aes(t,ch7))
p8 = ggplot()+geom_line(data=d,aes(t,ch8))
p9 = ggplot()+geom_line(data=d,aes(t,ch9))
p10 = ggplot()+geom_line(data=d,aes(t,ch10))
p11 = ggplot()+geom_line(data=d,aes(t,ch11))
p12 = ggplot()+geom_line(data=d,aes(t,ch12))
p13 = ggplot()+geom_line(data=d,aes(t,ch13))
p14 = ggplot()+geom_line(data=d,aes(t,ch14))
p15 = ggplot()+geom_line(data=d,aes(t,ch15))
p16 = ggplot()+geom_line(data=d,aes(t,ch16))

#BOX jenkins
p10
ts = arma(d$ch10,order=c(20,0))
summary(ts)

p+geom_line(data=d,aes(t,ch2))
p
freq1 = d.mat[[1]][[2]]
freq2 = d.mat[[1]][[3]]
channel_names = matrix(unlist(d.mat[[1]][[4]]),nrow=16)

d = cbind(channel_names,d)

names(d) = c('ch', paste('x',1:(ncol(d)-1),sep=''))
ggplot()+geom_line(aes(x=1:(ncol(d)-1), y = d[1,2:ncol(d)]))

d.melt = melt(d.use, id.vars = 'ch', measure.vars = names(d.use)[2:ncol(d.use)])
head(d.melt)

require(tseries)


str(d.mat)

p1 = qplot(1:ncol(d),d[1,],geom='line')
p2 = qplot(1:ncol(d),d[2,],geom='line')
p3 = qplot(1:ncol(d),d[3,],geom='line')
p4 = qplot(1:ncol(d),d[4,],geom='line')
p5 = qplot(1:ncol(d),d[5,],geom='line')
p6 = qplot(1:ncol(d),d[6,],geom='line')
p7 = qplot(1:ncol(d),d[7,],geom='line')
p8 = qplot(1:ncol(d),d[8,],geom='line')
p9 = qplot(1:ncol(d),d[9,],geom='line')
p10 = qplot(1:ncol(d),d[10,],geom='line')
p11 = qplot(1:ncol(d),d[11,],geom='line')
p12 = qplot(1:ncol(d),d[12,],geom='line')
p13 = qplot(1:ncol(d),d[13,],geom='line')
p14 = qplot(1:ncol(d),d[14,],geom='line')
p15 = qplot(1:ncol(d),d[15,],geom='line')

p16 = qplot(1:ncol(d),d[16,],geom='line')

p16

qplot(1:ncol(d),abs(d[1,]),geom='line')

d = data.frame(cbind(1:16,d))

names(d.melt) = c('ch','index','amp')



qplot(data = d.melt,as.numeric(index),amp,geom='line',group=ch,color=ch)


str(d.mat[[1]])
#0804
#0097
#0990
.mat

require(ggplot2)
require(plyr)
require(tseries)

d = 
read.matrix(
