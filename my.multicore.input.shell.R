#!/usr/bin/env Rscript
setwd("~/Desktop/test/source")


#for (x1 in 1:9){
runs = 100
siteA = 100
freqX = 0
freqY = 0
neut = seq(0.01,125,0.1)
idx = which(neut > 0)
neut = neut[idx]



my.stat <- c(); 
my.summary.stat <- c();
sink(paste('my.model.input',1,sep=''));
cat('1')# number of runs
cat('\n')
cat('8000')#population size = 2N
cat('\n')
cat('0.01')#mutation rate
cat('\n')
cat('0.01')#recombination rate
cat('\n')
cat('0')#migration rate (4Nm)
cat('\n')
cat('11')#Number of epochs
cat('\n')
cat('1000 2000 3000 5000 7000 9000 12000 15000 20000 25000')#Epoch Breakpoints
cat('\n')
cat('1.0 0.90 0.80 0.70 0.60 0.50 0.40 0.30 0.15 0.20 0.40')#Scaling Factors
#cat('1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0')
#cat('1.0 0.20 0.10 0.30 0.50 0.50 0.50 0.50 0.5 0.5 0.5')#Scaling Factors
cat('\n')
cat('0 0')# age of siteA, Y chr: set to 0 assumes no sweep, infinite age
cat('\n')
cat('0')#SDR site >= 0 
cat('\n')
cat(siteA)#siteA site >= 0
cat('\n')
cat('1.0')#male female recombination ratio
cat('\n')
cat(freqX)
cat('\n')
cat(freqY)
cat('\n')
cat('0')# freqA in Y previous
cat('\n')
cat('0.001') #small number to prevent simulation very close to SDR
#cat('0')#number of carriers
cat('\n')
cat('2')#sampling scheme
cat('\n')
cat(paste(neut,collapse = ' ')) # Neutral Sites
cat('\n')
for(i in 1:16){
cat( paste(c(0,0,0),collapse = ' ') )#8chrX
cat('\n')
};
for(i in 1:16){
cat( paste(c(0,1,0),collapse = ' ') )
cat('\n')
}; #8chrY
sink();


system(paste("./run_simulation.sh",runs,siteA,freqX*100,freqY*100,sep=' '));



system(paste("./test.sh",runs,siteA,freqX*100,freqY*100,sep=' '));
















