
## Squarem implementation of the Poisson Random Batch model

library(SQUAREM)
library(optimx)
library(gtools)
library(parallel)
library(lineprof)

#######  Source the important files

source('poisson_loglik.R')
source('simplex_functions.R')
source('estimation_poisson_batch.R')
source('simulate_counts_poisson.R')

###  If you want to simulate from the batch model : try running the example_simulation.R file

source('example_simulation.R');

##### From now on we assume we only have the counts matrix available ##########################

source('poisson_topic_loglink.R')

### apply the main function for the modeling

res <- Poisson_topic.loglink(counts,n_clus=4,lab_batch)
barplot(t(res$omega),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
