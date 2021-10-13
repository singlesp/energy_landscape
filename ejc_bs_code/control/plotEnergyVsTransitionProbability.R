args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]
c <- args[4]

library(ggplot2)
library(R.matlab)
library(RColorBrewer)
library(lm.beta)
library(reshape2)
library(viridis)

numClusters = 4; #make sure this is reflected in the file name below!!!!!!!
a<-readMat("/Users/sps253/Documents/brain_states-master/results/example/TransProbsData_bp12_k4.mat") #transprob mats
b<-readMat("/Users/sps253/Documents/brain_states-master/results/example/EnergyData_bp12_k4.mat") #group avg E_calcs
c<-readMat("/Users/sps253/Documents/brain_states-master/sch454_subjcentroids.mat") #subject specific E_calcs
#b<-readMat("/Users/sps253/Documents/PTSD_MDMA1/clean/scatter_data.mat")
caps<-b$A[1:8,3]
limdiff<-b$limdiff[1:8,3]


masterdir <- paste(basedir,'results/',name_root,'/',sep='')
source(paste('/Users/sps253/Documents/brain_states-master/code/plottingfxns/plottingfxns.R',sep=''))
source(paste('/Users/sps253/Documents/brain_states-master/code/miscfxns/statfxns.R',sep=''))

clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])

clusterNames <- a$clusters
clusterColors <- getClusterColors(numClusters)
RNcolors <- c('#005C9F','#FF8400')  

savedir <- paste(masterdir,'analyses/control_energy/EvsTP_plots/',sep = '')
dir.create(path = savedir,recursive = TRUE)

method = 'spearman'

### load transition probabilities ###

restTP <- a$LSDTransitionProbabilityMats
restTP <- colMeans(restTP) # group average
nbackTP <- a$PLTransitionProbabilityMats
nbackTP <- colMeans(nbackTP,na.rm=T)

onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1))) # construct indices to select only persistence probabilites
offDiag = 1:(numClusters^2); offDiag <- offDiag[-onDiag] # construct indices to select only transition probabilites

### load control energies ###

E.wholebrain <- readMat(paste(masterdir,'analyses/control_energy/Tsweep/c',c,'_k',numClusters,
    '/GroupAverageTransitionEnergies_k',numClusters,'.mat',sep=''))
#T.idx <- which.min(abs(E.wholebrain$T.rng-T)) # get index of T value to use based on input
T.idx <- which.min(cor(E.wholebrain$E.full[offDiag,],restTP[offDiag],method=method))
T.opt.WB <- signif(E.wholebrain$T.rng[T.idx],2)
E.T.WB <- E.wholebrain$E.full[,T.idx]
E.wholebrain.BCTNull <- E.wholebrain$E.BCTnull[,,T.idx]

### load sps group average energies ###
E.wholebrain <- b$E.full
#T.idx <- which.min(abs(E.wholebrain$T.rng-T)) # get index of T value to use based on input
T.opt.WB <- 5
E.T.WB <- E.wholebrain

### load sps subject specific energies ###
E.wholebrain.LSD <- c$E.full[1:15,]
E.wholebrain.LSD <- matrix(colMeans(E.wholebrain.LSD),1,16)

E.wholebrain.PL <- c$E.full[16:30,]
E.wholebrain.PL <- matrix(colMeans(E.wholebrain.PL),1,16)


E.weighted.Null <-b$NullTransitionEnergy
E.weighted.Null <- matrix(colMeans(E.weighted.Null),1,16)

# this came out to 5.5 .. even though 5.5 gives max correlation at rest
# can stil use 5 for null testing, conservatively as it's not optimal
T.opt.WB <- 5


E.wholebrain.DLWNull <- readMat(paste(masterdir,'analyses/control_energy/DLWNullsTransitionEnergy_FA_k',
	numClusters,'c',c,'T',T.opt.WB,'.mat',sep=''))$DLWNullTransitionEnergy

### compute p-values ###

# whole brain control #
r.RealBrain.rest <- cor.test(E.T.WB[offDiag],matrix(t(restTP),1,16)[offDiag],method=method) #LSD: off diags / group avg E # correlation between real brain network energy and rest trans probs
r.RealBrain.rest <- cor.test(E.T.WB,matrix(t(restTP),1,16),method=method) #LSD: whole mat / group avg E
r.RealBrain.rest <- cor.test(E.wholebrain.LSD,matrix(t(restTP),1,16),method=method) #LSD: whole mat / subj specific E
r.RealBrain.rest <- cor.test(E.wholebrain.LSD[offDiag],matrix(t(restTP),1,16)[offDiag],method=method) #LSD: off diag / subj specific E

r.DLWNull.rest <- cor.test(E.weighted.Null[offDiag],restTP[offDiag],method=method) # correlation between DLW null energy and rest trans probs
r.DLWNull.rest <- cor.test(E.weighted.Null,matrix(restTP,1,16),method=method) 
p.DLWNull.rest <- pval.label.np(pvals=mean(r.RealBrain.rest >= r.DLWNull.rest),n=length(r.DLWNull.rest))

r.BCTNull.rest <- cor(E.wholebrain.BCTNull[offDiag,],restTP[offDiag],method=method) # correlation between BCT null energy and rest trans probs
p.BCTNull.rest <- pval.label.np(pvals=mean(r.RealBrain.rest >= r.BCTNull.rest),n=length(r.BCTNull.rest))

r.RealBrain.nback <- cor.test(E.T.WB[offDiag],matrix(t(nbackTP),1,16)[offDiag],method=method) #PL offdiags / grp avg E # correlation between real brain network energy and nback trans probs
r.RealBrain.nback <- cor.test(E.T.WB,matrix(t(nbackTP),1,16),method=method)  #PL: whole mat / grp avg E
r.RealBrain.nback <- cor.test(E.wholebrain.PL,matrix(t(nbackTP),1,16),method=method)  #PL: whole mat / subj specific E
r.RealBrain.nback <- cor.test(E.wholebrain.PL[offDiag],matrix(t(nbackTP),1,16)[offDiag],method=method)  #PL: off diag mat / subj specific E

r.DLWNull.nback <- cor(t(E.wholebrain.DLWNull[,offDiag]),nbackTP[offDiag],method=method) # correlation between DLW null energy and nback trans probs
p.DLWNull.nback <- pval.label.np(pvals=mean(r.RealBrain.nback >= r.DLWNull.nback),n=length(r.DLWNull.nback))

r.BCTNull.nback <- cor(E.wholebrain.BCTNull[offDiag,],nbackTP[offDiag],method=method) # correlation between BCT null energy and nback trans probs
p.BCTNull.nback <- pval.label.np(pvals=mean(r.RealBrain.nback >= r.BCTNull.nback),n=length(r.BCTNull.nback))

wholebrain.pvals.rest <- data.frame(test=c('DLW','BCT'),plabel=c(p.DLWNull.rest,p.BCTNull.rest))
wholebrain.pvals.nback <- data.frame(test=c('DLW','BCT'),plabel=c(p.DLWNull.nback,p.BCTNull.nback))

### plot energy vs. trans prob ###
# scatter.ETP() is in code/plottingfxns/plottingfxns.R

# convert to rank values so linear fit makes sense 
E.T.WB <- rank(E.T.WB[offDiag])
restTP <- rank(restTP[offDiag])
nbackTP <- rank(nbackTP[offDiag])
E.wholebrain.LSD <-rank(E.wholebrain.LSD[offDiag])
E.wholebrain.PL <- rank(E.wholebrain.PL[offDiag])

#or keep diagonals
E.T.WB <- rank(E.T.WB)
restTP <- rank(matrix(t(restTP),1,16))
nbackTP <- rank(matrix(t(nbackTP),1,16))
E.wholebrain.LSD <- rank(E.wholebrain.LSD)
E.wholebrain.PL <- rank(E.wholebrain.PL)

init <- rep(1:numClusters,each=numClusters-1) # starting cluster for each off-diagonal transition
term <- rep(1:numClusters,numClusters)[offDiag] # ending cluster for each off-diagonal transition

# note: rest looks like there are actualy several upward trends that on average form a downward trend
# I checked by coloring the points with init and term above, which revealed that
# these lines do not just reflect rows or columns of the transition matrix
p.rest <- scatter.ETP(E.T.WB,restTP,RNcolors[1],ttl='LSD',method=method) #,plabel.df=wholebrain.pvals.rest,method=method)
p.rest <- scatter.ETP(E.wholebrain.LSD,restTP,RNcolors[1],ttl='LSD',method=method) #,plabel.df=wholebrain.pvals.rest,method=method)
save(p.rest,E.T.WB,restTP,wholebrain.pvals.rest,file=paste0(savedir,'Fig5c__WholeBrainControlRest.RData'))
ggsave(plot = p.rest,filename = paste(savedir,'WholeBrainEnergyVsRestTP_c',c,'T',T.opt.WB,method,'_k',numClusters,'.pdf',sep =''),
	units = 'in',height = 1.5,width = 1.5)
p.nback <- scatter.ETP(E.T.WB,nbackTP,RNcolors[2],ttl='PL',method=method) #,plabel.df=wholebrain.pvals.nback,method=method)
p.nback <- scatter.ETP(E.wholebrain.PL,nbackTP,RNcolors[2],ttl='PL',method=method) 
save(p.nback,E.T.WB,nbackTP,wholebrain.pvals.nback,file=paste0(savedir,'Fig5c__WholeBrainControl2Back.RData'))
ggsave(plot = p.nback,filename = paste(savedir,'WholeBrainEnergyVs2backTP_c',c,'T',T.opt.WB,method,'_k',numClusters,'.pdf',sep =''),
	units = 'in',height = 1.5,width = 1.5)

p.rest
p.nback

p.caps <- scatter.ETP(limdiff,caps,RNcolors[1],ttl='LSD',method=method) #,plabel.df=wholebrain.pvals.rest,method=method)

p.caps
