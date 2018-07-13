setwd("/Users/Daphne/Documents/Yuan/Deconvolution/Deconvolution_Github/ISC")

#load packages and functions
source("../Deconvolution_Functions.R")

#######load data
#load bulk data
bulkData<-read.csv("data/GSE92377_FPKM_summary.csv")

#load single-cell data and labels
load("data/dataSC.RData")
load("data/trueLabels.RData")

labels<-trueLabels

#change to real labels
newcat<-c("NonCycISC","CycISC","TA","Ent","PreEnt","Goblet","Paneth","Tuft","EE")
for (i in 1:length(newcat)){
  labels[which(labels==(i-1))]<-newcat[i]
}


########deconvolution

Signature<-buildSignatureMatrixMAST(dataSC,labels,"results")

allCounts_DWLS<-NULL
allCounts_OLS<-NULL
allCounts_SVR<-NULL
for(j in 1:(dim(bulkData)[2]-1)){
  S<-Signature
  Bulk<-bulkData[,j+1]
  names(Bulk)<-bulkData[,1]
  Genes<-intersect(rownames(S),names(Bulk))
  B<-Bulk[Genes]
  S<-S[Genes,]
  solOLS<-solveOLS(S,B)
  solDWLS<-solveDampenedWLS(S,B)
  solSVR<-solveSVR(S,B)

  allCounts_DWLS<-cbind(allCounts_DWLS,solDWLS)
  allCounts_OLS<-cbind(allCounts_OLS,solOLS)
  allCounts_SVR<-cbind(allCounts_SVR,solSVR)
}

#save solutions

save(allCounts_SVR,file="results/allCounts_SVR.RData")
save(allCounts_DWLS,file="results/allCounts_DWLS.RData")
save(allCounts_OLS,file="results/allCounts_OLS.RData")


library(lattice)
pdf("results/Boxplot_DWLS.pdf",width=10)
boxplot(allCounts_DWLS[1,1:2],allCounts_DWLS[1,3:5],allCounts_DWLS[1,6:8],allCounts_DWLS[2,1:2],allCounts_DWLS[2,3:5],allCounts_DWLS[2,6:8],allCounts_DWLS[3,1:2],allCounts_DWLS[3,3:5],allCounts_DWLS[3,6:8],colSums(allCounts_DWLS[4:9,1:2]),colSums(allCounts_DWLS[4:9,3:5]),colSums(allCounts_DWLS[4:9,6:8]),names=c("Control","LOF","GOF","Control","LOF","GOF","Control","LOF","GOF","Control","LOF","GOF"),xlab="Condition",ylab="Proportion",ylim=c(0,1))
stripchart(c(allCounts_DWLS[1,1:2],allCounts_DWLS[1,3:5],allCounts_DWLS[1,6:8],allCounts_DWLS[2,1:2],allCounts_DWLS[2,3:5],allCounts_DWLS[2,6:8],allCounts_DWLS[3,1:2],allCounts_DWLS[3,3:5],allCounts_DWLS[3,6:8],colSums(allCounts_DWLS[4:9,1:2]),colSums(allCounts_DWLS[4:9,3:5]),colSums(allCounts_DWLS[4:9,6:8]))~c(1,1,2,2,2,3,3,3,4,4,5,5,5,6,6,6,7,7,8,8,8,9,9,9,10,10,11,11,11,12,12,12), vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, lwd=2,col = 'blue')
text(2,1,"Cycling ISC")
text(5,1,"Noncycling ISC")
text(8,1,"Transit Amplifying")
text(11,1,"Differentiated")
dev.off()

pdf("results/Boxplot_OLS.pdf",width=10)
boxplot(allCounts_OLS[1,1:2],allCounts_OLS[1,3:5],allCounts_OLS[1,6:8],allCounts_OLS[2,1:2],allCounts_OLS[2,3:5],allCounts_OLS[2,6:8],allCounts_OLS[3,1:2],allCounts_OLS[3,3:5],allCounts_OLS[3,6:8],colSums(allCounts_OLS[4:9,1:2]),colSums(allCounts_OLS[4:9,3:5]),colSums(allCounts_OLS[4:9,6:8]),names=c("Control","LOF","GOF","Control","LOF","GOF","Control","LOF","GOF","Control","LOF","GOF"),xlab="Condition",ylab="Proportion",ylim=c(0,1))
stripchart(c(allCounts_OLS[1,1:2],allCounts_OLS[1,3:5],allCounts_OLS[1,6:8],allCounts_OLS[2,1:2],allCounts_OLS[2,3:5],allCounts_OLS[2,6:8],allCounts_OLS[3,1:2],allCounts_OLS[3,3:5],allCounts_OLS[3,6:8],colSums(allCounts_OLS[4:9,1:2]),colSums(allCounts_OLS[4:9,3:5]),colSums(allCounts_OLS[4:9,6:8]))~c(1,1,2,2,2,3,3,3,4,4,5,5,5,6,6,6,7,7,8,8,8,9,9,9,10,10,11,11,11,12,12,12), vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, lwd=2,col = 'blue')
text(2,1,"Cycling ISC")
text(5,1,"Noncycling ISC")
text(8,1,"Transit Amplifying")
text(11,1,"Differentiated")
dev.off()

pdf("results/Boxplot_SVR.pdf",width=10)
boxplot(allCounts_SVR[1,1:2],allCounts_SVR[1,3:5],allCounts_SVR[1,6:8],allCounts_SVR[2,1:2],allCounts_SVR[2,3:5],allCounts_SVR[2,6:8],allCounts_SVR[3,1:2],allCounts_SVR[3,3:5],allCounts_SVR[3,6:8],colSums(allCounts_SVR[4:9,1:2]),colSums(allCounts_SVR[4:9,3:5]),colSums(allCounts_SVR[4:9,6:8]),names=c("Control","LOF","GOF","Control","LOF","GOF","Control","LOF","GOF","Control","LOF","GOF"),xlab="Condition",ylab="Proportion",ylim=c(0,1))
stripchart(c(allCounts_SVR[1,1:2],allCounts_SVR[1,3:5],allCounts_SVR[1,6:8],allCounts_SVR[2,1:2],allCounts_SVR[2,3:5],allCounts_SVR[2,6:8],allCounts_SVR[3,1:2],allCounts_SVR[3,3:5],allCounts_SVR[3,6:8],colSums(allCounts_SVR[4:9,1:2]),colSums(allCounts_SVR[4:9,3:5]),colSums(allCounts_SVR[4:9,6:8]))~c(1,1,2,2,2,3,3,3,4,4,5,5,5,6,6,6,7,7,8,8,8,9,9,9,10,10,11,11,11,12,12,12), vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, lwd=2,col = 'blue')
text(2,1,"Cycling ISC")
text(5,1,"Noncycling ISC")
text(8,1,"Transit Amplifying")
text(11,1,"Differentiated")
