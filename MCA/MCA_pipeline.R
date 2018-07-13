setwd("/Users/Daphne/Documents/Yuan/Deconvolution/Deconvolution_Github/MCA")

#load packages and functions
source("../Deconvolution_Functions.R")

#load single-cell data
X=load("data/MCA_33tissue_reduced_dge_61K.Rdata")
#load single-cell labels
load("data/labels.RData")

#create signature matrix
Sig<-buildSignatureMatrixUsingSeurat(merged_filter,labels,"results") #will take a while, perform in cluster, not PC


#load signature matrix
load("results/Sig.RData") #create signature matrix in cluster, not PC

#load bulk samples
bulkIds<-c(1:4,9:12)
for (id in bulkIds){
  load(paste("data/bulkData",id,".RData",sep=""))
}
#tissue identities of bulk samples
bulkIdTissues<-c("Lung","Liver","Lung","Liver","SmallIntestine","Kidney","SmallIntestine","Kidney")

#apply deconvolution methods to each bulk data set
allCounts_DWLS<-NULL
allCounts_OLS<-NULL
allCounts_SVR<-NULL

for (id in bulkIds){
  bulkData<-eval(parse(text=paste("bulkData",id,sep="")))
  trimmed<-trimData(Sig,bulkData)
  S<-trimmed$sig
  B<-trimmed$bulk
  solOLS<-solveOLS(S,B)
  solDWLS<-solveDampenedWLS(S,B)
  solSVR<-solveSVR(S,B)
  allCounts_DWLS<-cbind(allCounts_DWLS,solDWLS)
  allCounts_OLS<-cbind(allCounts_OLS,solOLS)
  allCounts_SVR<-cbind(allCounts_SVR,solSVR)
}

#round and save solutions

allCounts_SVR<-round(allCounts_SVR,5)
allCounts_DWLS<-round(allCounts_DWLS,5)
allCounts_OLS<-round(allCounts_OLS,5)

save(allCounts_SVR,file="results/allCounts_SVR.RData")
save(allCounts_DWLS,file="results/allCounts_DWLS.RData")
save(allCounts_OLS,file="results/allCounts_OLS.RData")

#visualize results


#for kidney only

#correlation plots
load("data/goldStandard_Kidney.RData")

correlationOLS<-cor(allCounts_OLS[,6],goldStandard_Kidney)
correlationDWLS<-cor(allCounts_DWLS[,6],goldStandard_Kidney)
correlationSVR<-cor(allCounts_SVR[,6],goldStandard_Kidney)

colours<-c("pink","lightblue","plum3","grey")
pdf("results/Correlations_kidney.pdf",width=10,height=3.5)
par(mfrow=c(1,3))
plot(as.vector(goldStandard_Kidney),allCounts_OLS[,6],main="QP",ylab="inferred",xlab="true",xlim=c(0,1),ylim=c(0,1),col=colours[1])
usr <- par( "usr" )
abline(0,1,col=colours[1])
text(usr[ 2 ], usr[ 4 ],adj = c( 1.5, 1.5 ),labels=paste("Cor=",round(mean(correlationOLS),3),sep=""))
plot(as.vector(goldStandard_Kidney),allCounts_DWLS[,6],main="DWLS",ylab="inferred",xlab="true",xlim=c(0,0.6),ylim=c(0,0.6),col=colours[2])
usr <- par( "usr" )
text(usr[ 2 ], usr[ 4 ],adj = c( 1.5, 1.5 ),labels=paste("Cor=",round(mean(correlationDWLS),3),sep=""))
abline(0,1,col=colours[2])
plot(as.vector(goldStandard_Kidney),allCounts_SVR[,6],main="SVR",ylab="inferred",xlab="true",xlim=c(0,0.6),ylim=c(0,0.6),col=colours[3])
usr <- par( "usr" )
text(usr[ 2 ], usr[ 4 ],adj = c( 1.5, 1.5 ),labels=paste("Cor=",round(mean(correlationSVR),3),sep=""))
abline(0,1,col=colours[3])
dev.off()

#plot results as a heatmap
counts<-cbind(allCounts_OLS[,6],allCounts_DWLS[,6],allCounts_SVR[,6],goldStandard_Kidney)

colnames(counts)<-c("QP","DWLS","SVR","Truth")
library(RColorBrewer)
pdf("results/Heatmap_kidney.pdf")
heatmap.2(log(counts+.01),scale="none",col=brewer.pal(9,"Blues"),Colv=NA,Rowv=NA,density.info="none", trace="none",labRow = FALSE,margins = c(6, 12))
dev.off()


#for all
#calculate correlation, sensitivity, specificity

#calculate correlation
corOLS<-c()
corDWLS<-c()
corSVR<-c()
for (i in 1:length(bulkIds)){
  #for each bulk data set, calculate correlation
  #first, find gold standard for corresponding tissue
  tissue<-bulkIdTissues[i]
  load(paste("data/goldStandard_",tissue,".RData",sep=""))
  goldStandard<-eval(parse(text=paste("goldStandard_",tissue,sep="")))
  corOLS<-c(corOLS,cor(allCounts_OLS[,i],goldStandard))
  corDWLS<-c(corDWLS,cor(allCounts_DWLS[,i],goldStandard))
  corSVR<-c(corSVR,cor(allCounts_SVR[,i],goldStandard))
}

#boxplot of correlation
pdf("results/Boxplot_correlations.pdf")
boxplot(corOLS,corDWLS,corSVR,col=colours,names=c("QP","DWLS","SVR"))
dev.off()

#calculate specificity and sensitivity
p<-0.02
specOLS<-c()
specDWLS<-c()
specSVR<-c()
sensOLS<-c()
sensDWLS<-c()
sensSVR<-c()
for (i in 1:length(bulkIds)){
  #for each bulk data set, calculate correlation
  #first, find gold standard for corresponding tissue
  tissue<-bulkIdTissues[i]
  load(paste("data/goldStandard_",tissue,".RData",sep=""))
  goldStandard<-eval(parse(text=paste("goldStandard_",tissue,sep="")))
  correctTissues<-names(goldStandard)[which(goldStandard>p)]
  incorrectTissues<-names(goldStandard)[which(goldStandard<=p)]
  specDWLS<-c(specDWLS,length(intersect(names(allCounts_DWLS[,i])[which(allCounts_DWLS[,i]<=p)],incorrectTissues))/length(incorrectTissues))
  specOLS<-c(specOLS,length(intersect(names(allCounts_OLS[,i])[which(allCounts_OLS[,i]<=p)],incorrectTissues))/length(incorrectTissues))
  specSVR<-c(specSVR,length(intersect(names(allCounts_SVR[,i])[which(allCounts_SVR[,i]<=p)],incorrectTissues))/length(incorrectTissues))
  sensDWLS<-c(sensDWLS,length(intersect(names(allCounts_DWLS[,i])[which(allCounts_DWLS[,i]>p)],correctTissues))/length(correctTissues))
  sensOLS<-c(sensOLS,length(intersect(names(allCounts_OLS[,i])[which(allCounts_OLS[,i]>p)],correctTissues))/length(correctTissues))
  sensSVR<-c(sensSVR,length(intersect(names(allCounts_SVR[,i])[which(allCounts_SVR[,i]>p)],correctTissues))/length(correctTissues))               
}

#boxplot of specificity
pdf("results/Boxplot_specificity.pdf")
boxplot(specOLS,specDWLS,specSVR,col=colours,names=c("QP","DWLS","SVR"))
dev.off()

#boxplot of sensitivity
pdf("results/Boxplot_sensitivity.pdf")
boxplot(sensOLS,sensDWLS,sensSVR,col=colours,names=c("QP","DWLS","SVR"))
dev.off()
