####read in all data
setwd("/Users/Daphne/Documents/Yuan/Deconvolution/Deconvolution_Github/Simulation_Schelker")

#load packages and functions
source("../Deconvolution_Functions.R")

#set colors for plotting
cols<-c("grey","navy","yellow","darkgreen","green","lightgreen","purple","blue","maroon","red")

#true cell identities
cellnames<-read.table("data/cellnames.csv",sep=",",col.names=FALSE)
cellnames<-cellnames[,1]
cellname<-c(as.character(cellnames[1:2]),as.character(cellnames[6:13])) #collapse T cells into one group

#genes to include in signature, based on lm22 signature matrix from Newman et al.
lm22genes<-read.table("data/selected_genes_lm22.csv",sep=",")
lm22genes<-as.list(lm22genes)$V1

#cell type classification for pbmc data using clustering results from Schelker et al.
celltype<-read.table("data/celltype.csv",sep=",",col.names=FALSE)
celltype<-celltype[,1]

#all single-cell data
scdata<-read.table("data/singlecelldata.csv",sep=",",header=TRUE)

#consider only lm22 signature matrix genes
scdata_reduced<-scdata[lm22genes,]
data<-scdata_reduced

#to make signature, use only lm22 genes, all cells
#deconvolve each of the bulk data sets, created by adding up ALL patient samples (one bulk set for each patient)

#sample_id is the single-cell data set name to which the cell pertains
sample_id<-read.table("data/sampleid.csv",sep=",",col.names=FALSE)
sample_id<-sample_id[,1]

#create bulk data sets- should be 27 in total

#f_true contains the true cell type proportions for each data set
f_true = matrix(nrow=12,ncol=length(unique(sample_id)))
#m contains simulated bulk data for each of 27 single-cell data sets
m = matrix(nrow=length(lm22genes),ncol=length(unique(sample_id)))


for (j in 1:length(unique(sample_id))){
  ix<-which(sample_id==unique(sample_id)[j])
  for (i in 1:length(unique(celltype))){
    iy<-which(celltype==unique(celltype)[i])
    f_true[i,j]<-length(intersect(iy,ix))/length(ix)
  }
  m[,j]<-rowMeans(scdata_reduced[,ix])
}

#reorder cell types
reorder<-order(unique(celltype))
f_true_ordered<-f_true[reorder,]

#collapse T cell groups
F_true_ordered<-matrix(nrow=10,ncol=27)
F_true_ordered[1,]<-f_true_ordered[1,]
F_true_ordered[2,]<-colSums(f_true_ordered[2:4,])
F_true_ordered[3:10,]<-f_true_ordered[5:12,]

#create signature matrix from all samples
B<-matrix(nrow=dim(scdata_reduced)[1],ncol=length(unique(celltype)))
for (i in 1:length(unique(celltype))){
  B[,i]<-rowMeans(scdata_reduced[,which(celltype==unique(celltype)[i])])
}
#remove the profile for unknown cells (celltype 0)
unknownindex<-which(unique(celltype)==0)
B[,unknownindex]<-rep(NA,dim(B)[1])
#reorder
B<-B[,reorder]
############

####perform deconvolution
solutionsOLS<-NULL
solutionsDWLSIterative<-NULL
solutionsSVR<-NULL
for (i in 1:27){
  m_j<-m[,i]
  B_j<-B[,-1]
  #
  sol<-solveOLS(B_j,m_j)
  sol<-round(sol,5)
  sol_full<-c(0,sol) #no unknown type
  #
  solDWLSIterative<-solveDampenedWLS(B_j,m_j)
  solDWLSIterative<-round(solDWLSIterative,5)
  solDWLSIterative_full<-c(0,solDWLSIterative) #no unknown type
  #
  solSVR<-solveSVR(B_j,m_j)
  solSVR<-round(solSVR,5)
  solSVR_full<-c(0,solSVR)
  #
  solutionsOLS<-cbind(solutionsOLS,sol_full)
  solutionsDWLSIterative<-cbind(solutionsDWLSIterative,solDWLSIterative_full)
  solutionsSVR<-cbind(solutionsSVR,solSVR_full)
}

#order solutions and collapse T cell groups
SolOLS_ordered<-matrix(nrow=10,ncol=27)
SolOLS_ordered[1,]<-solutionsOLS[1,]
SolOLS_ordered[2,]<-colSums(solutionsOLS[2:4,])
SolOLS_ordered[3:10,]<-solutionsOLS[5:12,]

SolDWLS_ordered<-matrix(nrow=10,ncol=27)
SolDWLS_ordered[1,]<-solutionsDWLSIterative[1,]
SolDWLS_ordered[2,]<-colSums(solutionsDWLSIterative[2:4,])
SolDWLS_ordered[3:10,]<-solutionsDWLSIterative[5:12,]

SolSVR_ordered<-matrix(nrow=10,ncol=27)
SolSVR_ordered[1,]<-solutionsSVR[1,]
SolSVR_ordered[2,]<-colSums(solutionsSVR[2:4,])
SolSVR_ordered[3:10,]<-solutionsSVR[5:12,]

#save deconvolution results
save(SolOLS_ordered,file="results/SolOLS_ordered.RData")
save(SolDWLS_ordered,file="results/SolDWLS_ordered.RData")
save(SolSVR_ordered,file="results/SolSVR_ordered.RData")

#calculate correlations between truth and estimates

correlationsOLS<-c()
correlationsDWLS<-c()
correlationsSVR<-c()
for(j in 1:10){
  correlationsOLS<-c(correlationsOLS,cor(F_true_ordered[j,],SolOLS_ordered[j,]))
  correlationsDWLS<-c(correlationsDWLS,cor(F_true_ordered[j,],SolDWLS_ordered[j,]))
  correlationsSVR<-c(correlationsSVR,cor(F_true_ordered[j,],SolSVR_ordered[j,]))
}

correlationsOLS
correlationsDWLS
correlationsSVR

#overall correlations
correlationsSVRoverall<-cor(as.vector(F_true_ordered),as.vector(SolSVR_ordered))
correlationsOLSoverall<-cor(as.vector(F_true_ordered),as.vector(SolOLS_ordered))
correlationsDWLSoverall<-cor(as.vector(F_true_ordered),as.vector(SolDWLS_ordered))

#plot results of estimated vs. true

types<-c("SVR","OLS","DWLS")

for (type in types){
  pdf(paste("results/Correlations_",type,".pdf",sep=""),width=10)
  op <- par(mfrow = c(3,4),
            oma = c(5,4,0,0) + 1,
            mar = c(0,0,2,2) + 1)
  estimates<-eval(parse(text=paste("Sol",type,"_ordered",sep="")))
  plot(as.vector(F_true_ordered),as.vector(estimates),main="Overall",ylim=c(0,1),xlim=c(0,1),col="black")
  abline(0,1)
  text(.75,0.25,labels=paste("cor=",round(eval(parse(text=paste("correlations",type,"overall",sep=""))),2),sep=""))
  for(i in 1:10){
    plot(F_true_ordered[i,],estimates[i,],main=cellname[i],ylim=c(0,1),xlim=c(0,1),col=cols[i])
    abline(0,1)
    text(.75,0.25,labels=paste("cor=",round(eval(parse(text=paste("correlations",type,sep="")))[i],2),sep=""))
  }
  dev.off()
}


#calculate modified relative error
pseudo<-0.005
relativeErrorDWLS<-NULL
for(j in 1:10){
  relativeErrorDWLS<-cbind(relativeErrorDWLS,abs((F_true_ordered[j,]+pseudo)-(SolDWLS_ordered[j,]+pseudo))*100/(((F_true_ordered[j,]+pseudo)+(SolDWLS_ordered[j,]+pseudo))/2))
}
relativeErrorDWLS[which(is.na(relativeErrorDWLS))]<-0

relativeErrorOLS<-NULL
for(j in 1:10){
  relativeErrorOLS<-cbind(relativeErrorOLS,abs((F_true_ordered[j,]+pseudo)-(SolOLS_ordered[j,]+pseudo))*100/((F_true_ordered[j,]+SolOLS_ordered[j,]+2*pseudo)/2))
}
relativeErrorOLS[which(is.na(relativeErrorOLS))]<-0

relativeErrorSVR<-NULL
for(j in 1:10){
  relativeErrorSVR<-cbind(relativeErrorSVR,abs(F_true_ordered[j,]+pseudo-(SolSVR_ordered[j,]+pseudo))*100/((F_true_ordered[j,]+SolSVR_ordered[j,]+2*pseudo)/2))
}
relativeErrorSVR[which(is.na(relativeErrorSVR))]<-0

#plot average modified relative error for each cell type for each data set
pdf("results/Relative_Error_vs_Cell_Type_Proportion.pdf")
par(mfrow=c(1,1))
plot(rowMeans(F_true_ordered),colMeans(relativeErrorDWLS),ylab="Relative percent error",xlab="Average proportion of cell type",col="black",pch=21,bg="lightblue")
points(rowMeans(F_true_ordered),colMeans(relativeErrorOLS),col="black",pch=21,bg="pink")
points(rowMeans(F_true_ordered),colMeans(relativeErrorSVR),col="black",pch=21,bg="plum3")
abline(lm(colMeans(relativeErrorDWLS)~rowMeans(F_true_ordered)),col="lightblue",lwd=2)
abline(lm(colMeans(relativeErrorOLS)~rowMeans(F_true_ordered)),col="pink",lwd=2)
abline(lm(colMeans(relativeErrorSVR)~rowMeans(F_true_ordered)),col="plum3",lwd=2)
legend("topright",c("DWLS","OLS","SVR"),col=c("lightblue","pink","plum3"),pch="-",lwd=2)
dev.off()

#calculate overall modified relative error
relativeErrorSVRoverall<-abs(as.vector(F_true_ordered+pseudo)-as.vector(SolSVR_ordered+pseudo))*100/((as.vector(F_true_ordered)+as.vector(SolSVR_ordered)+2*pseudo)/2)
relativeErrorSVRoverall[which(is.na(relativeErrorSVRoverall))]<-0
mean(relativeErrorSVRoverall)

relativeErrorDWLSoverall<-abs(as.vector(F_true_ordered+pseudo)-as.vector(SolDWLS_ordered+pseudo))*100/((as.vector(F_true_ordered)+as.vector(SolDWLS_ordered)+2*pseudo)/2)
relativeErrorDWLSoverall[which(is.na(relativeErrorDWLSoverall))]<-0
mean(relativeErrorDWLSoverall)

relativeErrorOLSoverall<-abs(as.vector(F_true_ordered+pseudo)-as.vector(SolOLS_ordered+pseudo))*100/((as.vector(F_true_ordered)+as.vector(SolOLS_ordered)+2*pseudo)/2)
relativeErrorOLSoverall[which(is.na(relativeErrorOLSoverall))]<-0
mean(relativeErrorOLSoverall)

#calculate absolute error
absoluteErrorSVR<-NULL
for(j in 1:10){
  absoluteErrorSVR<-cbind(absoluteErrorSVR,abs((F_true_ordered[j,])-(SolSVR_ordered[j,])))
}
AverageAbsoluteErrorSVR<-colMeans(absoluteErrorSVR)

absoluteErrorSVRoverall<-abs(as.vector(F_true_ordered)-as.vector(SolSVR_ordered))
mean(absoluteErrorSVRoverall)

absoluteErrorOLS<-NULL
for(j in 1:10){
  absoluteErrorOLS<-cbind(absoluteErrorOLS,abs((F_true_ordered[j,])-(SolOLS_ordered[j,])))
}
AverageAbsoluteErrorOLS<-colMeans(absoluteErrorOLS)

absoluteErrorOLSoverall<-abs(as.vector(F_true_ordered)-as.vector(SolOLS_ordered))
mean(absoluteErrorOLSoverall)

absoluteErrorDWLS<-NULL
for(j in 1:10){
  absoluteErrorDWLS<-cbind(absoluteErrorDWLS,abs((F_true_ordered[j,])-(SolDWLS_ordered[j,])))
}
AverageAbsoluteErrorDWLS<-colMeans(absoluteErrorDWLS)

absoluteErrorDWLSoverall<-abs(as.vector(F_true_ordered)-as.vector(SolDWLS_ordered))
mean(absoluteErrorDWLSoverall)

