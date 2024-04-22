
library(GenomicRanges)
library(glmnet)
#source("divide.data.r")
study <- "BRCA"
library(R.utils)
load(paste(study,".data.Rdata", sep=""))
load(paste("probe.per.promoter_",study,".Rdata", sep=""))
sample.id<-colnames(data$meth)
patient.id <- substr(sample.id, 1, 12)

dup.tf <- duplicated(patient.id)|duplicated(patient.id, fromLast=TRUE)
dup.id <- patient.id[dup.tf]

dup.unique.id <- unique(dup.id)
#test<-foldids[patient.id%in%dup.unique.id ,]
#test[order(test$patient.id),]


promoter<-read.delim("hg19_promoter.txt")
promoter.range <- GRanges(seqnames = promoter$chrID, ranges = IRanges(start=promoter$start, end=promoter$end), strand = promoter$strand, gene.name =promoter$gene.name)

gene.lists <-promoter.range$gene.name[promoter.range$gene.name %in% rownames(data$gene.exp)]

gene.list2 <- gene.lists[gene.lists %in% names(probe.per.promoter)]

n <- sum(ncol(data$gene.exp))
#curi<-2897##MKRN3
r2 <- c()
r2.name <- c()
#curi=2
nfolds<-5
lambda.rule <- "lambda.min"
enhancer.range=1e7
library(parallel)
seq.num <- 1: length(gene.list2)
sub.gene.list <- gene.list2[seq.num]
for(curi in c(1:(length(seq.num)))){
promoter.range1 <- promoter.range[promoter.range$gene.name == sub.gene.list[curi]]
candi.enhancer.range1<- GRanges(seqnames=seqnames(promoter.range1), 
                                IRanges(start=start(ranges(promoter.range1)) -enhancer.range,end=end(ranges(promoter.range1)) + enhancer.range),
                                ene.name=promoter.range1$gene.name)
probe.id <- unique(queryHits(findOverlaps(data$meth.range, candi.enhancer.range1)))
#probe.name.unique <- unique(c(promoter.probe.name, gene.body.probe.name))
y=t(data$gene.exp[rownames(data$gene.exp) ==  sub.gene.list[curi],,drop=FALSE])
x=t(data$meth[probe.id,])



p<-ncol(data)-1
# x1 = data[,-ncol(data)]
set.seed(1234)

p<-ncol(data)-1
# x1 = data[,-ncol(data)]
set.seed(1234)
sam1<- sample(rep(seq(nfolds), length = length(dup.unique.id)))

dup.unique.foldid<- data.frame(patient.id=dup.unique.id, foldid=sam1)
sam2<- sample(rep(nfolds:1, length = sum(!dup.tf)))

rest.foldid<- data.frame(patient.id=patient.id[!dup.tf], foldid=sam2)

mergy.foldid <- rbind(dup.unique.foldid, rest.foldid)

#foldids <- data.frame(sample.id=sample.id, mergy.foldid[match(patient.id, mergy.foldid$patient.id),])
foldid <- mergy.foldid[match(patient.id, mergy.foldid$patient.id),"foldid"]
x=apply(x,2,as.double)
ny <- ncol(y)
if(ny ==1) {
  y <- as.numeric(y)
  cv.fit <- cv.glmnet(x, y, keep=TRUE, alpha=0.5, foldid=foldid)
  #coef(cv.fit, s=lambda.rule)
  id<-which(cv.fit$lambda == cv.fit$lambda.min)
  if(var(cv.fit$fit.preval[,id])!=0) {
    R2 <- cor(cv.fit$fit.preval[,id],y)^2
  } else {
    R2 <- 0
  }
} else {
  
  sR2 <- rep(0, ny)
  for(i in 1:ny) {
    y.m <- as.numeric(y[,i])
    cv.fit <- cv.glmnet(x, y.m, keep=TRUE, alpha=0.5, foldid=foldid)
    id<-which(cv.fit$lambda == cv.fit$lambda.min)
    
    if(var(cv.fit$fit.preval[,id])!=0) {
      sR2[i] <- cor(cv.fit$fit.preval[,id],y.m)^2
    } else {
      sR2[i] <- 0
    }
    
  } 
  R2 <- max(sR2)
}
r2.name <-c(r2.name,sub.gene.list[curi])
r2=c(r2,R2)
}

bb=as.data.frame(r2)
bb$name=r2.name
library(data.table)
write.csv(bb,"../result_r2/geneExplorer_5fold.csv")

