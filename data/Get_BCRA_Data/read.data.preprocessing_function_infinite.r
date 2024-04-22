# 12/11/19
library(data.table)
library(REMP)
#study <- "GBMLGG"
data.processing <- function(study){

meth.beta <- fread(paste("HumanMethylation450_",study, sep=""),header=TRUE, check.names=FALSE)
sj.names <- read.table(paste("HumanMethylation450_",study, sep=""),header=FALSE, row.names=1,check.names=FALSE,nrows=1,stringsAsFactors=FALSE)[1,]

meth.beta <- data.frame(meth.beta)
rownames(meth.beta) <-meth.beta[,1]

meth.beta <- meth.beta[,-1]
colnames(meth.beta) <- sj.names
#meth.beta <- read.table("HumanMethylation450_HNSC",header=TRUE, row.names=1,check.names=FALSE)
 # 485577    580

meth.m <- log2(meth.beta/(1-meth.beta))

meth.m[meth.beta == 1] <- 10 # remove infite case
meth.m[meth.beta == 0] <- -10
#save(meth.m, file="meth.m.Rdata")
#any(meth.beta == 1, na.rm = TRUE)
# any(is.infinite(meth.m))
# #load("meth.m.Rdata")
# length(is.infinite(meth.m)) # 440418339
# sum(is.infinite(meth.m)) #25508
# meth.m.org <- meth.m
# meth.m[is.infinite(meth.m.org)] <- NA
#meth.m["cg18106189",]
#any.infinite<- apply(meth.m, 2, function(x) any(is.infinite(x)))
# any.na<- apply(meth.m, 1, function(x) any(is.na(x)))
# 
# sum(any.na)/nrow(meth.m) # 0.2276323 missing or infinite value
# 
n <- ncol(meth.m)
percent.na<- apply(meth.m, 1, function(x) sum(is.na(x))/n)
# hist(percent.na[percent.na!=0])
# 
# with.na<- percent.na[percent.na!=0]
# sum(with.na <=0.3)
# sum(with.na <=0.2)
impute.ind <- percent.na!=0&percent.na<=0.2
if(sum(impute.ind) > 30000) {
  temp <- grooMethy(meth.m[impute.ind,]) # impute.ind 
  imputed<-getM(temp)
  
  
  # hist(imputed[which.na2])
  # hist(imputed[!which.na2])
  # 
  # any(is.na(imputed))
  
  meth.m.imputed <- meth.m
  meth.m.imputed[impute.ind,] <- imputed
} else {
  impute.ind2 <- percent.na<=0.2
  temp <- grooMethy(meth.m[impute.ind2,]) # impute.ind 
  imputed<-getM(temp)
  
  
  # hist(imputed[which.na2])
  # hist(imputed[!which.na2])
  # 
  # any(is.na(imputed))
  
  meth.m.imputed <- meth.m
  meth.m.imputed[impute.ind2,] <- imputed
}
# which.na<-apply(meth.m[impute.ind,], 1, function(x) which(is.na(x)))
# which.na2<-is.na(meth.m[impute.ind,])





any.na<- apply(meth.m.imputed, 1, function(x) any(is.na(x)))

sum(any.na)/nrow(meth.m.imputed) #18.5% missing

#save(meth.m.imputed, file="meth.m.imputed.Rdata")

meth.m.imputed.no.na <- meth.m.imputed[!any.na,] # no missing value only
#save(meth.m.imputed.no.na, file="meth.m.imputed.no.na.Rdata")

##################################################
# gene expression
gene.exp0<- read.table(paste("HiSeqV2_",study, sep=""),header=TRUE, row.names=1,check.names=FALSE)

gene.exp.mean <- rowMeans(gene.exp0)
sum(gene.exp.mean<=1) #3305
high.gene.exp <- gene.exp0[gene.exp.mean>1,]
meth0<-meth.m.imputed.no.na[,colnames(meth.m.imputed.no.na) %in% colnames(high.gene.exp)]
gene.exp <- high.gene.exp[,colnames(meth0)] #856 individual
dim(gene.exp) #17225   542
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probe <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19, 
                       what = "Locations")
probe.range <- GRanges(seqnames = probe$chr, ranges = IRanges(probe$pos, 
                                                              width = 1, names = rownames(probe)), strand = probe$strand, 
                       name = rownames(probe))

id <- match(rownames(meth0), names(probe.range))
sort.non.na.id <- sort(id[!is.na(id)])
meth.range <- probe.range[sort.non.na.id]

meth<- meth0[names(meth.range),]

data <- list(meth=meth, gene.exp=gene.exp, meth.range=meth.range)
#save(data, file="data.Rdata" )

#synonym search

promoter<-read.delim("hg19_promoter.txt")
promoter.range <- GRanges(seqnames = promoter$chrID, ranges = IRanges(start=promoter$start, end=promoter$end), strand = promoter$strand, gene.name =promoter$gene.name)


source("synonym.search.r")
syn.name<- synonym.search(promoter.range$gene.name,rownames(data$gene.exp))

org.name<-rownames(data$gene.exp)


data$gene.exp<- as.matrix(data$gene.exp)
rownames(data$gene.exp) <- syn.name
#save(data,file=paste(study,".data.Rdata", sep=""))
######################################
# promoter.per.probe
promoter<-read.delim("hg19_promoter.txt")
promoter.range <- GRanges(seqnames = promoter$chrID, ranges = IRanges(start=promoter$start, end=promoter$end), strand = promoter$strand, gene.name =promoter$gene.name)

#load("HNSC.data2.Rdata")

probe.range <- data$meth.range
ids <- unique(subjectHits(findOverlaps(probe.range, promoter.range)))
sids <- sort(ids)
library(parallel)
probe.per.promoter<-mclapply(1:length(sids), function(curi){
  num<-sids[curi]
  probe.range[queryHits(findOverlaps(probe.range, 
                                     promoter.range[num]))]
},mc.cores=20)
names(probe.per.promoter) <- promoter$gene.name[sids]
save(probe.per.promoter, file=paste("probe.per.promoter_",study,".Rdata", sep=""))
}
