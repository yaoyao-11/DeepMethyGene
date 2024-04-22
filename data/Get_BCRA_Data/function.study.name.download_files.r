# 01/23/21
# install GenomicRanges, glmnet, data.table, R.utils, and REMP
# input study abbreviation
# default range is +-10Mb from promoter regions (enhancer.range)

library(GenomicRanges)
library(glmnet)
#source("divide.data.r")
study <- "BRCA"

file <- paste("https://tcga.xenahubs.net/download/TCGA.",study,".sampleMap/HumanMethylation450.gz",sep="")
download.file(file, paste("HumanMethylation450_",study,".gz", sep=""))
library(R.utils)
gunzip(paste("HumanMethylation450_",study,".gz", sep=""))
file2 <- paste("https://tcga.xenahubs.net/download/TCGA.",study,".sampleMap/HiSeqV2.gz",sep="")

download.file(file2, paste("HiSeqV2_",study,".gz", sep=""))
gunzip(paste("HiSeqV2_",study,".gz", sep=""))

source("read.data.preprocessing_function_infinite.r")
data.processing(study)

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

write.csv(data$meth, "BRCA_data_meth.csv", row.names = FALSE)
promoter.range.df <- as.data.frame(promoter.range)
write.csv(promoter.range.df, "BRCA_data_meth_range.csv", row.names = FALSE)
n <- sum(ncol(data$gene.exp))
curi<-1
nfolds<-10
lambda.rule <- "lambda.min"
#save(gene.list2, file="genome.list.Rdata")
#save(gene.list2, file="gene.list2.Rdata")
#gene.list2 <-gene.list2[1:20]
enhancer.range=1e7
library(parallel)
sub.gene.list <- gene.list2[seq.num]
promoter.range1 <- promoter.range[promoter.range$gene.name == sub.gene.list[curi]]
candi.enhancer.range1<- GRanges(seqnames=seqnames(promoter.range1), IRanges(start=start(ranges(promoter.range1)) -enhancer.range,end=end(ranges(promoter.range1)) + enhancer.range),gene.name=promoter.range1$gene.name)
probe.id <- unique(queryHits(findOverlaps(data$meth.range, candi.enhancer.range1)))
    #probe.name.unique <- unique(c(promoter.probe.name, gene.body.probe.name))
need_gene_exp_13988=t(data$gene.exp[rownames(data$gene.exp) ==  sub.gene.list[curi],,drop=FALSE])
write.csv(need_gene_exp_13988, "need_gene_exp_13988(13982).csv", row.names = FALSE)

