#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
#biocLite("CAGEr")
library(CAGEr)
library(BSgenome.Hsapiens.UCSC.hg19)

## workspace start
pathsToInputFiles <- "/media/nikos/data1/mnt/raid0/georgaki/projects/miRNA_genes/DATASETS/FANTOM/FANTOM5/Human/H1_ESC/nucleus/ENCFF000TXQ.bam"
#pathsToInputFiles <- list.files(inputDir, full.names = TRUE)
myCAGEset <- new("CAGEset",genomeName="BSgenome.Hsapiens.UCSC.hg19",inputFiles=pathsToInputFiles,inputFilesType = "bam",sampleLabels = c("ENCFF000TXQ"))
#getCTSS(myCAGEset) 
getCTSS(myCAGEset, sequencingQualityThreshold = 10, mappingQualityThreshold = 5)

ctss <- CTSStagCount(myCAGEset) #(get count of read for each CTSS)
head(ctss)
sampleLabels(myCAGEset)
librarySizes(myCAGEset)
#normalizeTagCount(myCAGEset, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 5*10^4)
normalizeTagCount(myCAGEset, method = "simpleTpm")
exportCTSStoBedGraph(myCAGEset, values = "normalized", oneFile = TRUE)
#save.image("~/biothesis/CAGE_TOOLS_TEST/data/CAGER_OUT/cageR_after_create_myCAGEset_and_getCTSS__normilized_15_6.RData")
#save.image("~/biothesis/CAGE_TOOLS_TEST/data/CAGER_OUT/cageR_after_create_myCAGEset_and_getCTSS__normilized_15_6_1.RData")
#clusterCTSS(object = myCAGEset, threshold = 1, thresholdIsTpm = TRUE,nrPassThreshold = 1, method = "distclu", maxDist = 20,removeSingletons = TRUE, keepSingletonsAbove = 5)
#clusterCTSS(object = myCAGEset, threshold = 1, thresholdIsTpm = TRUE,nrPassThreshold = 1, method = "paraclu", maxDist = 20,removeSingletons = TRUE, keepSingletonsAbove = 5)
#clusterCTSS(object = myCAGEset, threshold = 0, thresholdIsTpm = TRUE,nrPassThreshold = 1, method = "distclu", maxDist = 20,removeSingletons = TRUE, keepSingletonsAbove = 5)
#clusterCTSS(object = myCAGEset, threshold = 0.5, thresholdIsTpm = TRUE,nrPassThreshold = 1, method = "distclu", maxDist = 20,removeSingletons = TRUE, keepSingletonsAbove = 5)
save.image("~/biothesis/CAGE_TOOLS_TEST/data/CAGER_OUT/cageR_after_create_myCAGEset_and_getCTSS_SQ10_MQ5_SIMPLETMP.RData")

clusterCTSS(object = myCAGEset, threshold = 1, nrPassThreshold = 1,thresholdIsTpm = TRUE, method = "paraclu", maxDist = 20,removeSingletons = TRUE, keepSingletonsAbove = 5)


tc <- tagClusters(myCAGEset, sample = "ENCFF000TXQ")
save.image("~/biothesis/CAGE_TOOLS_TEST/data/CAGER_OUT/cageR_after_create_myCAGEset_and_getCTSS_SQ10_MQ5_SIMPLETMP_AFTER_CLUSTERING.RData")
tc
tc1<-tc[c("chr","start","end","cluster","tpm","strand")]
rownames(tc1) <- NULL
head(tc1)
#write.table(tc1, file="~/biothesis/CAGE_TOOLS_TEST/data/CAGER_OUT/cage_TC_threshold_0.5.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(tc1, file="~/biothesis/CAGE_TOOLS_TEST/data/CAGER_OUT/cage_TC_threshold_SQ10_MQ5.bed", quote=F, sep="\t", row.names=F, col.names=F)

## workspace end

