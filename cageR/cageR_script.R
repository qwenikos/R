#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome")
#biocLite("BSgenome.Drerio.UCSC.danRer7")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
#biocLite("CAGEr")
library(CAGEr)
library(BSgenome.Drerio.UCSC.danRer7)
library(BSgenome.Hsapiens.UCSC.hg19)

myCAGEsetcopy<-myCAGEset
mycageCoors<-myCAGEsetcopy@CTSScoordinates

inputDir <- system.file("extdata", package = "CAGEr")
pathsToInputFiles <- list.files(inputDir, full.names = TRUE)
basename(pathsToInputFiles)
myCAGEset <- new("CAGEset", genomeName = "BSgenome.Drerio.UCSC.danRer7",inputFiles = pathsToInputFiles, inputFilesType = "ctss",sampleLabels = c("zf_30p_dome", "zf_high", "zf_prim6_rep1", "zf_prim6_rep2","zf_unfertilized_egg"))
myCAGEset
TSS.df <- read.table(system.file("extdata/Zf.unfertilized.egg.chr17.ctss", package = "CAGEr"))
colnames(TSS.df) <- c("chr", "pos", "strand", "zf_unfertilized_egg")
TSS.df$chr <- as.character(TSS.df$chr)
TSS.df$pos <- as.integer(TSS.df$pos)
TSS.df$strand <- as.character(TSS.df$strand)
TSS.df$zf_unfertilized_egg <- as.integer(TSS.df$zf_unfertilized_egg)
head(TSS.df)
myCAGEset.coerced <- as(TSS.df, "CAGEset")
myCAGEset.coerced

getCTSS(myCAGEset)
ctss <- CTSStagCount(myCAGEset)
head(ctss)
sampleLabels(myCAGEset)
corr.m <- plotCorrelation(myCAGEset, samples = "all", method = "pearson")
mergeSamples(myCAGEset, mergeIndex = c(3,2,4,4,1),mergedSampleLabels = c("zf_unfertilized_egg","zf_high", "zf_30p_dome", "zf_prim6"))
librarySizes(myCAGEset)
y = -1 * alpha * x + beta
plotReverseCumulatives(myCAGEset, fitInRange = c(5, 1000), onePlot = TRUE)
normalizeTagCount(myCAGEset, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 5*10^4)
exportCTSStoBedGraph(myCAGEset, values = "normalized",format = "bedGraph", oneFile = TRUE)
exportCTSStoBedGraph(myCAGEset, values = "normalized", format = "BigWig")
clusterCTSS(object = myCAGEset, threshold = 1, thresholdIsTpm = TRUE,nrPassThreshold = 1, method = "distclu", maxDist = 20,removeSingletons = TRUE, keepSingletonsAbove = 5)
tc <- tagClusters(myCAGEset, sample = "zf_unfertilized_egg")
head(tc)
cumulativeCTSSdistribution(myCAGEset, clusters = "tagClusters")
quantilePositions(myCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
tc <- tagClusters(myCAGEset, sample = "zf_unfertilized_egg",returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
head(tc)
exportToBed(object = myCAGEset, what = "tagClusters",qLow = 0.1, qUp = 0.9, oneFile = TRUE)
plotInterquantileWidth(myCAGEset, clusters = "tagClusters",tpmThreshold = 3, qLow = 0.1, qUp = 0.9)
aggregateTagClusters(myCAGEset, tpmThreshold = 5,qLow = 0.1, qUp = 0.9, maxDist = 100)
consensusCl <- consensusClusters(myCAGEset)
head(consensusCl)
consensusCl <- consensusClusters(myCAGEset, sample = "zf_unfertilized_egg",returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
getExpressionProfiles(myCAGEset, what = "consensusClusters", tpmThreshold = 10,nrPassThreshold = 1, method = "som", xDim = 4, yDim = 2)
plotExpressionProfiles(myCAGEset, what = "consensusClusters")
class3_1 <- extractExpressionClass(myCAGEset,what = "consensusClusters", which = "3_1")
head(class3_1)
exportToBed(myCAGEset, what = "consensusClusters", colorByExpressionProfile = TRUE)
cumulativeCTSSdistribution(myCAGEset, clusters = "consensusClusters")
scoreShift(myCAGEset, groupX = "zf_unfertilized_egg", groupY = "zf_prim6",testKS = TRUE, useTpmKS = FALSE)
shifting.promoters <- getShiftingPromoters(myCAGEset,tpmThreshold = 5, scoreThreshold = 0.6,fdrThreshold = 0.01)
head(shifting.promoters)


data(FANTOM5humanSamples)
head(FANTOM5humanSamples)
astrocyteSamples <- FANTOM5humanSamples[grep("Astrocyte", FANTOM5humanSamples[,"description"]),]
astrocyteSamples
astrocyteSamples[,"sample"]
astrocyteCAGEset <- importPublicData(source = "FANTOM5", dataset = "human",sample = astrocyteSamples[1:3,"sample"])
mixedCAGEset <- importPublicData(source = "FANTOM3and4",
                                   dataset = c("FANTOMtissueCAGEmouse", "FANTOMtissueCAGEmouse",
                                                 "FANTOMtimecourseCAGEmouse"), group = c("liver", "liver",
                                                                                           "liver_under_constant_darkness"),
                                   sample = c("cloned_mouse", "control_mouse", "4_hr"))
mixedCAGEset
sessionInfo()


library(BSgenome.Hsapiens.UCSC.hg19)
pathsToInputFiles <- "/media/nikos/data1/mnt/raid0/georgaki/projects/miRNA_genes/DATASETS/FANTOM/FANTOM5/Human/H1_ESC/nucleus/ENCFF000TXQ.bam"
#pathsToInputFiles <- list.files(inputDir, full.names = TRUE)
basename(pathsToInputFiles)

myCAGEset <- new("CAGEset",genomeName="BSgenome.Hsapiens.UCSC.hg19",inputFiles=pathsToInputFiles,inputFilesType = "bam",sampleLabels = c("ENCFF000TXQ"))
getCTSS(myCAGEset)       #read the file in structure 
ifctss <- CTSStagCount(myCAGEset)
myCAGEset
ctss <- CTSStagCount(myCAGEset) #(get count of read for each CTSS)
ctss
librarySizes(myCAGEset) # get the total number of cage sample
plotReverseCumulatives(myCAGEset, fitInRange = c(5, 1000), onePlot = TRUE)
myCAGEsetcopy<-myCAGEset
##normalization tag per milion
normalizeTagCount(myCAGEsetcopy, method = "simpleTpm")
normalizeTagCount(myCAGEsetcopy, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 5*10^4)

clusterCTSS(object = myCAGEsetcopy, threshold = 1, thresholdIsTpm = TRUE,
              nrPassThreshold = 1, method = "distclu", maxDist = 20,
              removeSingletons = TRUE, keepSingletonsAbove = 5)

exportCTSStoBedGraph(myCAGEsetcopy, values = "normalized", format = "bedGraph", oneFile = TRUE)
exportToBed(myCAGEsetcopy, "tagClusters", qLow = NULL, qUp = NULL, colorByExpressionProfile = FALSE, oneFile = TRUE)
exportToBed(myCAGEsetcopy, "tagClusters", qLow = NULL, qUp = NULL, colorByExpressionProfile = FALSE, oneFile = TRUE)
#load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#exportCTSStoBedGraph(exampleCAGEset, values = "normalized")
head(myCAGEset@tagCountMatrix)
