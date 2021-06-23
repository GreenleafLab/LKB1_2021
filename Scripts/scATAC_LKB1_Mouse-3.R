library(ArchR)

devtools::install_github("GreenleafLab/ArchR", 
  ref = "dev-jg",
  auth_token = "6e052aa859b315361bf5afed59a6061fef4c0eaf", 
  repos = BiocManager::repositories(),
  force =TRUE
)

wd <- "/Volumes/JG_SSD_2/Data/Analysis/LKB1_Mouse2"
dir.create(wd)
setwd(wd)


#Get Input Files
validBC <- getValidBarcodes(
	csvFiles = list.files("/Volumes/JG_SSD_2/Analysis/LKB1_Mouse", pattern = ".csv", full.names = TRUE),
	sampleNames = gsub(".singlecell.csv", "", list.files("/Volumes/JG_SSD_2/Analysis/LKB1_Mouse", pattern = ".csv", full.names = FALSE))
)
inputFiles <- getInputFiles("/Volumes/JG_SSD_2/Analysis/LKB1_Mouse/")

#Set Defaults
addArchRGenome("mm10")
addArchRThreads(10)

#Create Arrows
ArrowFiles <- createArrowFiles(
	inputFiles = inputFiles,
	validBarcodes = validBC
)

#Identify Doublets
doubletScores <- addDoubletScores(list.files(pattern = ".arrow"))

#Create ArchR Project
proj <- ArchRProject(list.files(pattern = ".arrow"))

#Filter Doublets
proj <- filterDoublets(proj)

#Filter Cells
proj <- proj[proj$TSSEnrichment >= 6]
proj <- proj[proj$nFrags <= 50000]

#QC Plots
p1 <- plotGroups(proj, colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.35, pal = paletteDiscrete(proj$Sample, set="grove")) + 
	scale_color_manual(values = paletteDiscrete(proj$Sample, set = "grove")) + scale_fill_manual(values = paletteDiscrete(proj$Sample, set = "grove"))
p2 <- plotGroups(proj, colorBy = "cellColData", name = "nFrags", plotAs = "violin", alpha = 0.35, pal = paletteDiscrete(proj$Sample, set="grove")) + 
	scale_color_manual(values = paletteDiscrete(proj$Sample, set = "grove")) + scale_fill_manual(values = paletteDiscrete(proj$Sample, set = "grove"))
plotPDF(p1, p2, name = "Plot-QC-Stats.pdf", width = 4, height = 4, ArchRProj = proj)

#Reduce Dims
proj <- addIterativeLSI(proj,
	iterations = 2,
	clusterParams = list(
		resolution = c(0.1), 
        sampleCells = 10000, 
        n.start = 10
    ),
    varFeatures = 25000,
    excludeChr = c("chrX", "chrY"),
    force = TRUE
)

#Add Clusters
proj <- addClusters(proj, resolution = 0.4, force = TRUE, reducedDims = "IterativeLSI")

#Add UMAP
proj <- addUMAP(proj, force = TRUE, reducedDims = "IterativeLSI", minDist = 0.6, name = "UMAP4")

#Add TSNE
proj <- addTSNE(proj, force = TRUE, reducedDims = "IterativeLSI")

#Add Genotype Cluster Info
proj$Genotype <- stringr::str_split(proj$Sample, pattern = "\\_", simplify = TRUE)[,2]
proj$GenotypeCluster <- paste0(proj$Clusters, ".", proj$Genotype)

#Rename Clusters
labels <- paste0("Cluster", 1:12)
names(labels) <- paste0("Cluster", c(7, 6, 8, 9, 11, 5, 4, 10, 3, 1, 2, 12))
proj$Clusters2 <- mapLabels(proj$Clusters, labels, names(labels))

#Plot UMAP
p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", size = 0.5, discreteSet = "grove")
p2 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", size = 0.5, discreteSet = "bear")
p3 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", size = 0.5, discreteSet = "calm")
p4 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", size = 0.5, discreteSet = "kelly")
p5 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", size = 0.5, discreteSet = "paired")
p6 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", size = 0.5, discreteSet = "circus")
plotPDF(p1, p2, p3, p4, p5, p6, name = "Plot-UMAP", width = 4, height = 4, ArchRProj = proj)



p2 <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters2", size = 0.5, pal = palLKB1, labelMeans = FALSE, embedding = "UMAP2")
p3 <- plotEmbedding(proj, colorBy = "cellColData", name = "Genotype", pal = c("KPL" = "firebrick3", "KP" = "dodgerblue3"), size = 0.5)
p4 <- plotEmbedding(proj, colorBy = "cellColData", name = "log10(nFrags)")
p5 <- plotEmbedding(proj, colorBy = "cellColData", name = "TSSEnrichment")
p6 <- plotEmbedding(proj, colorBy = "cellColData", name = "DoubletEnrichment")
plotPDF(p1, p2, p3, p4, p5, p6, name = "Plot-UMAP", width = 4, height = 4, ArchRProj = proj)

#Impute Weights
proj <- addImputeWeights(proj)


#######################################################################################
# Analysis
#######################################################################################

#Plot Genes
markerGenes <- c("Sox9", "Sox17", "Sox2", "Nr4a1", "Nr4a2", "Nr4a3", "Il33", "Sik1", "Nkx2-1", "Hmga2", "Sox21", "Sox18")
p1 <- plotEmbedding(
	proj, 
	bins = 150, 
	continuousSet = "comet",
	colorBy = "GeneScoreMatrix", 
	name = markerGenes, 
	imputeWeights = getImputeWeights(proj),
	quantCut = c(0, 1)
)
plotPDF(p1, name = "Plot-UMAP-Markers-Comet", width = 4, height = 4, ArchRProj = proj)

markerGenes <- c("Sox9", "Sox17", "Sox2", "Nr4a1", "Nr4a2", "Nr4a3", "Il33", "Sik1", "Nkx2-1", "Hmga2", "Sox21", "Sox18")
p1 <- plotEmbedding(
	proj, 
	bins = 150, 
	continuousSet = "horizonExtra",
	colorBy = "GeneScoreMatrix", 
	name = markerGenes, 
	imputeWeights = getImputeWeights(proj)
)
plotPDF(p1, name = "Plot-UMAP-Markers-Horizon", width = 4, height = 4, ArchRProj = proj)

#######################################################################################
# Save Results
#######################################################################################
#saveRDS(proj, "Save-Project.rds")
#save.image("Save-ArchR-Feb21-807.Rdata")


#######################################################################################
# Call Peaks
#######################################################################################
proj <- addGroupCoverages(proj, groupBy = "Clusters2")
proj <- addReproduciblePeakSet(proj, groupBy = "Clusters2", pathToMacs2 = "/Users/jeffreygranja/Library/Python/2.7/bin/macs2")
proj <- addPeakMatrix(proj)


#######################################################################################
# Add Peak Annotation
#######################################################################################
proj <- addMotifAnnotations(proj, name = "Motif1", version = 1, force = TRUE)


#######################################################################################
# Add Motif Deviations
#######################################################################################
proj <- addDeviationsMatrix(proj, peakAnnotation = "Motif1", matrixName = "Motif1Matrix", force = TRUE)


#######################################################################################
# Motif Deviations
#######################################################################################

motifs <- grep("z:",getFeatures(proj, select = "Sox17|Nkx2|Runx|Tead|Smarcc|Jun|Fos", useMatrix = "Motif1Matrix"),value=TRUE)
motifs <- sort(motifs)

markerGenes <- c("Sox9", "Sox17", "Sox2", "Nr4a1", "Nr4a2", "Nr4a3", "Il33", "Sik1", "Nkx2-1", "Hmga2", "Sox21", "Sox18")
p1 <- plotEmbedding(
	proj, 
	bins = 150, 
	continuousSet = "solarExtra",
	colorBy = "Motif1Matrix", 
	name = motifs, 
	imputeWeights = getImputeWeights(proj),
	quantCut = c(0, 1)	
)
plotPDF(p1, name = "Plot-UMAP-MotifMarkers", width = 4, height = 4, ArchRProj = proj)


#######################################################################################
# Motif Footprint
#######################################################################################

seNone <- plotFootprints(
	proj, 
	inputSE = seNone,
	pal = palLKB1,
	groupBy = "Clusters2",
	positions = getPositions(proj, "Motif1")[c("Sox17_843", "Onecut2_1054", "Nkx21_446")],
	width = 3,
	height = 5
)


pal <- paletteDiscrete(gtools::mixedsort(unique(proj$Clusters2)))
pal[1] <- pal["Cluster12"]
pal["Cluster12"] <- "#3D3D3D"
palLKB1 <- pal

seDivide <- plotFootprints(
	proj, 
	inputSE = seNone,
	normMethod = "divide",
	groupBy = "Clusters2",
	positions = getPositions(proj, "Motif1")[c("Sox17_843", "Onecut2_1054", "Nkx21_446")]
)

seSubtract <- plotFootprints(
  proj, 
  inputSE = seNone,
  pal = palLKB1,
  normMethod = "subtract",
  groupBy = "Clusters2",
  positions = getPositions(proj, "Motif1")[c("Sox17_843", "Onecut2_1054", "Nkx21_446")]
)


#######################################################################################
# Heatmap Time
#######################################################################################


p1 <- plotFragmentSizes(proj) + scale_color_manual(values = paletteDiscrete(proj$Sample, set = "grove"))
p2 <- plotTSSEnrichment(proj) + scale_color_manual(values = paletteDiscrete(proj$Sample, set = "grove"))
plotPDF(p1, p2, name = "QC-Plot", width = 5, height = 4)


cm <- ArchR:::.confusionMatrix(i = paste0(proj$Clusters2), j = paste0(proj$Genotype))
cm <- cm / rowSums(cm)
colnames(cm)

library(pheatmap)

cm2 <- ArchR:::.binarySort(t(as.matrix(cm[gtools::mixedsort(rownames(cm)),])), clusterCols = FALSE, cutOff = 0.7)

p <- pheatmap(t(cm2[[1]]), border_color = "black", color = paletteContinuous("whiteBlue"), cluster_cols = FALSE, cluster_rows = FALSE)

plotPDF(p, name = "ConfusionHeatmap", width = 5, height = 7)


markersPeaks <- markerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "Clusters2")
markersGenes <- markerFeatures(ArchRProj = proj, useMatrix = "GeneScoreMatrix", groupBy = "Clusters2")


p <- markerPlot(markersGenes, "Cluster12", plotAs = "MA")
plotPDF(p, name = "Cluster12-Markers-Genes", width = 5, height = 7)

p <- markerPlot(markersPeaks, "Cluster12", plotAs = "MA")
plotPDF(p, name = "Cluster12-Markers-Peaks", width = 5, height = 7)

motifsUp <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif1",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,"Cluster12"])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

#Ggplot
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

plotPDF(ggUp, name = "Cluster12-Markers-Motifs-Enriched", width = 5, height = 5)

motifsDo <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif1",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )

df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,"Cluster12"])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

#Ggplot
ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

plotPDF(ggDo, name = "Cluster12-Markers-Motifs-DisEnriched", width = 5, height = 5)




p <- markerPlot(markersGenes, "Cluster12", cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1" , plotAs = "volcano") +
	scale_color_manual(values = 
		c("Not-Differential" = "lightgrey", 
			"Down-Regulated" = "darkorchid", 
			"Up-Regulated" = "forestgreen"
	)) + #scale_x_continuous(limits = c(-5, 5), breaks = seq(-8,8, 2)) +
	ggrepel::geom_text_repel(data = df2, aes(Log2FC, -log10(FDR), label = name), color = "black")

plotPDF(p, name = "Cluster12-Markers-Genes", width = 5, height = 7)
markerDF <- getMarkers(markersGenes, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1")[["Cluster12"]]

df2 <- data.frame(markerDF[markerDF$name %in% c("Sox17", "Nkx2-1"), ])

head(data.frame(markerDF), 25)

heatmapPeaks2 <- markerHeatmap(
  seMarker = markersPeaks[,colnames(cm2[[1]])], 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  clusterCols = FALSE,
  transpose = TRUE
)
plotPDF(heatmapPeaks2, name = "MarkerHeatmap", width = 11, height = 7)


#Identify Motif Enrichments
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif1",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

heatmapEM <- enrichHeatmap(enrichMotifs[,colnames(cm2[[1]])], n = 10, transpose = TRUE)


genes <- getFeatures(proj)


df <- lapply(seq_along(genes), function(x){
  if(x%%1000==0){
    print(x)
  }
  ma <- grep(paste0(tolower(genes[x])), tolower(rownames(enrichMotifs)),ignore.case=TRUE)
  if(length(ma) > 0){
    DataFrame(Genes = genes[x], Motifs = rownames(enrichMotifs)[ma])
  }else{
    NULL
  }
}) %>% Reduce("rbind", .)

gMat <- .getMatrixValues(ArchRProj = proj, name = df[,1], matrixName = "GeneScoreMatrix")

gMat2 <- gMat[which(matrixStats::rowVars(as.matrix(gMat))!=0),]

imputeWeights <- proj@imputeWeights$Weights[colnames(gMat2), colnames(gMat2)]
gMat2 <- (imputeWeights %*% t(gMat2))


gMat3 <- ArchR:::.groupMeans(t(gMat2), groups = getCellColData(proj)[rownames(gMat2),"Clusters2"], sparse = TRUE)

df2 <- df[df[,1] %in% rownames(gMat3),]
m2 <- as.matrix(assay(enrichMotifs[,colnames(gMat3)])[df2[,2],])
m2[m2 < 0] <- 0
df2$cor <- ArchR:::rowCorCpp(X = as.matrix(gMat3[df2[,1],]), Y = m2, idxX = seq_len(nrow(df2)), idxY = seq_len(nrow(df2)))

df3 <- df2[which(df2$cor > 0.5),]

df4 <- df3[which(df3$Genes==stringr::str_split(df3$Motifs,pattern="\\_",simplify=TRUE)),]

heatmapEM <- enrichHeatmap(enrichMotifs[df4[,2],colnames(cm2[[1]])], n = 7, transpose = TRUE, clusterCols = FALSE, rastr = FALSE)

plotPDF(heatmapEM, name = "MotifHeatmap", width = 11, height = 7)


heatmapEM2 <- enrichHeatmap(enrichMotifs[df4[,2],colnames(cm2[[1]])], n = 7, transpose = TRUE, clusterCols = FALSE, returnMatrix = TRUE)
motfis <- stringr::str_split(colnames(heatmapEM2), pattern = " ", simplify=TRUE)[,1]

genesPlot <- df4[match(motfis, df4[,2]),1]

plotM <- t(.rowZscores(gMat3[genesPlot, rownames(heatmapEM2
  )]))
plotM2 <- list(plotM)

ht <- .ArchRHeatmap(
  mat = plotM2[[1]] * 4 / 3,
  scale = FALSE,
  limits = c(-3,3),
  color = paletteContinuous("blueYellow"), 
  clusterCols = FALSE, 
  clusterRows = FALSE,
  labelRows = TRUE,
  labelCols = TRUE,
  useRaster=FALSE,
  draw = FALSE,
  borderColor="black",
  name = ""
  )

plotPDF(ht, name = "GeneHeatmap", width = 11, height = 7)


########################################
# Summarized Experiment
########################################


seATAC <- readRDS("/Volumes/JG_SSD_1/Data/LKB1/KP_KPL_Tumors/output/ATAC_Analysis/Main/ArchRSE_Log2Norm.rds")
assay(seATAC)


regions <- proj@reducedDims[[1]]$LSIFeatures
LSI <- getReducedDims(proj)
embed <- getEmbedding(proj, "UMAP4", returnDF = FALSE)

gr <- GRanges(regions$seqnames, IRanges((regions$idx -1) * 500, width = 500))

subATAC <- subsetByOverlaps(seATAC, gr, ignore.strand = TRUE)

subATAC <- subATAC[order(rowSums(.getAssay(subATAC, "counts")), decreasing = TRUE), ]
o <- DataFrame(findOverlaps(subATAC, gr, ignore.strand = TRUE))
o <- o[!duplicated(o$subjectHits),]
subATAC <- subATAC[o$queryHits, ]
rownames(subATAC) <- paste0("f", o$subjectHits)

bulkMat <- .safeSubset(
	mat = .getAssay(subATAC, "counts"), 
	subsetRows = paste0("f", seq_len(nrow(proj@reducedDims[[1]]$LSIFeatures)))
)



depthN <- quantile(LSI$colSm, 0.5)
n <- 100

simLSI <- .safelapply(seq_len(ncol(bulkMat)), function(x){
	print(x)
	counts <- bulkMat[, x]
	simMat <- lapply(seq_len(n), function(x){
		as(tabulate(sample(x = rep(seq_along(counts), counts), size = depthN), nbins = length(counts)), "sparseMatrix")
	}) %>% Reduce("cbind", .)
	simLSI <- as.matrix(.projectLSI(simMat, LSI = LSI, verbose = FALSE))
	rownames(simLSI) <- paste0(colnames(bulkMat)[x], seq_len(n))
	simLSI
}, threads = 10) %>% Reduce("rbind", .)

if(is.null(embed$params$scaleDims)){
	scaleDims <- LSI$scaleDims
	if(is.na(scaleDims)){
		scaleDims <- FALSE
	}
}else{
	scaleDims <- embed$params$scaleDims
}
if(scaleDims){
	simLSI <- .scaleDims(simLSI)
} 
if(!is.null(embed$params$dimsToUse)){
	simLSI <- simLSI[, embed$params$dimsToUse]
}



proj <- addUMAP(proj, force = TRUE, reducedDims = "IterativeLSI", minDist = 0.6, name = "UMAP4")
embed <- getEmbedding(proj, "UMAP4", returnDF = FALSE)



model <- .loadUWOT(embed$params$uwotModel, embed$params$nc)
simUMAP <- uwot::umap_transform(X = rbind(LSI, simLSI), model = model, verbose = TRUE, n_threads = 8)

pdf("tes4.pdf")
ggPoint(embed[[1]][,1], embed[[1]][,2], paste0(proj$Genotype))
#ggPoint(simUMAP[,1], simUMAP[,2])
ggPoint(simUMAP[,1], simUMAP[,2], 
	c(rep("ref", nrow(embed[[1]])), stringr::str_split(rownames(simLSI), pattern = "\\_", simplify = TRUE)[,2]),
	size = 0.2, rastr = TRUE)

ggPoint(simUMAP[,1], simUMAP[,2], 
	c(rep("ref", nrow(embed[[1]])), stringr::str_split(rownames(simLSI), pattern = "\\_", simplify = TRUE)[,3]),
	size = 0.2, rastr = TRUE)

ggPoint(simUMAP[,1], simUMAP[,2], 
	c(rep("ref", nrow(embed[[1]])), stringr::str_split(rownames(simLSI), pattern = "_Rep", simplify = TRUE)[,1]),
	size = 0.2, rastr = TRUE)

dev.off()









ggPoint(c(embed[[1]][,1], simUMAP[,1]), c(embed[[1]][,2], simUMAP[,2]), 
	c(rep("ref", nrow(embed[[1]])), stringr::str_split(rownames(simLSI), pattern = "\\_", simplify = TRUE)[,2]),
	size = 0.2, rastr = TRUE)
ggPoint(c(embed[[1]][,1], simUMAP[,1]), c(embed[[1]][,2], simUMAP[,2]), 
	c(rep("ref", nrow(embed[[1]])), stringr::str_split(rownames(simLSI), pattern = "\\_", simplify = TRUE)[,3]),
	size = 0.2, rastr = TRUE)




simUMAP <- uwot::umap_transform(X = rbind(LSI[[1]], simLSI), model = umap, verbose = TRUE, n_threads = 8)
.saveUWOT(umap, "test.tar")

getOutputDirectory()

mat <- getReducedDims(proj)
umap <- uwot::umap(mat, metric = "euclidean", ret_nn = TRUE, ret_model = TRUE, min_dist = 0.4, n_neighbors = 40, verbose = TRUE)

.saveUWOT(umap, "test3.tar")
save_uwot(umap, file2)
model <- .loadUWOT("test3.tar", 30)
model <- load_uwot("test.tar")

file <- "test.tar"
file2 <- file.path(normalizePath(dirname(file)), basename(file))



proj <- addUMAP(proj, force = TRUE, reducedDims = "IterativeLSI", minDist = 0.1, name = "UMAP4")
embed <- getEmbedding(proj, "UMAP4", returnDF = FALSE)



p <- ArchRBrowserTrack(proj, region = GRanges("chr1", IRanges(4347353,4597354)), groupBy = "Clusters2")
plotPDF(p, name = "test")






cm <- ArchR:::.confusionMatrix(proj$Sample, proj$Clusters2)
cm1 <- as.matrix(cm/ Matrix::rowSums(cm))
cm2 <- as.matrix(t(t(cm) / Matrix::colSums(cm)))

d1 <- reshape2::melt(cm1)
d1$Var1 <- factor(paste0(d1$Var1), levels = gtools::mixedsort(unique(paste0(d1$Var1))))
d1$Var2 <- factor(paste0(d1$Var2), levels = gtools::mixedsort(unique(paste0(d1$Var2))))

pdf("Barplot1.pdf", width = 5, height = 5)
ggplot(d1, aes(x = Var1, y = value, fill = Var2)) + 
	geom_bar(stat = "identity") +
	scale_fill_manual(values = palLKB1) +
	theme_ArchR(xText90 = TRUE) +
	scale_y_continuous(limits = c(0,1), breaks = seq(0,10)/5, expand = c(0, 0))
dev.off()

d2 <- reshape2::melt(t(cm2))
d2$Var1 <- factor(paste0(d2$Var1), levels = gtools::mixedsort(unique(paste0(d2$Var1))))
d2$Var2 <- factor(paste0(d2$Var2), levels = gtools::mixedsort(unique(paste0(d2$Var2))))

pdf("Barplot2.pdf", width = 5, height = 5)
ggplot(d2, aes(x = Var1, y = value, fill = Var2)) + 
	geom_bar(stat = "identity") +
	scale_fill_manual(values=paletteDiscrete(proj$Sample, set = "grove")) +
	theme_ArchR(xText90 = TRUE) +
	scale_y_continuous(limits = c(0,1), breaks = seq(0,10)/5, expand = c(0, 0))
dev.off()



seATAC <- exportGroupSE(proj, useMatrix = "PeakMatrix", groupBy = "Clusters2")

log2cpm <- log2(t(t(assay(seATAC)) / colSums(assay(seATAC))) * 10^6 + 1)

cm <- cor(log2cpm)
cm[cm < 0.6] <- 0.6
library(pheatmap)
ht <- pheatmap(cm, color = paletteContinuous(), border_color = "black")
plotPDF(ht, name = "CorMap")

save.image("Save-LKB1-Mouse-Feb23.Rdata")


varPeaks <- head(order(matrixStats::rowVars(log2cpm), decreasing = TRUE), 10000)
pcaVar <- prcomp(t(log2cpm[varPeaks, ]))

pVar <- round(100*pcaVar$sdev^2/sum(pcaVar$sdev^2),3)
pVar <- paste0(paste0("PC", seq_along(pVar))," ",pVar,"%")
plotPCA <- data.frame(pcaVar$x)
plotPCA$Group <- colnames(log2cpm)


pdf("PCA-10k-Variable-Peaks.pdf", width = 4, height = 4, useDingbats = FALSE)
p <- ggplot(plotPCA, aes(PC1,PC2,color=Group)) + 
	geom_point(size=1.5) + theme_ArchR() + 
	xlab(pVar[1]) + ylab(pVar[2]) +
	scale_color_manual(values = palLKB1)
ArchR:::.fixPlotSize(p, plotWidth = 4, plotHeight = 4, height=1)
dev.off()



.getMatrixValues(ArrowFiles = getArrowFiles(proj))




markersTF <- markerFeatures(
	ArchRProj = proj, 
	useMatrix = "Motif1Matrix", 
	groupBy = "Clusters2",
	useGroups = c("Cluster12"),
	bgdGroups = paste0("Cluster", 6:11),
	testMethod = "ttest",
	useSeqnames = "z"
)


p <- markerPlot(markersTF, , name = "Cluster12", plotAs = "VolcanoDiff", cutOff = "abs(MeanDiff) > 3 & FDR <= 0.01")
plotPDF(p, name = "Test-Diff-CV")

df <- getMarkers(markersTF, cutOff = "abs(MeanDiff) > 3 & FDR <= 0.01")[[1]]


df <- getMarkers(markersTF, cutOff = "FDR <= 10")[[1]]

devB <- readRDS("/Volumes/JG_SSD_1/Data/LKB1/KP_KPL_Tumors/output/ATAC_Analysis/Main/CisBP-chromVAR.rds")

idx1 <- grep("KPL_Met", colnames(devB))
idx2 <- grep("KPL_Primary", colnames(devB))

diffCV <- rowMeans(assays(devB[,idx1])[["z"]]) - rowMeans(assays(devB[,idx2])[["z"]])

df$diffBulk <- diffCV[df$name]

#plotDF <- readRDS("Save_diffDev_MM_HT.rds")[[1]]

plotDF <- df
plotDF$diff1 <- plotDF$diffBulk
plotDF$diff2 <- plotDF$MeanDiff
plotDF <- plotDF[which(rowSums(is.na(plotDF[,c("diff1","diff2")]))==0),]
plotDF$dist <- plotDF$diff1^2 + plotDF$diff2^2
plotDF$Family <- "NA"
plotDF$Sig <- FALSE
minCut <- 1
plotDF[abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Sig"] <- TRUE
plotDF[abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Family"] <- "SigOther"
plotDF$Common <- plotDF$name

for(i in seq_along(colorList)){
  idx <- unlist(lapply(colorList[[i]], function(x)grep(x,plotDF$Common, value = TRUE, ignore.case = TRUE)))
  plotDF$Family[which(plotDF$Common %in% idx)] <- names(colorList)[i]
}

plotDF <- plotDF[which(plotDF$Family != "FOXOther"), ]
plotDF[which(!plotDF$Sig), "Family"] <- "NA"

idxCenter <- which(abs(plotDF$diff1) < 3 & abs(plotDF$diff2) < 3)
plotCenter <- plotDF[idxCenter, ]
plotCenter <- plotCenter[order(plotCenter$dist, decreasing = TRUE), ]
plotCenter <- plotCenter[sample(seq_len(nrow(plotCenter)), nrow(plotCenter)), ]
plotDF2 <- plotDF[-idxCenter,]
plotDF2 <- plotDF2[order(plotDF2$dist, decreasing = FALSE), ]
plotDF <- rbind(plotDF2 , plotCenter)
plotDF <- rbind(plotDF[plotDF$Family %in% c("NA", "SigOther"),], plotDF[plotDF$Family %ni% c("NA", "SigOther"),])

max <- 40
plotDF$diff1[plotDF$diff1 > max] <- max
plotDF$diff1[plotDF$diff1 < -max] <- -max
plotDF$diff2[plotDF$diff2 > max] <- max
plotDF$diff2[plotDF$diff2 < -max] <- -max

pdf("Plot-diffDev_KPL_BSC.pdf", width = 5, height = 5, useDingbats = FALSE)

xlim <- extendrange(c(-max,max), f = 0)
ylim <- extendrange(c(-max,max), f = 0)

xlim <- extendrange(plotDF$diff2, f = 0)
ylim <- extendrange(plotDF$diff1, f = 0)

p <- ggplot(data = data.frame(plotDF), aes(diff2, diff1, color = Family, fill = Family)) +
  geom_point(size = 3, shape = 21, stroke = 0.2, color = "black") + 
  scale_fill_manual(values = pal) + 
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_vline(xintercept = 0, lty = "dashed") +
  xlab("ChromVAR scATAC KPL Met - Tumor") +
  ylab("ChromVAR bulk KPL Met - Tumor") + 
  coord_equal(ratio = diff(xlim)/diff(ylim), xlim = xlim, ylim = ylim) + 
  theme_ArchR() + 
  ggtitle(max)

ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5)

dev.off()


save.image("Save-Feb24.Rdata")




colorList <- list(
  "NKX" = c("NKX"),
  "NFI" = c("NFI"),
  "MEF2" = c("MEF2"),
  "TEAD" = c("TEAD"),
  "RUNX" = c("RUNX", "RUNT"),
  "SOX" = c("SOX", "SRY"),
  #"MESP" = c("MESP"),
  "ZEB" = c("ZEB", "SNAI", "MESP"),
  "FOXOther" = c("FOX"),
  "FOXA" = c("FOXA1", "FOXA2"), 
  "NR4A" = c("NR4A"),
  "BZIP" = c("JUN","FOS","BATF","FOSL", "ATF3", "JDP", "BACH"),
  "ONECUT" = c("ONECUT"),
  "ETS" = c("ETS", "ERG")
)

pal <- paletteDiscrete(values = names(colorList))
pal["SigOther"] <- "gray47"
pal["NA"] <- "lightgrey"
pal["NKX"] <- "#00ADEE"
pal["NFI"] <- "#C06CAB"
pal["MEF2"] <- "#FFDE17"
pal["TEAD"] <- "#F15A29"
pal["RUNX"] <- "#8B5E3C"
pal["SOX"] <- "#00A551"
pal["ZEB"] <- "#92278F"
pal["FOXOther"] <- "bisque"
pal["FOXA"] <- "#272E6A"
pal["NR4A"] <- "#EC008C"
pal["BZIP"] <- "#D51F26"
pal["ONECUT"] <- "slategray1"
pal["ETS"] <- "slateblue1"






gMat <- .getMatrixValues(ArchRProj = proj, 
	name = c("Nkx2-1", "Sox17", "Sox9", "Sox2", "Sox18", "Hmga2", "Sik1", "Nr4a1"), 
	matrixName = "GeneScoreMatrix", log2Norm = FALSE)
gMat2 <- gMat[which(matrixStats::rowVars(as.matrix(gMat))!=0),]

imputeWeights <- proj@imputeWeights$Weights[colnames(gMat2), colnames(gMat2)]
gMat2 <- (imputeWeights %*% t(gMat2))

gMat3 <- log2(gMat2 + 1)

pdf("Test.pdf")
ggGroup(proj$Clusters2, gMat3[, 1], title = "NKX2-1")
ggGroup(proj$Clusters2, gMat3[, 2], title = "Sox17")
dev.off()


for(i in seq_len(ncol(gMat3))){

	top <- quantile(gMat3[, i], 0.99)
	groupHi <- top * 2/3
	groupMed <- top * 1/3
	idCells <- rep("Low", nrow(gMat3))
	idCells[gMat3[, i] > groupMed] <- "Med"
	idCells[gMat3[, i] > groupHi] <- "High"

	cm <- .confusionMatrix(proj$Sample, idCells)
	cm <- cm / rowSums(cm)


	d1 <- reshape2::melt(as.matrix(cm))
	d1$Var1 <- factor(paste0(d1$Var1), levels = gtools::mixedsort(unique(paste0(d1$Var1))))
	d1$Var2 <- factor(paste0(d1$Var2), levels = rev(c("Low", "Med", "High")))

	pal <- c("Low" = "#FF7900", "Med" = "#7600D3", "High" = "#2FAEED")

	pdf(paste0("Barplot-",colnames(gMat3)[i],".pdf"), width = 5, height = 5)
	print(ggGroup(proj$Sample, gMat3[, i], title = colnames(gMat3)[i]))
	print(ggGroup(proj$Clusters2, gMat3[, i], title = colnames(gMat3)[i]))
	p <- ggplot(d1, aes(x = Var1, y = value, fill = Var2)) + 
		geom_bar(stat = "identity") +
		scale_fill_manual(values = pal) +
		theme_ArchR(xText90 = TRUE) + ggtitle(colnames(gMat3)[i])
	print(p)
	dev.off()

}






top <- quantile(gMat3[, 1], 0.99)
groupHi <- top * 2/3
groupMed <- top * 1/3
idCells <- rep("Low", nrow(gMat3))
idCells[gMat3[, 1] > groupMed] <- "Med"
idCells[gMat3[, 1] > groupHi] <- "High"

cm <- .confusionMatrix(proj$Sample, idCells)
cm <- cm / rowSums(cm)


d1 <- reshape2::melt(as.matrix(cm))
d1$Var1 <- factor(paste0(d1$Var1), levels = gtools::mixedsort(unique(paste0(d1$Var1))))
d1$Var2 <- factor(paste0(d1$Var2), levels = rev(c("Low", "Med", "High")))

pal <- c("Low" = "#FF7900", "Med" = "#7600D3", "High" = "#2FAEED")

pdf("Barplot-NKX21.pdf", width = 5, height = 5)
ggplot(d1, aes(x = Var1, y = value, fill = Var2)) + 
	geom_bar(stat = "identity") +
	scale_fill_manual(values = pal) +
	theme_ArchR(xText90 = TRUE)
dev.off()




top <- quantile(gMat3[, 2], 0.99)
groupHi <- top * 2/3
groupMed <- top * 1/3
idCells <- rep("Low", nrow(gMat3))
idCells[gMat3[, 2] > groupMed] <- "Med"
idCells[gMat3[, 2] > groupHi] <- "High"

cm <- .confusionMatrix(proj$Sample, idCells)
cm <- cm / rowSums(cm)


d1 <- reshape2::melt(as.matrix(cm))
d1$Var1 <- factor(paste0(d1$Var1), levels = gtools::mixedsort(unique(paste0(d1$Var1))))
d1$Var2 <- factor(paste0(d1$Var2), levels = rev(c("Low", "Med", "High")))

pal <- c("Low" = "#FF7900", "Med" = "#7600D3", "High" = "#2FAEED")

pdf("Barplot-Sox17.pdf", width = 5, height = 5)
ggplot(d1, aes(x = Var1, y = value, fill = Var2)) + 
	geom_bar(stat = "identity") +
	scale_fill_manual(values = pal) +
	theme_ArchR(xText90 = TRUE)
dev.off()









